// ss_gmp_bench.c
// Fast Solovay–Strassen primality test with GMP on 512-bit integers.
// Prints primes found; times each test using Apple Silicon-friendly cycle counter.
//
// Build (Apple Silicon / Homebrew):
//   clang -O3 -flto -fomit-frame-pointer -mcpu=apple-m1 -std=c11 \
//     -I/opt/homebrew/include -L/opt/homebrew/lib \
//     -o ss_gmp_bench ss_gmp_bench.c -lgmp
//
// Examples:
//   ./ss_gmp_bench
//   ./ss_gmp_bench --count 20000 --bits 512 --rounds 10
//   ./ss_gmp_bench --use-gmp --count 50000 --bits 512 --no-print-primes
//
// Notes:
// - Solovay–Strassen uses the Jacobi symbol: a^( (n-1)/2 ) ≡ (a/n) (mod n) for odd prime n.
// - In practice, MR or GMP’s mpz_probab_prime_p are stronger/faster; SS here is for parity.

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(__APPLE__)
  #include <mach/mach_time.h>
#endif

// ===================== Apple Silicon "cycle" counter =====================

static inline unsigned long long read_cycles(void) {
#if defined(__aarch64__) && defined(__APPLE__) && defined(__has_builtin)
  #if __has_builtin(__builtin_readcyclecounter)
    return __builtin_readcyclecounter(); // very low overhead, close to cycles
  #else
    return mach_continuous_time();       // ticks (monotonic), not literal cycles
  #endif
#else
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (unsigned long long)ts.tv_sec * 1000000000ull + (unsigned long long)ts.tv_nsec;
#endif
}

// ===================== Small-prime sieve (reject early) =====================

static const unsigned SMALL_PRIMES[] = {
  3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
  101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,
  211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,
  337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,
  461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,
  601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,
  739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,
  881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997
};
static const size_t N_SMALL = sizeof(SMALL_PRIMES)/sizeof(SMALL_PRIMES[0]);

static inline int small_sieve_composite(const mpz_t n) {
    // n < 2 => composite (not prime), even (≠2) => composite
    if (mpz_cmp_ui(n, 2) < 0) return 1;
    if (mpz_even_p(n)) return mpz_cmp_ui(n, 2) != 0;

    for (size_t i = 0; i < N_SMALL; ++i) {
        unsigned p = SMALL_PRIMES[i];
        int cmp = mpz_cmp_ui(n, p);
        if (cmp == 0) return 0;               // equal to small prime → prime
        if (mpz_tdiv_ui(n, p) == 0) return 1; // divisible → composite
    }
    return 0; // inconclusive
}

// ===================== Solovay–Strassen core =====================
//
// For odd n > 2:
//   pick base a ∈ [2, n-2]
//   j = Jacobi(a, n) ∈ {-1,0,1}
//   x = a^((n-1)/2) mod n
//   if j == 0 → composite
//   else compare x with j mod n: x == 1  (if j= 1) or x == n-1 (if j=-1)
// Repeat for R rounds; if all pass → probable prime.

typedef struct {
    mpz_t t, a, x, nm1, n_minus_3;
} ss_ctx;

static inline void ss_ctx_init(ss_ctx* c) {
    mpz_init(c->t);
    mpz_init(c->a);
    mpz_init(c->x);
    mpz_init(c->nm1);
    mpz_init(c->n_minus_3);
}
static inline void ss_ctx_clear(ss_ctx* c) {
    mpz_clear(c->t);
    mpz_clear(c->a);
    mpz_clear(c->x);
    mpz_clear(c->nm1);
    mpz_clear(c->n_minus_3);
}

// One Solovay–Strassen round with base 'a' (in ctx->a). Returns 1 pass, 0 composite.
static inline int ss_round(const mpz_t n, ss_ctx* c) {
    // j = Jacobi(a, n)
    int j = mpz_jacobi(c->a, n);
    if (j == 0) return 0; // a shares factor with n → composite

    // t = (n-1)/2
    mpz_sub_ui(c->t, n, 1);
    mpz_fdiv_q_2exp(c->t, c->t, 1);

    // x = a^t mod n
    mpz_powm(c->x, c->a, c->t, n);

    // Convert j ∈ {-1,1} into modulus representative in {1, n-1}
    if (j == 1) {
        // expect x == 1 (mod n)
        return mpz_cmp_ui(c->x, 1) == 0;
    } else { // j == -1
        // expect x == n-1
        return mpz_cmp(c->x, c->nm1) == 0;
    }
}

// Returns 1 for probable prime, 0 for composite.
static inline int is_probable_prime_ss(const mpz_t n, int rounds, gmp_randstate_t rng, ss_ctx* c) {
    // Quick rejects & exact small primes
    if (small_sieve_composite(n)) return 0;

    // SS requires odd n > 2; we’ve filtered evens and n<2.
    // Prepare constants once per n
    mpz_sub_ui(c->nm1, n, 1);
    mpz_sub_ui(c->n_minus_3, n, 3);

    // A few cheap fixed bases first (optional)
    static const unsigned FIXED_BASES[] = { 2u, 3u, 5u, 7u, 11u };
    int k = (rounds < (int)(sizeof(FIXED_BASES)/sizeof(FIXED_BASES[0])) ?
             rounds : (int)(sizeof(FIXED_BASES)/sizeof(FIXED_BASES[0])));
    for (int i = 0; i < k; ++i) {
        if (mpz_cmp_ui(n, FIXED_BASES[i]) <= 0) continue;
        mpz_set_ui(c->a, FIXED_BASES[i]);
        if (!ss_round(n, c)) return 0;
    }
    if (rounds <= k) return 1;

    // Random bases in [2, n-2]
    for (int i = k; i < rounds; ++i) {
        mpz_urandomm(c->a, rng, c->n_minus_3);
        mpz_add_ui(c->a, c->a, 2);
        if (!ss_round(n, c)) return 0;
    }
    return 1;
}

// ===================== Candidate generation, CLI, printing =====================

static inline void rand_odd_bigint(mpz_t x, gmp_randstate_t rng, unsigned bits) {
    mpz_urandomb(x, rng, bits);
    mpz_setbit(x, bits-1); // exact bit-length
    mpz_setbit(x, 0);      // odd
}

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [--count N] [--bits B] [--rounds R] [--use-gmp] [--no-print-primes]\n"
        "  --count N        number of random odd candidates (default 10000)\n"
        "  --bits B         bit-length of candidates (default 512)\n"
        "  --rounds R       SS rounds (default 12)\n"
        "  --use-gmp        use GMP's mpz_probab_prime_p instead (faster/stronger baseline)\n"
        "  --no-print-primes  do not print primes found (faster)\n",
        prog);
}

int main(int argc, char** argv) {
    int count = 10000;
    unsigned bits = 512;
    int rounds = 12;
    int use_gmp = 0;
    int print_primes = 1;

    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--help")) { usage(argv[0]); return 0; }
        else if (!strcmp(argv[i], "--count") && i+1 < argc) { count = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--bits") && i+1 < argc) { bits = (unsigned)atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--rounds") && i+1 < argc) { rounds = atoi(argv[++i]); }
        else if (!strcmp(argv[i], "--use-gmp")) { use_gmp = 1; }
        else if (!strcmp(argv[i], "--no-print-primes")) { print_primes = 0; }
        else { usage(argv[0]); return 1; }
    }
    if (count <= 0 || bits < 16 || rounds < 1) { usage(argv[0]); return 1; }

    // RNG
    gmp_randstate_t rng;
    gmp_randinit_default(rng);
    gmp_randseed_ui(rng, (unsigned long)time(NULL));

    mpz_t n; mpz_init(n);
    ss_ctx ctx; ss_ctx_init(&ctx);

    unsigned long long* t = (unsigned long long*)malloc((size_t)count * sizeof(*t));
    if (!t) { fprintf(stderr, "OOM\n"); return 1; }

    int primes = 0, prime_no = 0;
    unsigned long long sum = 0, minv = ~0ull, maxv = 0;
    char* hexbuf = NULL;

    for (int i = 0; i < count; ++i) {
        rand_odd_bigint(n, rng, bits);

        unsigned long long t0 = read_cycles();
        int is_pp;
        if (use_gmp) {
            // GMP’s tuned probabilistic test (often faster/stronger than SS)
            is_pp = mpz_probab_prime_p(n, rounds) > 0;
        } else {
            is_pp = is_probable_prime_ss(n, rounds, rng, &ctx);
        }
        unsigned long long t1 = read_cycles();

        unsigned long long dt = t1 - t0;
        t[i] = dt;
        sum += dt;
        if (dt < minv) minv = dt;
        if (dt > maxv) maxv = dt;

        if (is_pp) {
            primes++;
            if (print_primes) {
                free(hexbuf);
                hexbuf = mpz_get_str(NULL, 16, n);
                printf("Prime #%d (candidate index %d)\n  hex: %s\n", ++prime_no, i, hexbuf);
            }
        }
    }

    const double ln2 = 0.693147180559945309417232121458176568;
    int expected = (int)((double)count / ((double)bits * ln2) + 0.5);
    double avg = (double)sum / (double)count;

    printf("\nTested %d random %u-bit odd integers.\n", count, bits);
    printf("Probable primes found: %d (expected ~ %d)\n", primes, expected);
    printf("Per-test counter: avg = %.2f, min = %llu, max = %llu\n",
           avg, (unsigned long long)minv, (unsigned long long)maxv);

    free(hexbuf);
    free(t);
    ss_ctx_clear(&ctx);
    mpz_clear(n);
    gmp_randclear(rng);
    return 0;
}