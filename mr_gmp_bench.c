// mr_gmp_bench.c
// Fast Miller–Rabin (MR) primality test with GMP on 512-bit ints.
// Prints any primes found; times each test with a low-overhead counter on Apple Silicon.
//
// Build (Apple Silicon / Homebrew):
//   clang -O3 -flto -mcpu=apple-m1 -std=c11 \
//     -I/opt/homebrew/include -L/opt/homebrew/lib \
//     -o mr_gmp_bench mr_gmp_bench.c -lgmp
//
// Examples:
//   ./mr_gmp_bench
//   ./mr_gmp_bench --count 20000 --bits 512 --rounds 8
//   ./mr_gmp_bench --use-gmp --count 20000 --bits 512
//
// Notes:
// - "cycles" uses __builtin_readcyclecounter() when Clang exposes it; else mach_continuous_time() ticks.
// - Printing primes is on by default; use --no-print-primes to suppress for cleaner timing.

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(__APPLE__)
  #include <mach/mach_time.h>
#endif

// ===================== Cycle counter (Apple Silicon friendly) =====================

static inline unsigned long long read_cycles(void) {
#if defined(__aarch64__) && defined(__APPLE__) && defined(__has_builtin)
  #if __has_builtin(__builtin_readcyclecounter)
    return __builtin_readcyclecounter(); // very low overhead, close to cycles
  #else
    return mach_continuous_time();       // ticks, monotonic
  #endif
#else
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (unsigned long long)ts.tv_sec * 1000000000ull + (unsigned long long)ts.tv_nsec;
#endif
}

// ===================== Small-prime sieve (bigger to kill composites early) =====================

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
    if (mpz_cmp_ui(n, 2) < 0) return 1;               // n < 2 → composite (not prime)
    if (mpz_even_p(n)) return mpz_cmp_ui(n, 2) != 0;  // even and != 2 → composite

    for (size_t i = 0; i < N_SMALL; ++i) {
        unsigned p = SMALL_PRIMES[i];
        int cmp = mpz_cmp_ui(n, p);
        if (cmp == 0) return 0;                       // exactly a small prime → prime
        if (mpz_tdiv_ui(n, p) == 0) return 1;         // divisible by small prime → composite
    }
    return 0; // inconclusive
}

// ===================== MR core (inline helpers, reused temporaries) =====================

typedef struct {
    mpz_t d, x, nm1, a, n_minus_3;
    unsigned s;
} mr_ctx;

static inline void mr_ctx_init(mr_ctx* c) {
    mpz_init(c->d);
    mpz_init(c->x);
    mpz_init(c->nm1);
    mpz_init(c->a);
    mpz_init(c->n_minus_3);
    c->s = 0;
}
static inline void mr_ctx_clear(mr_ctx* c) {
    mpz_clear(c->d);
    mpz_clear(c->x);
    mpz_clear(c->nm1);
    mpz_clear(c->a);
    mpz_clear(c->n_minus_3);
}

static inline void split_n_minus_1(const mpz_t n, mr_ctx* c) {
    mpz_sub_ui(c->d, n, 1);
    c->s = 0;
    while (mpz_even_p(c->d)) {
        mpz_fdiv_q_2exp(c->d, c->d, 1); // d >>= 1
        c->s++;
    }
    mpz_sub_ui(c->nm1, n, 1);
    mpz_sub_ui(c->n_minus_3, n, 3);
}

// strong test to base c->a; returns 1 pass, 0 composite.
static inline int mr_strong_test_base(const mpz_t n, mr_ctx* c) {
    // x = a^d mod n
    mpz_powm(c->x, c->a, c->d, n);
    if (mpz_cmp_ui(c->x, 1) == 0 || mpz_cmp(c->x, c->nm1) == 0) return 1;

    // repeat s-1 times: x = x^2 mod n; if x == n-1 pass
    for (unsigned r = 1; r < c->s; ++r) {
        // Squaring explicitly is fine; GMP may special-case powm_ui(...,2) too.
        mpz_mul(c->x, c->x, c->x);
        mpz_mod(c->x, c->x, n);
        if (mpz_cmp(c->x, c->nm1) == 0) return 1;
    }
    return 0;
}

// Probable-prime test with k fixed small bases then (rounds-k) random bases.
// Returns 1 probable prime, 0 composite.
static inline int is_probable_prime_mr(const mpz_t n, int rounds, gmp_randstate_t rng, mr_ctx* c) {
    if (small_sieve_composite(n)) return 0;

    split_n_minus_1(n, c);

    // Cheap deterministic bases first (reduce RNG & variance).
    static const unsigned FIXED_BASES[] = {2u,3u,5u,7u,11u,13u,17u,19u};
    const int K = (rounds < (int)(sizeof(FIXED_BASES)/sizeof(FIXED_BASES[0])) ?
                   rounds : (int)(sizeof(FIXED_BASES)/sizeof(FIXED_BASES[0])));
    for (int i = 0; i < K; ++i) {
        if (mpz_cmp_ui(n, FIXED_BASES[i]) <= 0) continue; // if n <= base, skip
        mpz_set_ui(c->a, FIXED_BASES[i]);
        if (!mr_strong_test_base(n, c)) return 0;
    }
    if (rounds <= K) return 1;

    // Random bases in [2, n-2]
    for (int i = K; i < rounds; ++i) {
        mpz_urandomm(c->a, rng, c->n_minus_3);
        mpz_add_ui(c->a, c->a, 2);
        if (!mr_strong_test_base(n, c)) return 0;
    }
    return 1;
}

// ===================== Random candidates, CLI, printing =====================

static inline void rand_odd_bigint(mpz_t x, gmp_randstate_t rng, unsigned bits) {
    mpz_urandomb(x, rng, bits);
    mpz_setbit(x, bits-1); // ensure exact bitlength
    mpz_setbit(x, 0);      // odd
}

static void usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [--count N] [--bits B] [--rounds R] [--use-gmp] [--no-print-primes]\n"
        "  --count N        number of random odd candidates (default 10000)\n"
        "  --bits B         bit-length of candidates (default 512)\n"
        "  --rounds R       MR rounds (default 12; used when not --use-gmp)\n"
        "  --use-gmp        use GMP's mpz_probab_prime_p instead of custom MR (fast)\n"
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
    mr_ctx ctx; mr_ctx_init(&ctx);

    unsigned long long* t = (unsigned long long*)malloc((size_t)count*sizeof(*t));
    if (!t) { fprintf(stderr, "OOM\n"); return 1; }

    int primes = 0;
    int prime_no = 0;
    unsigned long long sum = 0, minv = ~0ull, maxv = 0;

    char* hexbuf = NULL;

    for (int i = 0; i < count; ++i) {
        rand_odd_bigint(n, rng, bits);

        unsigned long long t0 = read_cycles();
        int is_pp;
        if (use_gmp) {
            // reps=rounds maps to MR repetitions inside GMP; 12 is strong for 512-bit.
            is_pp = mpz_probab_prime_p(n, rounds) > 0;
        } else {
            is_pp = is_probable_prime_mr(n, rounds, rng, &ctx);
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

    double avg = (double)sum / (double)count;
    // expected ≈ count / ln(2^bits) = count / (bits * ln 2)
    const double ln2 = 0.693147180559945309417232121458176568;
    int expected = (int)( (double)count / ( (double)bits * ln2 ) + 0.5 );

    printf("\nTested %d random %u-bit odd integers.\n", count, bits);
    printf("Probable primes found: %d (expected ~ %d)\n", primes, expected);
    printf("Per-test counter: avg = %.2f, min = %llu, max = %llu\n",
           avg, (unsigned long long)minv, (unsigned long long)maxv);

    free(hexbuf);
    free(t);
    mr_ctx_clear(&ctx);
    mpz_clear(n);
    gmp_randclear(rng);
    return 0;
}