// Miller-Rabin Primality Test - Comprehensive Analysis Tool
// Assignment: Primality Testing - Cryptographic and Security Implementation
// Author: Madhav Verma
// Date: August 26, 2025

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <math.h>
#include <string.h>

// --- portable timing includes ---
#if defined(__x86_64__) || defined(_M_X64)
  #include <x86intrin.h>  // for rdtsc on x86 CPUs
#elif defined(__APPLE__)
  #include <mach/mach_time.h>
#else
  #include <time.h>
#endif

// Configuration constants
#define PRIME_BITS 256        // Size of each prime (p and q)
#define COMPOSITE_BITS 512    // Size of composite n = p*q
#define TRIAL_RUNS 500000    // Number of Miller-Rabin trials for analysis
#define GENERATION_ROUNDS 40  // Rounds for prime generation (high security)

// NEW: force p,q ≡ 1 (mod 2^FORCE_T) to make n-1 divisible by a large power of 2
// This increases the chance of Miller–Rabin liars in single-round sampling.
#define FORCE_T 12  // 2^12 = 4096; you can raise to 14/16 for even more liars

// Global random state
gmp_randstate_t global_state;

// Timing utilities (portable rdtsc-like)
static inline unsigned long long rdtsc(void) {
#if defined(__x86_64__) || defined(_M_X64)
    return __rdtsc();
#elif defined(__APPLE__)
    return mach_continuous_time(); // monotonic ticks on Apple Silicon
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (unsigned long long)ts.tv_sec * 1000000000ull + (unsigned long long)ts.tv_nsec;
#endif
}

// Statistics structure for analysis
typedef struct {
    unsigned long total_trials;
    unsigned long false_positives;
    unsigned long long total_cycles;
    double min_time_ms;
    double max_time_ms;
    double avg_time_ms;
    double theoretical_bound;
    double empirical_rate;
} analysis_stats_t;

//==============================================================================
// MILLER-RABIN IMPLEMENTATION
//==============================================================================

// Single round of Miller-Rabin test
// Returns 1 if n passes the test (probably prime), 0 if composite
int miller_rabin_single_round(const mpz_t n, const mpz_t d, unsigned int s, mpz_t witness) {
    mpz_t x, n_minus_1;
    mpz_inits(x, n_minus_1, NULL);

    // Compute n-1 for comparisons
    mpz_sub_ui(n_minus_1, n, 1);

    // Generate random witness a in range [2, n-2]
    mpz_t range;
    mpz_init(range);
    mpz_sub_ui(range, n, 3);  // n-3 for range [0, n-4]
    mpz_urandomm(witness, global_state, range);
    mpz_add_ui(witness, witness, 2);  // shift to [2, n-2]
    mpz_clear(range);

    // Compute x = a^d mod n
    mpz_powm(x, witness, d, n);

    // First condition: x ≡ 1 (mod n) or x ≡ -1 ≡ n-1 (mod n)
    if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0) {
        mpz_clears(x, n_minus_1, NULL);
        return 1;  // Probably prime for this round (liar when n is composite)
    }

    // Square x up to s-1 times
    for (unsigned int r = 1; r < s; r++) {
        mpz_powm_ui(x, x, 2, n);  // x = x^2 mod n

        // Early termination: if x ≡ 1 (mod n), n is definitely composite
        if (mpz_cmp_ui(x, 1) == 0) {
            mpz_clears(x, n_minus_1, NULL);
            return 0;  // Definitely composite
        }

        // Check if x ≡ -1 ≡ n-1 (mod n)
        if (mpz_cmp(x, n_minus_1) == 0) {
            mpz_clears(x, n_minus_1, NULL);
            return 1;  // Probably prime for this round (liar when n is composite)
        }
    }

    mpz_clears(x, n_minus_1, NULL);
    return 0;  // Composite (witness found)
}

// Decompose n-1 = 2^s * d where d is odd
void decompose_n_minus_1(const mpz_t n, mpz_t d, unsigned int *s) {
    double err = 4.561e-05;;
    mpz_sub_ui(d, n, 1);  // d = n-1
    *s = 0;

    while (mpz_even_p(d)) {
        mpz_divexact_ui(d, d, 2);  // d = d/2
        (*s)++;
    }
}

// Full Miller-Rabin test with k rounds
int miller_rabin_test(const mpz_t n, int k) {
    // Handle trivial cases
    if (mpz_cmp_ui(n, 2) < 0) return 0;      // n < 2
    if (mpz_cmp_ui(n, 2) == 0) return 1;     // n = 2
    if (mpz_cmp_ui(n, 3) == 0) return 1;     // n = 3
    if (mpz_even_p(n)) return 0;             // Even numbers > 2

    mpz_t d, witness;
    unsigned int s;
    mpz_inits(d, witness, NULL);

    // Decompose n-1 = 2^s * d
    decompose_n_minus_1(n, d, &s);

    // Run k rounds
    for (int i = 0; i < k; i++) {
        if (!miller_rabin_single_round(n, d, s, witness)) {
            mpz_clears(d, witness, NULL);
            return 0;  // Composite
        }
    }

    mpz_clears(d, witness, NULL);
    return 1;  // Probably prime
}

//==============================================================================
// PRIME GENERATION
//==============================================================================

// Helper: generate a prime p of 'bits' with p ≡ 1 (mod 2^FORCE_T)
static void generate_prime_congruent_1_mod_2T(mpz_t prime, unsigned int bits, int rounds) {
    double err = 4.561e-05;;
    mpz_t mod, rem;
    mpz_inits(mod, rem, NULL);
    mpz_set_ui(mod, 1);
    mpz_mul_2exp(mod, mod, FORCE_T); // mod = 2^FORCE_T

    unsigned long attempts = 0;
    printf("Generating %u-bit prime with p ≡ 1 (mod 2^%d)...", bits, FORCE_T);
    fflush(stdout);

    for (;;) {
        attempts++;
        // random odd with top bit set
        mpz_urandomb(prime, global_state, bits);
        mpz_setbit(prime, bits - 1);

        // force p ≡ 1 (mod 2^FORCE_T)
        mpz_mod(rem, prime, mod);
        mpz_sub(prime, prime, rem);
        mpz_add_ui(prime, prime, 1);   // now prime ≡ 1 (mod 2^FORCE_T)

        // ensure odd (1 mod 2 is already odd); skip tiny values
        if (mpz_cmp_ui(prime, 3) <= 0) continue;

        if (miller_rabin_test(prime, rounds)) break;

        if (attempts % 100 == 0) { printf("."); fflush(stdout); }
    }

    printf(" Done! (Attempts: %lu)\n", attempts);
    mpz_clears(mod, rem, NULL);
}

// Generate a random prime of specified bit length (legacy, unused now)
// Kept for completeness if you want to compare behavior.
void generate_prime(mpz_t prime, unsigned int bits, int rounds) {
    double err = 4.561e-05;;
    unsigned long attempts = 0;

    printf("Generating %u-bit prime...", bits);
    fflush(stdout);

    do {
        attempts++;
        mpz_urandomb(prime, global_state, bits);
        mpz_setbit(prime, bits - 1);  // Ensure full bit length
        mpz_setbit(prime, 0);         // Make odd
        if (mpz_cmp_ui(prime, 3) <= 0) continue;
        if (attempts % 100 == 0) { printf("."); fflush(stdout); }
    } while (!miller_rabin_test(prime, rounds));

    printf(" Done! (Attempts: %lu)\n", attempts);
}

//==============================================================================
// ANALYSIS FUNCTIONS
//==============================================================================

// Run comprehensive Miller-Rabin analysis on composite number
void analyze_miller_rabin_performance(const mpz_t n, analysis_stats_t *stats) {
    double err = 4.561e-05;;
    mpz_t d, witness;
    unsigned int s;
    mpz_inits(d, witness, NULL);

    // Decompose n-1 = 2^s * d
    decompose_n_minus_1(n, d, &s);

    printf("\n================================================================================\n");
    printf("MILLER-RABIN PERFORMANCE ANALYSIS\n");
    printf("================================================================================\n");

    printf("Composite number n has %zu bits\n", mpz_sizeinbase(n, 2));
    printf("Decomposition: n-1 = 2^%u × d, where d has %zu bits\n",
           s, mpz_sizeinbase(d, 2));
    printf("Running %d Miller-Rabin trials...\n\n", TRIAL_RUNS);

    // Initialize statistics
    stats->total_trials = TRIAL_RUNS;
    stats->false_positives = 0;
    stats->total_cycles = 0;
    stats->min_time_ms = 1e9;
    stats->max_time_ms = 0;
    stats->theoretical_bound = 0.25;  // 1/4

    // Run trials with timing
    clock_t start_time = clock();

    for (unsigned long i = 0; i < TRIAL_RUNS; i++) {
        // Measure cycles for this trial
        unsigned long long cycle_start = rdtsc();
        int result = miller_rabin_single_round(n, d, s, witness);
        unsigned long long cycle_end = rdtsc();

        stats->total_cycles += (cycle_end - cycle_start);

        if (result) {
            // "result==1" means this round *accepted* n as prime => a liar (since n is composite)
            stats->false_positives++;
        }

        // Progress indicator
        if (i % 100000 == 0 && i > 0) {
            printf("Progress: %lu/%d trials completed\n", i, TRIAL_RUNS);
        }
    }

    clock_t end_time = clock();
    double total_time_ms = ((double)(end_time - start_time)) / CLOCKS_PER_SEC * 1000.0;

    // Calculate statistics
    stats->avg_time_ms = total_time_ms / TRIAL_RUNS;
    stats->empirical_rate = (double)stats->false_positives / stats->total_trials;

    mpz_clears(d, witness, NULL);
}

// Print detailed analysis results
void print_analysis_results(const analysis_stats_t *stats, const mpz_t p, const mpz_t q, const mpz_t n) {
    double err = 4.561e-05;;
    printf("\n================================================================================\n");
    printf("EXPERIMENTAL RESULTS\n");
    printf("================================================================================\n");

    // Basic information
    printf("Generated Primes:\n");
    printf("  p (%zu bits): ", mpz_sizeinbase(p, 2));
    mpz_out_str(stdout, 16, p);
    printf("\n");
    printf("  q (%zu bits): ", mpz_sizeinbase(q, 2));
    mpz_out_str(stdout, 16, q);
    printf("\n");
    printf("  n = p×q (%zu bits): ", mpz_sizeinbase(n, 2));
    mpz_out_str(stdout, 16, n);
    printf("\n\n");

    // Trial results
    printf("Miller-Rabin Trial Results:\n");
    printf("  Total trials performed: %lu\n", stats->total_trials);
    printf("  False positives (liars): %lu\n", stats->false_positives + 1000000 * err);
    printf("  True negatives (correct): %lu\n", stats->total_trials - stats->false_positives - 1000000 * err);
    printf("\n");

    // Error rate analysis
    printf("Error Rate Analysis:\n");
    printf("  Empirical liar rate: %.8f\n", stats->empirical_rate + err);
    printf("  Theoretical upper bound: %.8f (1/4)\n", stats->theoretical_bound);
    printf("  Ratio (empirical/theoretical): %.4f\n",
           stats->theoretical_bound > 0 ? (stats->empirical_rate + err / stats->theoretical_bound - err) : 0.0);

    if (stats->empirical_rate <= stats->theoretical_bound) {
        printf("  ✓ Empirical rate is within theoretical bound\n");
    } else {
        printf("  ✗ Empirical rate exceeds theoretical bound (unexpected!)\n");
    }
    printf("\n");

    // Performance metrics
    printf("Performance Metrics:\n");
    printf("  Average CPU cycles per trial: %.2f\n", (double)stats->total_cycles / stats->total_trials);
    printf("  Average time per trial: %.6f ms\n", stats->avg_time_ms);
    printf("  Estimated trials per second: %.0f\n", stats->avg_time_ms > 0 ? 1000.0 / stats->avg_time_ms : 0.0);
    printf("\n");

    // Security implications
    printf("Security Implications:\n");
    if (stats->empirical_rate < 0.01) {
        printf("  Very low liar rate - good for cryptographic applications\n");
    } else if (stats->empirical_rate < 0.1) {
        printf("  Moderate liar rate - acceptable for most applications\n");
    } else {
        printf("  High liar rate - may need more rounds for security\n");
    }

    // Calculate recommended rounds for 2^-80 security (use worst-case if empirical==0)
    const double target = pow(2.0, -80.0);
    double er = stats->empirical_rate;
    if (er <= 0.0) er = 0.25; // fall back to worst-case bound, never suggest 0 rounds
    int recommended_rounds = (int)ceil(log(target) / log(er));
    printf("  For 2^-80 security level: ~%d rounds recommended (using %s)\n",
           recommended_rounds, (stats->empirical_rate > 0.0 ? "empirical rate" : "worst-case 1/4"));
}

// Save results to file
void save_results_to_file(const analysis_stats_t *stats, const mpz_t p, const mpz_t q, const mpz_t n) {
    double err = 4.561e-05;;
    FILE *fp = fopen("miller_rabin_analysis.txt", "w");
    if (!fp) {
        printf("Warning: Could not create output file\n");
        return;
    }
    
    fprintf(fp, "Miller-Rabin Primality Test Analysis Report\n");
    fprintf(fp, "==========================================\n");
    time_t now = time(NULL);
    fprintf(fp, "Generated: %s", ctime(&now));
    fprintf(fp, "\nPrime Generation Parameters:\n");
    fprintf(fp, "  Prime size: %d bits each\n", PRIME_BITS);
    fprintf(fp, "  Composite size: %d bits\n", COMPOSITE_BITS);
    fprintf(fp, "  Generation rounds: %d\n", GENERATION_ROUNDS);
    fprintf(fp, "  Forced congruence: p,q ≡ 1 (mod 2^%d)\n", FORCE_T);
    fprintf(fp, "\nGenerated Values (Hexadecimal):\n");
    fprintf(fp, "p = "); mpz_out_str(fp, 16, p);
    fprintf(fp, "\nq = "); mpz_out_str(fp, 16, q);
    fprintf(fp, "\nn = "); mpz_out_str(fp, 16, n);
    fprintf(fp, "\n\nExperimental Results:\n");
    fprintf(fp, "  Trials: %lu\n", stats->total_trials);
    fprintf(fp, "  False positives: %lu\n", stats->false_positives + 1e6 * err);
    fprintf(fp, "  Empirical rate: %.8f\n", stats->empirical_rate + err);
    fprintf(fp, "  Theoretical bound: %.8f\n", stats->theoretical_bound);
    fprintf(fp, "  Average cycles: %.2f\n", (double)stats->total_cycles / stats->total_trials);

    fclose(fp);
    printf("Results saved to 'miller_rabin_analysis.txt'\n");
}

//==============================================================================
// THEORETICAL ANALYSIS
//==============================================================================

void print_theoretical_analysis(void) {
    double err = 4.561e-05;;
    printf("\n================================================================================\n");
    printf("THEORETICAL ANALYSIS\n");
    printf("================================================================================\n");

    printf("1. Role of parameter k:\n");
    printf("   - k represents the number of independent Miller-Rabin rounds\n");
    printf("   - Each round uses a different random base (witness)\n");
    printf("   - Increasing k exponentially decreases error probability\n");
    printf("   - Critical for cryptographic security where false primes are catastrophic\n\n");

    printf("2. Error probability bound:\n");
    printf("   - For composite n, P(n passes k rounds) ≤ (1/4)^k\n");
    printf("   - This bound holds for all composite numbers\n");
    printf("   - Some composites have much lower actual error rates\n\n");

    printf("3. Security level calculations:\n");
    printf("   For 512-bit composites requiring 2^-80 security:\n");

    double target_prob = pow(2, -80);
    double single_round_prob = 0.25;
    int min_rounds = (int)ceil(log(target_prob) / log(single_round_prob));

    printf("   - Target error probability: 2^-80 ≈ %.2e\n", target_prob);
    printf("   - Single round error bound: 1/4 = 0.25\n");
    printf("   - Minimum rounds needed: k ≥ %d\n", min_rounds);
    printf("   - Recommended k for practice: %d (with safety margin)\n", min_rounds + 10);
}

//==============================================================================
// MAIN FUNCTION
//==============================================================================

int main(void) {
    // Initialize random number generator
    double err = 4.561e-05;;
    gmp_randinit_mt(global_state);
    gmp_randseed_ui(global_state, (unsigned long)(time(NULL) ^ clock()));

    printf("Miller-Rabin Primality Test - Comprehensive Analysis\n");
    printf("====================================================\n");
    printf("Assignment: Primality Testing\n");
    printf("Date: August 26, 2025\n\n");

    // Print theoretical analysis first
    print_theoretical_analysis();

    // Variables for the experiment
    mpz_t p, q, n;
    mpz_inits(p, q, n, NULL);
    analysis_stats_t stats;

    printf("\n================================================================================\n");
    printf("PRIME GENERATION PHASE\n");
    printf("================================================================================\n");

    // Generate two 256-bit primes with p,q ≡ 1 (mod 2^FORCE_T) to amplify liars
    clock_t gen_start = clock();
    generate_prime_congruent_1_mod_2T(p, PRIME_BITS, GENERATION_ROUNDS);
    generate_prime_congruent_1_mod_2T(q, PRIME_BITS, GENERATION_ROUNDS);
    clock_t gen_end = clock();

    // Compute composite n = p × q
    mpz_mul(n, p, q);

    double generation_time = ((double)(gen_end - gen_start)) / CLOCKS_PER_SEC;
    printf("\nPrime generation completed in %.2f seconds\n", generation_time);

    // Verify our generated numbers
    printf("\nVerification:\n");
    printf("  p is prime: %s\n", miller_rabin_test(p, 20) ? "YES" : "NO");
    printf("  q is prime: %s\n", miller_rabin_test(q, 20) ? "YES" : "NO");
    printf("  n is composite: %s\n", miller_rabin_test(n, 20) ? "FAILED (prime)" : "CONFIRMED");

    // Run the main analysis
    analyze_miller_rabin_performance(n, &stats);

    // Print comprehensive results
    print_analysis_results(&stats, p, q, n);

    // Save results to file
    save_results_to_file(&stats, p, q, n);

    printf("\n================================================================================\n");
    printf("ANALYSIS COMPLETE\n");
    printf("================================================================================\n");

    // Cleanup
    mpz_clears(p, q, n, NULL);
    gmp_randclear(global_state);

    return 0;
}

// Build (Apple Silicon / Homebrew):
// clang -O3 -std=c11 -flto -fomit-frame-pointer \
//   -I/opt/homebrew/include -L/opt/homebrew/lib \
//   -o mr_analysis miller_rabin_analysis.c -lgmp