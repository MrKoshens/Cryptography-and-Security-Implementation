// Solovay–Strassen Primality Test - Comprehensive Analysis Tool
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
#define TRIAL_RUNS 1000000    // Number of Solovay–Strassen trials for analysis
#define GENERATION_ROUNDS 40  // Rounds for prime generation (high security)

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
// SOLOVAY–STRASSEN IMPLEMENTATION
//==============================================================================
//
// One SS round for odd n>2 with random base a in [2, n-2]:
//   j = Jacobi(a, n) ∈ {-1,0,1}; if j==0 → composite
//   t = (n-1)/2
//   x = a^t mod n
//   pass iff (j== 1 and x==1) or (j==-1 and x==n-1)
// Returns 1 if n passes the round (probable prime), 0 if composite.
//

// Single round of Solovay–Strassen
// Returns 1 if n passes the test (probably prime), 0 if composite
int solovay_strassen_single_round(const mpz_t n, mpz_t scratch_a, mpz_t scratch_t, mpz_t scratch_x) {
    // choose random a in [2, n-2]
    mpz_t range, nm1;
    mpz_inits(range, nm1, NULL);
    mpz_sub_ui(range, n, 3);     // range size: n-3 → values 0..n-4
    mpz_urandomm(scratch_a, global_state, range);
    mpz_add_ui(scratch_a, scratch_a, 2); // shift to 2..n-2
    mpz_sub_ui(nm1, n, 1);

    // Jacobi(a, n)
    int j = mpz_jacobi(scratch_a, n);
    if (j == 0) { mpz_clears(range, nm1, NULL); return 0; } // composite

    // t = (n-1)/2
    mpz_fdiv_q_2exp(scratch_t, nm1, 1);

    // x = a^t mod n
    mpz_powm(scratch_x, scratch_a, scratch_t, n);

    int pass;
    if (j == 1) {
        pass = (mpz_cmp_ui(scratch_x, 1) == 0);
    } else { // j == -1
        pass = (mpz_cmp(scratch_x, nm1) == 0);
    }

    mpz_clears(range, nm1, NULL);
    return pass;
}

// Full Solovay–Strassen test with k rounds
int solovay_strassen_test(const mpz_t n, int k) {
    // Handle trivial cases
    if (mpz_cmp_ui(n, 2) < 0) return 0;      // n < 2
    if (mpz_cmp_ui(n, 2) == 0) return 1;     // n = 2
    if (mpz_cmp_ui(n, 3) == 0) return 1;     // n = 3
    if (mpz_even_p(n)) return 0;             // even > 2 → composite

    mpz_t a, t, x;
    mpz_inits(a, t, x, NULL);

    for (int i = 0; i < k; i++) {
        if (!solovay_strassen_single_round(n, a, t, x)) {
            mpz_clears(a, t, x, NULL);
            return 0; // composite
        }
    }

    mpz_clears(a, t, x, NULL);
    return 1; // probably prime
}

//==============================================================================
// PRIME GENERATION
//==============================================================================

// Generate a random prime of specified bit length (using SS test)
void generate_prime(mpz_t prime, unsigned int bits, int rounds) {
    unsigned long attempts = 0;

    printf("Generating %u-bit prime...", bits);
    fflush(stdout);

    do {
        attempts++;
        // Generate random odd number with MSB set
        mpz_urandomb(prime, global_state, bits);
        mpz_setbit(prime, bits - 1);  // Ensure full bit length
        mpz_setbit(prime, 0);         // Make odd

        // Skip trivial cases
        if (mpz_cmp_ui(prime, 3) <= 0) continue;

        if (attempts % 100 == 0) {
            printf(".");
            fflush(stdout);
        }

    } while (!solovay_strassen_test(prime, rounds));

    printf(" Done! (Attempts: %lu)\n", attempts);
}

//==============================================================================
// ANALYSIS FUNCTIONS
//==============================================================================

// Run comprehensive Solovay–Strassen analysis on composite number
void analyze_solovay_strassen_performance(const mpz_t n, analysis_stats_t *stats) {
    mpz_t a, t, x, nm1;
    mpz_inits(a, t, x, nm1, NULL);
    mpz_sub_ui(nm1, n, 1);

    printf("\n================================================================================\n");
    printf("SOLOVAY–STRASSEN PERFORMANCE ANALYSIS\n");
    printf("================================================================================\n");

    // Show basic exponent size for SS ((n-1)/2)
    mpz_t half_exp; mpz_init(half_exp);
    mpz_fdiv_q_2exp(half_exp, nm1, 1);
    printf("Composite number n has %zu bits\n", mpz_sizeinbase(n, 2));
    printf("Exponent (n-1)/2 has %zu bits\n", mpz_sizeinbase(half_exp, 2));
    mpz_clear(half_exp);

    printf("Running %d Solovay–Strassen trials...\n\n", TRIAL_RUNS);

    // Initialize statistics
    stats->total_trials = TRIAL_RUNS;
    stats->false_positives = 0;
    stats->total_cycles = 0;
    stats->min_time_ms = 1e9;
    stats->max_time_ms = 0;
    stats->theoretical_bound = 0.5;  // ≤ 1/2 per round for composites

    // Run trials with timing
    clock_t start_time = clock();

    for (unsigned long i = 0; i < TRIAL_RUNS; i++) {
        unsigned long long cycle_start = rdtsc();
        int result = solovay_strassen_single_round(n, a, t, x);
        unsigned long long cycle_end = rdtsc();

        stats->total_cycles += (cycle_end - cycle_start);

        if (result) {
            stats->false_positives++;
        }

        if (i % 100000 == 0 && i > 0) {
            printf("Progress: %lu/%d trials completed\n", i, TRIAL_RUNS);
        }
    }

    clock_t end_time = clock();
    double total_time_ms = ((double)(end_time - start_time)) / CLOCKS_PER_SEC * 1000.0;

    // Calculate statistics
    stats->avg_time_ms = total_time_ms / TRIAL_RUNS;
    stats->empirical_rate = (double)stats->false_positives / stats->total_trials;

    mpz_clears(a, t, x, nm1, NULL);
}

// Print detailed analysis results
void print_analysis_results(const analysis_stats_t *stats, const mpz_t p, const mpz_t q, const mpz_t n) {
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
    printf("Solovay–Strassen Trial Results:\n");
    printf("  Total trials performed: %lu\n", stats->total_trials);
    printf("  False positives (liars): %lu\n", stats->false_positives);
    printf("  True negatives (correct): %lu\n", stats->total_trials - stats->false_positives);
    printf("\n");

    // Error rate analysis
    printf("Error Rate Analysis:\n");
    printf("  Empirical liar rate: %.8f\n", stats->empirical_rate);
    printf("  Theoretical upper bound: %.8f (1/2)\n", stats->theoretical_bound);
    printf("  Ratio (empirical/theoretical): %.4f\n", stats->empirical_rate / stats->theoretical_bound);

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
    printf("  Estimated trials per second: %.0f\n", 1000.0 / stats->avg_time_ms);
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

    // Calculate recommended rounds for 2^-80 security (using empirical rate)
    double security_level = pow(2, -80);
    int recommended_rounds = (int)ceil(log(security_level) / log(stats->empirical_rate));
    printf("  For 2^-80 security level: ~%d rounds recommended\n", recommended_rounds);
}

// Save results to file
void save_results_to_file(const analysis_stats_t *stats, const mpz_t p, const mpz_t q, const mpz_t n) {
    FILE *fp = fopen("solovay_strassen_analysis.txt", "w");
    if (!fp) {
        printf("Warning: Could not create output file\n");
        return;
    }

    fprintf(fp, "Solovay–Strassen Primality Test Analysis Report\n");
    fprintf(fp, "================================================\n");
    fprintf(fp, "Generated: %s\n", ctime(&(time_t){time(NULL)}));
    fprintf(fp, "\nPrime Generation Parameters:\n");
    fprintf(fp, "  Prime size: %d bits each\n", PRIME_BITS);
    fprintf(fp, "  Composite size: %d bits\n", COMPOSITE_BITS);
    fprintf(fp, "  Generation rounds: %d\n", GENERATION_ROUNDS);
    fprintf(fp, "\nGenerated Values (Hexadecimal):\n");
    fprintf(fp, "p = ");
    mpz_out_str(fp, 16, p);
    fprintf(fp, "\nq = ");
    mpz_out_str(fp, 16, q);
    fprintf(fp, "\nn = ");
    mpz_out_str(fp, 16, n);
    fprintf(fp, "\n\nExperimental Results:\n");
    fprintf(fp, "  Trials: %lu\n", stats->total_trials);
    fprintf(fp, "  False positives: %lu\n", stats->false_positives);
    fprintf(fp, "  Empirical rate: %.8f\n", stats->empirical_rate);
    fprintf(fp, "  Theoretical bound: %.8f\n", stats->theoretical_bound);
    fprintf(fp, "  Average cycles: %.2f\n", (double)stats->total_cycles / stats->total_trials);

    fclose(fp);
    printf("Results saved to 'solovay_strassen_analysis.txt'\n");
}

//==============================================================================
// THEORETICAL ANALYSIS
//==============================================================================

void print_theoretical_analysis(void) {
    printf("\n================================================================================\n");
    printf("THEORETICAL ANALYSIS (Solovay–Strassen)\n");
    printf("================================================================================\n");

    printf("1. Role of parameter k:\n");
    printf("   - k is the number of independent Solovay–Strassen rounds (random bases)\n");
    printf("   - Each round checks a^( (n-1)/2 ) ≡ (a/n) (mod n) using the Jacobi symbol\n");
    printf("   - Increasing k exponentially decreases the error probability\n\n");

    printf("2. Error probability bound:\n");
    printf("   - For composite odd n, P(n passes one SS round) ≤ 1/2\n");
    printf("   - Therefore P(pass k rounds) ≤ (1/2)^k\n\n");

    printf("3. Security level calculations:\n");
    printf("   For 512-bit composites requiring 2^-80 security:\n");

    double target_prob = pow(2, -80);
    double single_round_prob = 0.5;
    int min_rounds = (int)ceil(log(target_prob) / log(single_round_prob));

    printf("   - Target error probability: 2^-80 ≈ %.2e\n", target_prob);
    printf("   - Single round error bound: 1/2 = 0.5\n");
    printf("   - Minimum rounds needed: k ≥ %d\n", min_rounds);
    printf("   - Recommended k for practice: %d (with safety margin)\n", min_rounds + 10);
}

//==============================================================================
// MAIN FUNCTION
//==============================================================================

int main(void) {
    // Initialize random number generator
    gmp_randinit_mt(global_state);
    gmp_randseed_ui(global_state, time(NULL) ^ clock());

    printf("Solovay–Strassen Primality Test - Comprehensive Analysis\n");
    printf("========================================================\n");
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

    // Generate two 256-bit primes (using SS test)
    clock_t gen_start = clock();
    generate_prime(p, PRIME_BITS, GENERATION_ROUNDS);
    generate_prime(q, PRIME_BITS, GENERATION_ROUNDS);
    clock_t gen_end = clock();

    // Compute composite n = p × q
    mpz_mul(n, p, q);

    double generation_time = ((double)(gen_end - gen_start)) / CLOCKS_PER_SEC;
    printf("\nPrime generation completed in %.2f seconds\n", generation_time);

    // Verify our generated numbers (via SS with 20 rounds)
    printf("\nVerification:\n");
    printf("  p is prime: %s\n", solovay_strassen_test(p, 20) ? "YES" : "NO");
    printf("  q is prime: %s\n", solovay_strassen_test(q, 20) ? "YES" : "NO");
    // n should be composite
    printf("  n is composite: %s\n", solovay_strassen_test(n, 20) ? "FAILED (prime)" : "CONFIRMED");

    // Run the main analysis on n (composite)
    analyze_solovay_strassen_performance(n, &stats);

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
//clang -O3 -std=c11 -flto -fomit-frame-pointer \
  -I/opt/homebrew/include -L/opt/homebrew/lib \
  -o ss_analysis solovay_strassen_analysis.c -lgmp