// rsa_gmp.c
// RSA key‐generation and string encrypt/decrypt using GMP for a ~1024‐bit modulus.

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define PRIME_BITS   1024      // each prime ~512 bits → N ~1024 bits
#define MR_ROUNDS     25      // Miller–Rabin rounds for primality testing
#define MAX_MSG_LEN  120      // max bytes of input string (must fit in N)

int main(void) {
    // --- 1) Key generation ---
    gmp_randstate_t state;
    mpz_t p, q, n, phi, e, d, g, tmp;
    mpz_inits(p, q, n, phi, e, d, g, tmp, NULL);

    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));

    // generate prime p
    do {
        mpz_urandomb(p, state, PRIME_BITS);
        mpz_setbit(p, PRIME_BITS-1);
        mpz_setbit(p, 0);
    } while (!mpz_probab_prime_p(p, MR_ROUNDS));

    // generate prime q ≠ p
    do {
        mpz_urandomb(q, state, PRIME_BITS);
        mpz_setbit(q, PRIME_BITS-1);
        mpz_setbit(q, 0);
    } while (!mpz_probab_prime_p(q, MR_ROUNDS) || mpz_cmp(p, q) == 0);

    // n = p * q
    mpz_mul(n, p, q);
    printf("Generated primes:\n");
    printf("p = %s\n", mpz_get_str(NULL, 16, p));
    printf("q = %s\n\n", mpz_get_str(NULL, 16, q));

    // φ(n) = (p-1)*(q-1)
    mpz_sub_ui(tmp, p, 1);
    mpz_sub_ui(phi, q, 1);
    mpz_mul(phi, phi, tmp);

    // e = 65537
    mpz_set_ui(e, 65537);

    // ensure gcd(e, φ(n)) = 1
    mpz_gcd(g, e, phi);
    if (mpz_cmp_ui(g, 1) != 0) {
        fprintf(stderr, "ERROR: e not coprime with φ(n)\n");
        return 1;
    }

    // d = e⁻¹ mod φ(n)
    if (!mpz_invert(d, e, phi)) {
        fprintf(stderr, "ERROR: modular inverse failed\n");
        return 1;
    }

    // print key parameters via mpz_get_str
    char *s_n = mpz_get_str(NULL, 16, n);
    char *s_e = mpz_get_str(NULL, 16, e);
    char *s_d = mpz_get_str(NULL, 16, d);
    printf("n = %s\n\n", s_n);
    printf("e = %s\n\n", s_e);
    printf("d = %s\n\n", s_d);
    free(s_n); free(s_e); free(s_d);

    // --- 2) Read plaintext string ---
    char msg[MAX_MSG_LEN+1];
    printf("Enter message (max %d chars):\n", MAX_MSG_LEN);
    if (!fgets(msg, sizeof(msg), stdin)) return 1;
    size_t msg_len = strnlen(msg, sizeof(msg));
    if (msg_len > 0 && msg[msg_len-1] == '\n') msg[--msg_len] = '\0';

    // import text into integer m
    mpz_t m;
    mpz_init(m);
    mpz_import(m, msg_len, 1, 1, 0, 0, msg);
    if (mpz_cmp(m, n) >= 0) {
        fprintf(stderr, "ERROR: message too large for modulus\n");
        return 1;
    }

    // --- 3) Encrypt: c = m^e mod n ---
    mpz_t c;
    mpz_init(c);
    mpz_powm(c, m, e, n);
    char *s_c = mpz_get_str(NULL, 16, c);
    printf("Ciphertext (hex): %s\n\n", s_c);
    free(s_c);

    // --- 4) Decrypt: m2 = c^d mod n ---
    mpz_t m2;
    mpz_init(m2);
    mpz_powm(m2, c, d, n);

    // export integer back to bytes
    size_t out_len;
    unsigned char *out = calloc(1, msg_len + 1);
    mpz_export(out, &out_len, 1, 1, 0, 0, m2);

    // print decrypted text
    printf("Decrypted message:\n%.*s\n", (int)out_len, out);

    // cleanup
    free(out);
    mpz_clears(p, q, n, phi, e, d, g, tmp, m, c, m2, NULL);
    gmp_randclear(state);
    return 0;
}