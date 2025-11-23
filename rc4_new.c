// rc4.c with inlining
// Fully-working, highly-optimized RC4 in C (KSA + PRGA)
// Uses only unsigned char so all arithmetic is mod 256 for free

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#ifdef __GNUC__
  #include <x86intrin.h>   // for __rdtsc() on x86_64
#else
  #include <time.h>
#endif

typedef uint8_t u8;

// RC4 context: 256-byte state array plus two 8-bit indices
typedef struct {
    u8 S[256];
    u8 i, j;
} RC4_CTX;

/**
 * rc4_init(ctx, key, keylen) KSA ROUND
 *   Perform RC4 Key-Scheduling Algorithm (KSA):
 *   - Initialize S to [0,1,2...255]
 *   - Mix in the key with wraparound via 8-bit adds (no modulus ops)
 */
void rc4_init(RC4_CTX *ctx, const u8 *key, size_t keylen) {
    for (int idx = 0; idx < 256; idx++) {
        ctx->S[idx] = (u8)idx;
    }
    u8 j = 0;
    for (int idx = 0; idx < 256; idx++) {
        j = (u8)(j + ctx->S[idx] + key[idx % keylen]);
        // swap S[idx] <- S[j]
        u8 tmp      = ctx->S[idx];
        ctx->S[idx] = ctx->S[j];
        ctx->S[j]   = tmp;
    }
    ctx->i = 0;
    ctx->j = 0;
}

/**
 * rc4_byte(ctx)
 *   Perform one PRGA step:
 *     i <- i + 1                // 1 ADD
 *     j <- j + S[i]             // 1 LOAD, 1 ADD
 *     swap S[i] and S[j]        // 2 LOAD, 2 STORE
 *     K <- S[tmp + S[j]]        // 1 ADD, 1 LOAD
 *   Keystream uses XOR:          // 1 LOAD, 1 XOR, 1 STORE (in rc4_crypt)
 *   Total per byte: 2 ADD, 3 LOAD, 2 STORE, 1 XOR
 *
 * Instruction breakdown (x86_64 estimate):
 *   inc    %al                  ; 1 cycle
 *   mov    al, S[i]             ; 2 cycles
 *   add    bl, al               ; 1 cycle
 *   mov    bl, S[j]             ; 2 cycles
 *   xchg   S[i], S[j]           ; 1 cycle
 *   lea    rdx, [tmp + S[j]]    ; 1 cycle
 *   mov    al, [S+rdx]          ; 2 cycles
 * ~10–12 cycles core + memory latency ≈18–22 cycles/byte
 * //try instruction parallelism
 * Marked always_inline to ensure compiler inlines it.
 */
static inline __attribute__((always_inline)) u8 rc4_byte(RC4_CTX *ctx) {
    u8 i = (u8)(ctx->i + 1);
    ctx->i = i;

    u8 j = (u8)(ctx->j + ctx->S[i]);
    ctx->j = j;

    u8 tmp    = ctx->S[i];
    ctx->S[i] = ctx->S[j];
    ctx->S[j] = tmp;
    //use xor based swapping to avoid extra temp variable
    //tmp = ctx->S[j];
    //ctx->S[j] = ctx->S[i];
    //ctx->S[i] = tmp;

    // K <- S[tmp + S[j]]
    // This is the keystream byte
    // tmp + S[j] is the index into S
    // S[tmp + S[j]] is the keystream byte
    // This is the core of the PRGA
    // It generates one byte of keystream per call
    // The result is used in rc4_crypt to XOR with plaintext/ciphertext
    // This is the final step of the PRGA
    // It produces the keystream byte that will be XORed with the input data
    // to produce the ciphertext or plaintext
    // This is a simple and efficient way to generate the keystream
    // It uses only 8-bit arithmetic and indexing, which is fast on modern CPUs

    return ctx->S[(u8)(tmp + ctx->S[j])];
}

/**
 * rc4_crypt(ctx, data, datalen) PRGA ROUND
 *   XOR-encrypt/decrypt data in place using rc4_byte for each byte.
 *   Measures cycles with __rdtsc() (negligible overhead).
 */
unsigned long long rc4_crypt(RC4_CTX *ctx, u8 *data, size_t datalen) {
#ifdef __GNUC__
    unsigned long long start = __rdtsc();
    for (size_t k = 0; k < datalen; k++) {
        data[k] ^= rc4_byte(ctx);
    }
    unsigned long long end = __rdtsc();
    return end - start;
#else
    clock_t start = clock();
    for (size_t k = 0; k < datalen; k++) {
        data[k] ^= rc4_byte(ctx);
    }
    clock_t end = clock();
    return (unsigned long long)(end - start);
#endif
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <key> <plaintext>\n", argv[0]);
        return 1;
    }

    const u8 *key = (const u8*)argv[1];
    size_t keylen = strlen(argv[1]);
    u8 *data      = (u8*)malloc(strlen(argv[2]) + 1);
    strcpy((char*)data, argv[2]);
    size_t datalen = strlen((char*)data);

    RC4_CTX ctx;
    rc4_init(&ctx, key, keylen);

    unsigned long long cycles = rc4_crypt(&ctx, data, datalen);

    printf("Ciphertext: ");
    for (size_t i = 0; i < datalen; i++) {
        printf("%02X", data[i]);
    }
    printf("\n");

    printf("PRGA cycles: %llu (≈ %.2f cycles/byte)\n",
           cycles, (double)cycles / (double)datalen);

    free(data);
    return 0;

    //compile with this command: $gcc -O3 -march=native -std=c11 -o rc4_new rc4_new.c /* That command invokes GCC to compile rc4_new.c with optimization level 3 targeting your native CPU and C11 standard and produces an executable named rc4_new. */
    //usage: $./rc4_new <key> <message>
}