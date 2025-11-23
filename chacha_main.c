#include <stdio.h>
#include <stdint.h>
#include <string.h>

// forward declaration
void chacha20_xor_best(uint8_t *out, const uint8_t *in, size_t len,
                       const uint8_t key[32], const uint8_t nonce[12], uint32_t counter);

int main(void) {
    // 256-bit key (here all zero for demo) and 96-bit nonce
    uint8_t key[32]   = {0};
    uint8_t nonce[12] = {0};

    // plaintext
    const char *pt = "Cryptography and Security Implementation!";
    size_t len = strlen(pt);
    uint8_t out[len];

    // encrypt with counter=1
    chacha20_xor_best(out, (const uint8_t*)pt, len, key, nonce, 1);

    // print ciphertext as hex
    for (size_t i = 0; i < len; i++) {
        printf("%02x", out[i]);
    }
    printf("\n");
    return 0;
}

//clang -O3 -mcpu=apple-m1 -std=c11 -o test_chacha chacha_main.c chacha20_simd.c