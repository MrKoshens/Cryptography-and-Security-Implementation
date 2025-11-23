#include <stdint.h>
#include <stddef.h>
#include <string.h>

// rotate left 32-bit
static inline uint32_t ROTL32(uint32_t x, int r) {
    return (x << r) | (x >> (32 - r));
}

// quarter round as defined in RFC 8439
//chacha20 quarter round macro
#define QUARTERROUND(a, b, c, d) \
    do {                         \
        a += b;                  \
        d ^= a;                  \
        d = ROTL32(d, 16);       \
        c += d;                  \
        b ^= c;                  \
        b = ROTL32(b, 12);       \
        a += b;                  \
        d ^= a;                  \
        d = ROTL32(d, 8);        \
        c += d;                  \
        b ^= c;                  \
        b = ROTL32(b, 7);        \
    } while (0)

// load 32-bit little-endian
static inline uint32_t load32_le(const uint8_t *p) {
    return (uint32_t)p[0]
         | ((uint32_t)p[1] << 8)
         | ((uint32_t)p[2] << 16)
         | ((uint32_t)p[3] << 24);
}

// store 32-bit little-endian
static inline void store32_le(uint8_t *p, uint32_t v) {
    p[0] = (uint8_t)(v);
    p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
}

// ChaCha20 block function: input is 16 words (constants, key, counter, nonce)
// output is 64 bytes (16 words little endian)
static void chacha20_block(uint8_t out[64], const uint32_t input[16]) {
    // copy into locals so compiler can keep in registers
    uint32_t x0  = input[0];
    uint32_t x1  = input[1];
    uint32_t x2  = input[2];
    uint32_t x3  = input[3];
    uint32_t x4  = input[4];
    uint32_t x5  = input[5];
    uint32_t x6  = input[6];
    uint32_t x7  = input[7];
    uint32_t x8  = input[8];
    uint32_t x9  = input[9];
    uint32_t x10 = input[10];
    uint32_t x11 = input[11];
    uint32_t x12 = input[12];
    uint32_t x13 = input[13];
    uint32_t x14 = input[14];
    uint32_t x15 = input[15];

    // 20 rounds = 10 iterations of double rounds
    for (int i = 0; i < 10; i++) {
        // column rounds
        QUARTERROUND(x0, x4, x8, x12);
        QUARTERROUND(x1, x5, x9, x13);
        QUARTERROUND(x2, x6, x10, x14);
        QUARTERROUND(x3, x7, x11, x15);
        // diagonal rounds
        QUARTERROUND(x0, x5, x10, x15);
        QUARTERROUND(x1, x6, x11, x12);
        QUARTERROUND(x2, x7, x8, x13);
        QUARTERROUND(x3, x4, x9, x14);
    }

    // add original state and serialize to output in little endian
    store32_le(out +  0, x0  + input[0]);
    store32_le(out +  4, x1  + input[1]);
    store32_le(out +  8, x2  + input[2]);
    store32_le(out + 12, x3  + input[3]);
    store32_le(out + 16, x4  + input[4]);
    store32_le(out + 20, x5  + input[5]);
    store32_le(out + 24, x6  + input[6]);
    store32_le(out + 28, x7  + input[7]);
    store32_le(out + 32, x8  + input[8]);
    store32_le(out + 36, x9  + input[9]);
    store32_le(out + 40, x10 + input[10]);
    store32_le(out + 44, x11 + input[11]);
    store32_le(out + 48, x12 + input[12]);
    store32_le(out + 52, x13 + input[13]);
    store32_le(out + 56, x14 + input[14]);
    store32_le(out + 60, x15 + input[15]);
}

// High-level ChaCha20 XOR: encrypt/decrypt `len` bytes of `in` into `out` using key, counter, nonce
// key: 32 bytes, nonce: 12 bytes, counter: 32-bit
void chacha20_xor(uint8_t *out, const uint8_t *in, size_t len,
                  const uint8_t key[32], const uint8_t nonce[12], uint32_t counter) {
    // constants "expand 32-byte k"
    const uint8_t *constants = (const uint8_t *)"expa"
                               "nd 3"
                               "2-by"
                               "te k"; // little endian pieces

    uint32_t state[16];

    // load constant words (little endian)
    state[0] = load32_le((const uint8_t *)"expa");
    state[1] = load32_le((const uint8_t *)"nd 3");
    state[2] = load32_le((const uint8_t *)"2-by");
    state[3] = load32_le((const uint8_t *)"te k");

    // load key
    for (int i = 0; i < 8; i++) {
        state[4 + i] = load32_le(key + 4 * i);
    }

    // counter and nonce
    state[12] = counter;
    state[13] = load32_le(nonce + 0);
    state[14] = load32_le(nonce + 4);
    state[15] = load32_le(nonce + 8);

    uint8_t block[64];
    size_t offset = 0;

    while (len > 0) {
        chacha20_block(block, state);

        size_t chunk = (len < 64) ? len : 64;
        for (size_t i = 0; i < chunk; i++) {
            out[offset + i] = in[offset + i] ^ block[i];
        }

        len -= chunk;
        offset += chunk;

        // increment counter (wraparound allowed)
        state[12]++;
    }
}