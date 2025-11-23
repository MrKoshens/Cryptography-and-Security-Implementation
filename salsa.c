// salsa20.c
// Simple Salsa20/20 implementation: core + single‐stream XOR

#include <stdint.h>
#include <stddef.h>
#include <string.h>

// 32‐bit left rotate
static inline uint32_t ROTL32(uint32_t x, int r) {
    return (x << r) | (x >> (32 - r));
}

// load/store little‐endian 32‐bit words
static inline uint32_t load32_le(const uint8_t *p) {
    return (uint32_t)p[0]
         | ((uint32_t)p[1] << 8)
         | ((uint32_t)p[2] << 16)
         | ((uint32_t)p[3] << 24);
}
static inline void store32_le(uint8_t *p, uint32_t v) {
    p[0] = (uint8_t)v;
    p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
}

/*
 * salsa20_core(out, in)
 *   input: 16 × 32-bit words  = constant + key + counter + nonce
 *   output: 64 bytes = 16 words little‐endian
 *   Performs 20 rounds = 10 double‐rounds of the Salsa quarter‐round.
 */
static void salsa20_core(uint8_t out[64], const uint32_t in[16]) {
    // copy input to working regs
    uint32_t x0 = in[0],  x1 = in[1],  x2 = in[2],  x3 = in[3];
    uint32_t x4 = in[4],  x5 = in[5],  x6 = in[6],  x7 = in[7];
    uint32_t x8 = in[8],  x9 = in[9],  x10 = in[10], x11 = in[11];
    uint32_t x12 = in[12], x13 = in[13], x14 = in[14], x15 = in[15];

    for (int i = 0; i < 10; i++) {
        // odd round: columns
        x4 ^= ROTL32(x0 + x12, 7);
        x8 ^= ROTL32(x4 + x0, 9);
        x12 ^= ROTL32(x8 + x4, 13);
        x0 ^= ROTL32(x12 + x8, 18);

        x9 ^= ROTL32(x5 + x1, 7);
        x13 ^= ROTL32(x9 + x5, 9);
        x1 ^= ROTL32(x13 + x9, 13);
        x5 ^= ROTL32(x1 + x13, 18);

        x14 ^= ROTL32(x10 + x6, 7);
        x2 ^= ROTL32(x14 + x10, 9);
        x6 ^= ROTL32(x2 + x14, 13);
        x10 ^= ROTL32(x6 + x2, 18);

        x3 ^= ROTL32(x15 + x11, 7);
        x7 ^= ROTL32(x3 + x15, 9);
        x11 ^= ROTL32(x7 + x3, 13);
        x15 ^= ROTL32(x11 + x7, 18);

        // even round: diagonals
        x1 ^= ROTL32(x0 + x3, 7);
        x2 ^= ROTL32(x1 + x0, 9);
        x3 ^= ROTL32(x2 + x1, 13);
        x0 ^= ROTL32(x3 + x2, 18);

        x6 ^= ROTL32(x5 + x4, 7);
        x7 ^= ROTL32(x6 + x5, 9);
        x4 ^= ROTL32(x7 + x6, 13);
        x5 ^= ROTL32(x4 + x7, 18);

        x11 ^= ROTL32(x10 + x9, 7);
        x8  ^= ROTL32(x11 + x10, 9);
        x9  ^= ROTL32(x8  + x11, 13);
        x10 ^= ROTL32(x9  + x8, 18);

        x12 ^= ROTL32(x15 + x14, 7);
        x13 ^= ROTL32(x12 + x15, 9);
        x14 ^= ROTL32(x13 + x12, 13);
        x15 ^= ROTL32(x14 + x13, 18);
    }

    // add original input (feed‐forward) and serialize
    store32_le(out +  0, x0  + in[0]);
    store32_le(out +  4, x1  + in[1]);
    store32_le(out +  8, x2  + in[2]);
    store32_le(out + 12, x3  + in[3]);
    store32_le(out + 16, x4  + in[4]);
    store32_le(out + 20, x5  + in[5]);
    store32_le(out + 24, x6  + in[6]);
    store32_le(out + 28, x7  + in[7]);
    store32_le(out + 32, x8  + in[8]);
    store32_le(out + 36, x9  + in[9]);
    store32_le(out + 40, x10 + in[10]);
    store32_le(out + 44, x11 + in[11]);
    store32_le(out + 48, x12 + in[12]);
    store32_le(out + 52, x13 + in[13]);
    store32_le(out + 56, x14 + in[14]);
    store32_le(out + 60, x15 + in[15]);
}

/*
 * salsa20_xor(out, in, len, key, nonce, counter)
 *   High‐level API: XOR‐encrypt/decrypt `len` bytes of `in` into `out`
 *   using a 256‐bit key, 96‐bit nonce, and 32‐bit counter.
 */
void salsa20_xor(uint8_t *out,
                 const uint8_t *in,
                 size_t len,
                 const uint8_t key[32],
                 const uint8_t nonce[12],
                 uint32_t counter)
{
    uint32_t state[16];
    uint8_t block[64];
    size_t off = 0;

    // constants "expand 32-byte k"
    const uint8_t *cstr = (const uint8_t *)"expand 32-byte k";

    // build initial state template
    state[0] = load32_le(cstr + 0);
    state[1] = load32_le(cstr + 4);
    state[2] = load32_le(cstr + 8);
    state[3] = load32_le(cstr + 12);
    // load key
    for (int i = 0; i < 8; i++) {
        state[4 + i] = load32_le(key + 4 * i);
    }
    // set nonce words
    state[13] = load32_le(nonce + 0);
    state[14] = load32_le(nonce + 4);
    state[15] = load32_le(nonce + 8);

    while (len > 0) {
        // set counter word
        state[12] = counter++;
        // generate keystream block
        salsa20_core(block, state);
        // XOR up to 64 bytes
        size_t chunk = (len < 64 ? len : 64);
        for (size_t i = 0; i < chunk; i++) {
            out[off + i] = in[off + i] ^ block[i];
        }
        off += chunk;
        len -= chunk;
    }
}