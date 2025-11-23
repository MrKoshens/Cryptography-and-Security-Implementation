// chacha20_simd.c
// A complete ChaCha20 implementation with scalar, interleaved-scalar,
// and NEON 4-way SIMD paths, with compilation fixed for Apple M1.

//chacha20_simd.c – comprehensive ChaCha20 stream‐cipher implementation:
//   chacha20_block: single 64-byte block core (20 rounds of mixing)
//   chacha20_xor_interleaved4: 4-block interleaved scalar path to hide instruction and memory latency
//   chacha20_xor_neon4: 4-way SIMD via ARM NEON for maximum throughput on M1
//   chacha20_xor_best: dispatch helper choosing the fastest available path

#include <stdint.h>
#include <stddef.h>
#include <string.h>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
  #include <arm_neon.h>
  #define HAVE_NEON 1
#else
  #define HAVE_NEON 0
#endif

//ROTL32(x, r):
//  Perform a 32-bit circular left rotation of x by r bits.
//  This is the basic bit-mixing primitive in ChaCha, allowing low-cost nonlinear diffusion
//  by repositioning bits within a word before XOR/add operations
// 32-bit left rotate
static inline uint32_t ROTL32(uint32_t x, int r) {
    return (x << r) | (x >> (32 - r));
}

// little-endian load/store
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

//QUARTERROUND(a, b, c, d):
  //The fundamental mixing step from RFC 8439. Takes four 32-bit words (a,b,c,d) and applies:
  //  1. a += b; d ^= a; d <<< 16
  // 2. c += d; b ^= c; b <<< 12
  // 3. a += b; d ^= a; d <<<  8
  //  4. c += d; b ^= c; b <<<  7
  //This sequence provides both inter-word mixing (adds/XORs) and intra-word diffusion (rotates).
//
// This is the core of ChaCha20's mixing function, applied in rounds to achieve diffusion.
// It operates on four 32-bit words, performing a series of additions, XORs, and rotations.
// The rotations are fixed to specific bit counts (16, 12, 8, 7) as per the RFC.
// The order of operations ensures that each word influences the others, creating a complex
// nonlinear transformation that is difficult to reverse without the key.

// Quarter round as per RFC 8439
// This macro implements the quarter-round operation on four 32-bit words
// using C macros. It takes four 32x4 vectors (a, b, c, d) and performs
// the same operations as the scalar QUARTERROUND but in parallel across all
// four lanes of each vector. This allows for high throughput on ARM64 CPUs
// with NEON support by processing four words simultaneously.


#define QUARTERROUND(a, b, c, d)       \
    do {                               \
        a += b;                        \
        d ^= a;                        \
        d = ROTL32(d, 16);             \
        c += d;                        \
        b ^= c;                        \
        b = ROTL32(b, 12);             \
        a += b;                        \
        d ^= a;                        \
        d = ROTL32(d, 8);              \
        c += d;                        \
        b ^= c;                        \
        b = ROTL32(b, 7);              \
    } while (0)

// ChaCha20 block function: input is 16 words (constants, key, counter, nonce)
// output is 64 bytes (16 words little endian)
// This function performs 20 rounds of mixing on the input state, which consists of
//chacha20_block(out, input):
//  - Copy the 16-word input into locals x0…x15 so the compiler keeps state in registers.
//  - Perform 10 double-rounds (20 total):
//      • Column rounds on (0,4,8,12)…(3,7,11,15)
//      • Diagonal rounds on (0,5,10,15)…(3,4,9,14)
//  - Feed-forward: add original input words back into x0…x15.
// - Serialize x0…x15 into the 64-byte output buffer in little-endian order.
//  This gives you one 64-byte keystream block
//  that can be XORed with plaintext/ciphertext to encrypt/decrypt.

// Single-block ChaCha20 core: outputs 64-byte keystream
static void chacha20_block(uint8_t out[64], const uint32_t input[16]) {
    uint32_t x0 = input[0],  x1 = input[1],  x2 = input[2],  x3 = input[3];
    uint32_t x4 = input[4],  x5 = input[5],  x6 = input[6],  x7 = input[7];
    uint32_t x8 = input[8],  x9 = input[9],  x10 = input[10], x11 = input[11];
    uint32_t x12 = input[12], x13 = input[13], x14 = input[14], x15 = input[15];

    for (int i = 0; i < 10; i++) {
        // Column rounds
        QUARTERROUND(x0, x4, x8,  x12);
        QUARTERROUND(x1, x5, x9,  x13);
        QUARTERROUND(x2, x6, x10, x14);
        QUARTERROUND(x3, x7, x11, x15);
        // Diagonal rounds
        QUARTERROUND(x0, x5, x10, x15);
        QUARTERROUND(x1, x6, x11, x12);
        QUARTERROUND(x2, x7, x8,  x13);
        QUARTERROUND(x3, x4, x9,  x14);
    }

    // Feed-forward and serialize to out in little-endian
    // This adds the original input state back to the mixed state
    // and stores the result in little-endian order
    // This is the final step that combines the original input with the mixed state
    // This ensures that the output is a function of both the input state and the mixing process
    store32_le(out +  0, x0 + input[0]);
    store32_le(out +  4, x1 + input[1]);
    store32_le(out +  8, x2 + input[2]);
    store32_le(out + 12, x3 + input[3]);
    store32_le(out + 16, x4 + input[4]);
    store32_le(out + 20, x5 + input[5]);
    store32_le(out + 24, x6 + input[6]);
    store32_le(out + 28, x7 + input[7]);
    store32_le(out + 32, x8 + input[8]);
    store32_le(out + 36, x9 + input[9]);
    store32_le(out + 40, x10 + input[10]);
    store32_le(out + 44, x11 + input[11]);
    store32_le(out + 48, x12 + input[12]);
    store32_le(out + 52, x13 + input[13]);
    store32_le(out + 56, x14 + input[14]);
    store32_le(out + 60, x15 + input[15]);
}

//chacha20_init_state(state, key, nonce, counter):
//  Initialize the 4×4 ChaCha matrix:
//    state[0..3]  = constant bytes of “expand 32-byte k”
//    state[4..11] = eight 32-bit words of the 256-bit key (loaded little-endian)
//    state[12]    = 32-bit block counter (increments each block)
//    state[13..15] = three 32-bit words of the 96-bit nonce
//  This layout matches RFC 8439 exactly.

// Initialize state: constants, key, counter, nonce
static void chacha20_init_state(uint32_t state[16],
                                const uint8_t key[32],
                                const uint8_t nonce[12],
                                uint32_t counter) {
    // Constants "expand 32-byte k"
    // These are the first four 32-bit words of the ChaCha20 state
    // They are fixed and defined in RFC 8439
    // The constants are little-endian encoded as:
    state[0]  = 0x61707865;
    state[1]  = 0x3320646e;
    state[2]  = 0x79622d32;
    state[3]  = 0x6b206574;
    // Key
    for (int i = 0; i < 8; i++) {
        state[4 + i] = load32_le(key + 4 * i);
    }
    // Counter
    state[12] = counter;
    // Nonce
    state[13] = load32_le(nonce + 0);
    state[14] = load32_le(nonce + 4);
    state[15] = load32_le(nonce + 8);
}

//chacha20_xor_interleaved4(out, in, len, …):
//  Process data in 256-byte chunks (4×64):
//    1. For k = 0…3: init state with counter+k, call chacha20_block -> block[k]
//    2. XOR each 64-byte block[k] with its corresponding input segment
//  This dead-simple interleaving hides the ~20-round latency by keeping four independent
//  ChaCha blocks in flight, improving throughput on out-of-order CPUs without SIMD.
// // This is a scalar path that interleaves four blocks to hide latency.
// 4-way interleaved scalar path
void chacha20_xor_interleaved4(uint8_t *out, const uint8_t *in, size_t len,
                               const uint8_t key[32], const uint8_t nonce[12],
                               uint32_t counter) {
    uint8_t block[4][64];
    size_t off = 0;
    while (len >= 256) {
        // Generate 4 blocks
        for (int k = 0; k < 4; k++) {
            uint32_t st[16];
            chacha20_init_state(st, key, nonce, counter + k);
            chacha20_block(block[k], st);
        }
        // XOR
        for (int k = 0; k < 4; k++) {
            for (int i = 0; i < 64; i++) {
                out[off + k * 64 + i] = in[off + k * 64 + i] ^ block[k][i];
            }
        }
        off += 256;
        len -= 256;
        counter += 4;
    }
    // Tail fallback
    // Tail fallback: if remaining data < 256 bytes, fall back to single-block scalar mode.
    // This loop handles the last 1–255 bytes correctly without overrun.
    while (len > 0) {
        uint32_t st[16];
        uint8_t buf[64];
        size_t c = (len < 64 ? len : 64);
        chacha20_init_state(st, key, nonce, counter);
        chacha20_block(buf, st);
        for (size_t i = 0; i < c; i++) {
            out[off + i] = in[off + i] ^ buf[i];
        }
        off += c;
        len -= c;
        counter++;
    }
}

#if HAVE_NEON
// NEON fixed rotates
// NEON fixed-rotate macros (ROTL16/12/8/7):
//   ARM NEON intrinsics require constant shift amounts. These macros implement
//   the four rotate constants used by the quarter‐round, all in vector form.
#define ROTL16(v) vorrq_u32(vshlq_n_u32(v, 16), vshrq_n_u32(v, 16))
#define ROTL12(v) vorrq_u32(vshlq_n_u32(v, 12), vshrq_n_u32(v, 20))
#define ROTL8(v)  vorrq_u32(vshlq_n_u32(v, 8),  vshrq_n_u32(v, 24))
#define ROTL7(v)  vorrq_u32(vshlq_n_u32(v, 7),  vshrq_n_u32(v, 25))

// NEON quarter round on 4 lanes
// NEON_QR(a, b, c, d):
//   This macro implements the quarter-round operation on four 32-bit words
//   using NEON intrinsics. It takes four 32x4 vectors (a, b, c, d) and performs
//   the same operations as the scalar QUARTERROUND but in parallel across all
//   four lanes of each vector. This allows for high throughput on ARM64 CPUs
//   with NEON support by processing four words simultaneously.
#define NEON_QR(a, b, c, d)            \
    do {                               \
        a = vaddq_u32(a, b);           \
        d = veorq_u32(d, a);           \
        d = ROTL16(d);                 \
        c = vaddq_u32(c, d);           \
        b = veorq_u32(b, c);           \
        b = ROTL12(b);                 \
        a = vaddq_u32(a, b);           \
        d = veorq_u32(d, a);           \
        d = ROTL8(d);                  \
        c = vaddq_u32(c, d);           \
        b = veorq_u32(b, c);           \
        b = ROTL7(b);                  \
    } while (0)


//chacha20_xor_neon4(out, in, len, …):
//  Vectorized 4-way ChaCha20 using NEON registers:
//    - Broadcast constants, key words, counter vector, and nonce into 16 uint32x4 vectors
//    - Perform 20 rounds by invoking NEON_QR on columns & diagonals
//    - Feed-forward: add the original input vectors back
//    - Gather lanes 0–3 from each vector to form four 64-byte keystream blocks
//    - XOR them with four input segments in lockstep
//  This achieves ~4× the per-byte throughput of scalar code on ARM64.    
// NEON 4-way SIMD path
void chacha20_xor_neon4(uint8_t *out, const uint8_t *in, size_t len,
                       const uint8_t key[32], const uint8_t nonce[12],
                       uint32_t counter) {
    size_t off = 0;
    uint8_t ks[4][64];
    while (len >= 256) {
        // Load input state into vectors
        // Broadcast constants, key, counter, nonce
        // 16 vectors: 4x4 matrix
        // Each vector holds 4 words (16 bytes)
        uint32x4_t in0  = vdupq_n_u32(0x61707865);
        uint32x4_t in1  = vdupq_n_u32(0x3320646e);
        uint32x4_t in2  = vdupq_n_u32(0x79622d32);
        uint32x4_t in3  = vdupq_n_u32(0x6b206574);
        uint32x4_t in4  = vdupq_n_u32(load32_le(key + 0));
        uint32x4_t in5  = vdupq_n_u32(load32_le(key + 4));
        uint32x4_t in6  = vdupq_n_u32(load32_le(key + 8));
        uint32x4_t in7  = vdupq_n_u32(load32_le(key + 12));
        uint32x4_t in8  = vdupq_n_u32(load32_le(key + 16));
        uint32x4_t in9  = vdupq_n_u32(load32_le(key + 20));
        uint32x4_t in10 = vdupq_n_u32(load32_le(key + 24));
        uint32x4_t in11 = vdupq_n_u32(load32_le(key + 28));
        uint32x4_t in12 = { counter, counter+1, counter+2, counter+3 };
        uint32x4_t in13 = vdupq_n_u32(load32_le(nonce + 0));
        uint32x4_t in14 = vdupq_n_u32(load32_le(nonce + 4));
        uint32x4_t in15 = vdupq_n_u32(load32_le(nonce + 8));
        // Working copy
        uint32x4_t x0=in0,  x1=in1,  x2=in2,  x3=in3;
        uint32x4_t x4=in4,  x5=in5,  x6=in6,  x7=in7;
        uint32x4_t x8=in8,  x9=in9,  x10=in10,x11=in11;
        uint32x4_t x12=in12,x13=in13,x14=in14,x15=in15;
        // 20 rounds
        for (int i = 0; i < 10; i++) {
            NEON_QR(x0, x4, x8,  x12);
            NEON_QR(x1, x5, x9,  x13);
            NEON_QR(x2, x6, x10, x14);
            NEON_QR(x3, x7, x11, x15);
            NEON_QR(x0, x5, x10, x15);
            NEON_QR(x1, x6, x11, x12);
            NEON_QR(x2, x7, x8,  x13);
            NEON_QR(x3, x4, x9,  x14);
        }
        // Feed-forward
        x0  = vaddq_u32(x0, in0);
        x1  = vaddq_u32(x1, in1);
        x2  = vaddq_u32(x2, in2);
        x3  = vaddq_u32(x3, in3);
        x4  = vaddq_u32(x4, in4);
        x5  = vaddq_u32(x5, in5);
        x6  = vaddq_u32(x6, in6);
        x7  = vaddq_u32(x7, in7);
        x8  = vaddq_u32(x8, in8);
        x9  = vaddq_u32(x9, in9);
        x10 = vaddq_u32(x10, in10);
        x11 = vaddq_u32(x11, in11);
        x12 = vaddq_u32(x12, in12);
        x13 = vaddq_u32(x13, in13);
        x14 = vaddq_u32(x14, in14);
        x15 = vaddq_u32(x15, in15);
        // Gather into array of 16 vectors
        uint32x4_t xs[16] = {x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15};
        // Extract lanes for each of 4 streams
        for (int w = 0; w < 16; w++) {
            store32_le(ks[0] + 4*w, vgetq_lane_u32(xs[w], 0));
            store32_le(ks[1] + 4*w, vgetq_lane_u32(xs[w], 1));
            store32_le(ks[2] + 4*w, vgetq_lane_u32(xs[w], 2));
            store32_le(ks[3] + 4*w, vgetq_lane_u32(xs[w], 3));
        }
        // XOR 4 blocks
        for (int k = 0; k < 4; k++) {
            for (int i = 0; i < 64; i++) {
                out[off + k*64 + i] = in[off + k*64 + i] ^ ks[k][i];
            }
        }
        off += 256;
        len -= 256;
        counter += 4;
    }
    // Tail fallback
    if (len > 0) {
        chacha20_xor_interleaved4(out + off, in + off, len, key, nonce, counter);
    }
}
#endif

//chacha20_xor_best(...):
//  Runtime dispatcher that selects the fastest path:
//    • NEON 4-way SIMD if available
//    • otherwise 4-way interleaved scalar
//  Provides a single API for encryption/decryption without burdening callers
//  with hardware-specific details.
// This function is the entry point for users, automatically choosing the best implementation.
// It checks for NEON support and calls the appropriate function based on availability.
// It abstracts away the complexity of choosing between SIMD and scalar paths,
// allowing users to simply call this function without worrying about the underlying details.
// Dispatcher: choose NEON path or scalar interleaved
void chacha20_xor_best(uint8_t *out, const uint8_t *in, size_t len,
                       const uint8_t key[32], const uint8_t nonce[12],
                       uint32_t counter) {
#if HAVE_NEON
    chacha20_xor_neon4(out, in, len, key, nonce, counter);
#else
    chacha20_xor_interleaved4(out, in, len, key, nonce, counter);
#endif
}


// compile with this on M1: $clang -O3 -mcpu=apple-m1 -std=c11 -o chacha20_simd chacha20_simd.c