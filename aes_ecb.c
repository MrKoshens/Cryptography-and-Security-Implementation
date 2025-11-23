/* aes_ecb_interactive.c
 *
 * AES-128 encryption and decryption primitives with an interactive CLI.
 * Supports:
 *   - Key expansion
 *   - SubBytes, ShiftRows, MixColumns, AddRoundKey
 *   - Single-block AES-128 encrypt/decrypt
 *   - Multi-block ECB encrypt with PKCS#7 padding
 *
 * What the program does
 *   1) Asks whether you want a single block test or multi-block ECB
 *   2) For single block
 *        - Asks plaintext hex and key hex
 *        - Encrypts one 16-byte block
 *        - If inputs match the FIPS-197 test vector, verifies the expected ciphertext
 *   3) For multi-block
 *        - Asks plaintext hex and key hex
 *        - Encrypts the full message in ECB with PKCS#7 padding
 *        - Prints ciphertext hex
 *
 * Build:  gcc -O3 -std=c11 -Wall -Wextra -o aes_ecb aes_ecb_interactive.c
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <ctype.h>

#define AES_BLOCK_BYTES 16
#define AES128_ROUNDS   10
#define AES128_KEY_BYTES 16
#define AES128_EXPANDED_KEY_BYTES (AES_BLOCK_BYTES * (AES128_ROUNDS + 1)) /* 176 */

static uint8_t SBOX[256];
static uint8_t INV_SBOX[256];

static inline uint8_t rotl8(uint8_t x, unsigned r) {
    return (uint8_t)((x << r) | (x >> (8 - r)));
}

/* GF(2^8) multiply with AES modulus x^8 + x^4 + x^3 + x + 1 (0x11B) */
static inline uint8_t gf_mul(uint8_t a, uint8_t b) {
    uint8_t res = 0;
    while (b) {
        if (b & 1) res ^= a;
        uint8_t hi = a & 0x80;
        a <<= 1;
        if (hi) a ^= 0x1B;
        b >>= 1;
    }
    return res;
}

static uint8_t gf_pow(uint8_t a, uint16_t e) {
    uint8_t r = 1;
    while (e) {
        if (e & 1) r = gf_mul(r, a);
        a = gf_mul(a, a);
        e >>= 1;
    }
    return r;
}

static uint8_t gf_inv(uint8_t a) {
    if (a == 0) return 0;
    /* In GF(2^8), a^(2^8 - 2) = a^254 is the multiplicative inverse */
    return gf_pow(a, 254);
}

/* Build SBOX and INV_SBOX using the canonical definition */
static void aes_init_sboxes(void) {
    for (int x = 0; x < 256; ++x) {
        uint8_t inv = gf_inv((uint8_t)x);
        uint8_t y = inv
                  ^ rotl8(inv, 1)
                  ^ rotl8(inv, 2)
                  ^ rotl8(inv, 3)
                  ^ rotl8(inv, 4)
                  ^ 0x63;
        SBOX[x] = y;
        INV_SBOX[y] = (uint8_t)x;
    }
}

/* Rcon for AES-128 rounds 1..10, first byte only. Others are zero. */
static const uint8_t RCON[11] = {
    0x00, /* unused */
    0x01, 0x02, 0x04, 0x08, 0x10,
    0x20, 0x40, 0x80, 0x1B, 0x36
};

/* Key expansion: input 16-byte key to 176-byte expanded round keys */
static void aes128_expand_key(const uint8_t key[16], uint8_t round_keys[176]) {
    memcpy(round_keys, key, 16);

    uint8_t temp[4];
    uint8_t *w = round_keys;

    int bytes_generated = 16;
    int rcon_idx = 1;

    while (bytes_generated < 176) {
        temp[0] = w[bytes_generated - 4];
        temp[1] = w[bytes_generated - 3];
        temp[2] = w[bytes_generated - 2];
        temp[3] = w[bytes_generated - 1];

        if (bytes_generated % 16 == 0) {
            /* RotWord */
            uint8_t t = temp[0];
            temp[0] = temp[1];
            temp[1] = temp[2];
            temp[2] = temp[3];
            temp[3] = t;
            /* SubWord */
            temp[0] = SBOX[temp[0]];
            temp[1] = SBOX[temp[1]];
            temp[2] = SBOX[temp[2]];
            temp[3] = SBOX[temp[3]];
            /* XOR Rcon */
            temp[0] ^= RCON[rcon_idx++];
        }

        for (int i = 0; i < 4; ++i) {
            w[bytes_generated] = w[bytes_generated - 16] ^ temp[i];
            ++bytes_generated;
        }
    }
}

/* State is 4x4 bytes in column-major order: state[row][col] */
typedef uint8_t state_t[4][4];

static void load_state(state_t s, const uint8_t in[16]) {
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            s[r][c] = in[4 * c + r];
}

static void store_state(uint8_t out[16], const state_t s) {
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            out[4 * c + r] = s[r][c];
}

static void add_round_key(state_t s, const uint8_t round_key[16]) {
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            s[r][c] ^= round_key[4 * c + r];
}

static void sub_bytes(state_t s) {
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            s[r][c] = SBOX[s[r][c]];
}

static void inv_sub_bytes(state_t s) {
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            s[r][c] = INV_SBOX[s[r][c]];
}

static void shift_rows(state_t s) {
    uint8_t t;

    t = s[1][0];
    s[1][0] = s[1][1];
    s[1][1] = s[1][2];
    s[1][2] = s[1][3];
    s[1][3] = t;

    t = s[2][0];
    s[2][0] = s[2][2];
    s[2][2] = t;
    t = s[2][1];
    s[2][1] = s[2][3];
    s[2][3] = t;

    t = s[3][3];
    s[3][3] = s[3][2];
    s[3][2] = s[3][1];
    s[3][1] = s[3][0];
    s[3][0] = t;
}

static void inv_shift_rows(state_t s) {
    uint8_t t;

    t = s[1][3];
    s[1][3] = s[1][2];
    s[1][2] = s[1][1];
    s[1][1] = s[1][0];
    s[1][0] = t;

    t = s[2][0];
    s[2][0] = s[2][2];
    s[2][2] = t;
    t = s[2][1];
    s[2][1] = s[2][3];
    s[2][3] = t;

    t = s[3][0];
    s[3][0] = s[3][1];
    s[3][1] = s[3][2];
    s[3][2] = s[3][3];
    s[3][3] = t;
}

static void mix_single_column(uint8_t *a0, uint8_t *a1, uint8_t *a2, uint8_t *a3) {
    uint8_t s0 = *a0, s1 = *a1, s2 = *a2, s3 = *a3;

    uint8_t r0 = gf_mul(0x02, s0) ^ gf_mul(0x03, s1) ^ s2 ^ s3;
    uint8_t r1 = s0 ^ gf_mul(0x02, s1) ^ gf_mul(0x03, s2) ^ s3;
    uint8_t r2 = s0 ^ s1 ^ gf_mul(0x02, s2) ^ gf_mul(0x03, s3);
    uint8_t r3 = gf_mul(0x03, s0) ^ s1 ^ s2 ^ gf_mul(0x02, s3);

    *a0 = r0; *a1 = r1; *a2 = r2; *a3 = r3;
}

static void mix_columns(state_t s) {
    for (int c = 0; c < 4; ++c) {
        mix_single_column(&s[0][c], &s[1][c], &s[2][c], &s[3][c]);
    }
}

static void inv_mix_single_column(uint8_t *a0, uint8_t *a1, uint8_t *a2, uint8_t *a3) {
    uint8_t s0 = *a0, s1 = *a1, s2 = *a2, s3 = *a3;

    uint8_t r0 = gf_mul(0x0e, s0) ^ gf_mul(0x0b, s1) ^ gf_mul(0x0d, s2) ^ gf_mul(0x09, s3);
    uint8_t r1 = gf_mul(0x09, s0) ^ gf_mul(0x0e, s1) ^ gf_mul(0x0b, s2) ^ gf_mul(0x0d, s3);
    uint8_t r2 = gf_mul(0x0d, s0) ^ gf_mul(0x09, s1) ^ gf_mul(0x0e, s2) ^ gf_mul(0x0b, s3);
    uint8_t r3 = gf_mul(0x0b, s0) ^ gf_mul(0x0d, s1) ^ gf_mul(0x09, s2) ^ gf_mul(0x0e, s3);

    *a0 = r0; *a1 = r1; *a2 = r2; *a3 = r3;
}

static void inv_mix_columns(state_t s) {
    for (int c = 0; c < 4; ++c) {
        inv_mix_single_column(&s[0][c], &s[1][c], &s[2][c], &s[3][c]);
    }
}

/* Encrypt one 16-byte block */
static void aes128_encrypt_block(const uint8_t in[16], uint8_t out[16], const uint8_t round_keys[176]) {
    state_t s;
    load_state(s, in);

    add_round_key(s, round_keys + 0);

    for (int r = 1; r <= 9; ++r) {
        sub_bytes(s);
        shift_rows(s);
        mix_columns(s);
        add_round_key(s, round_keys + 16 * r);
    }

    sub_bytes(s);
    shift_rows(s);
    add_round_key(s, round_keys + 16 * 10);

    store_state(out, s);
}

/* Decrypt one 16-byte block */
static void aes128_decrypt_block(const uint8_t in[16], uint8_t out[16], const uint8_t round_keys[176]) {
    state_t s;
    load_state(s, in);

    add_round_key(s, round_keys + 16 * 10);

    for (int r = 9; r >= 1; --r) {
        inv_shift_rows(s);
        inv_sub_bytes(s);
        add_round_key(s, round_keys + 16 * r);
        inv_mix_columns(s);
    }

    inv_shift_rows(s);
    inv_sub_bytes(s);
    add_round_key(s, round_keys + 0);

    store_state(out, s);
}

/* PKCS#7 padding helpers */
static uint8_t *pkcs7_pad(const uint8_t *in, size_t len, size_t *out_len) {
    size_t pad = AES_BLOCK_BYTES - (len % AES_BLOCK_BYTES);
    if (pad == 0) pad = AES_BLOCK_BYTES;
    *out_len = len + pad;
    uint8_t *buf = (uint8_t *)malloc(*out_len);
    if (!buf) return NULL;
    memcpy(buf, in, len);
    memset(buf + len, (int)pad, pad);
    return buf;
}

static bool pkcs7_unpad(uint8_t *buf, size_t *len_io) {
    if (*len_io == 0 || (*len_io % AES_BLOCK_BYTES) != 0) return false;
    uint8_t pad = buf[*len_io - 1];
    if (pad == 0 || pad > AES_BLOCK_BYTES) return false;
    if (pad > *len_io) return false;
    for (size_t i = 0; i < pad; ++i) {
        if (buf[*len_io - 1 - i] != pad) return false;
    }
    *len_io -= pad;
    return true;
}

static uint8_t *aes128_ecb_encrypt(const uint8_t *plaintext, size_t len, const uint8_t round_keys[176], size_t *out_len) {
    size_t padded_len;
    uint8_t *padded = pkcs7_pad(plaintext, len, &padded_len);
    if (!padded) return NULL;

    uint8_t *ct = (uint8_t *)malloc(padded_len);
    if (!ct) { free(padded); return NULL; }

    for (size_t i = 0; i < padded_len; i += AES_BLOCK_BYTES) {
        aes128_encrypt_block(padded + i, ct + i, round_keys);
    }

    free(padded);
    *out_len = padded_len;
    return ct;
}

/* Utility: hex */
static int hexval(int c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

static void strip_spaces(char *s) {
    size_t w = 0, n = strlen(s);
    for (size_t r = 0; r < n; ++r) {
        if (!isspace((unsigned char)s[r])) s[w++] = s[r];
    }
    s[w] = '\0';
}

static bool hex_to_bytes_strict(const char *hex, uint8_t **out, size_t *out_len) {
    size_t n = strlen(hex);
    if (n % 2 != 0) return false;
    size_t bytes = n / 2;
    uint8_t *buf = (uint8_t *)malloc(bytes);
    if (!buf) return false;

    for (size_t i = 0; i < bytes; ++i) {
        int hi = hexval(hex[2 * i]);
        int lo = hexval(hex[2 * i + 1]);
        if (hi < 0 || lo < 0) { free(buf); return false; }
        buf[i] = (uint8_t)((hi << 4) | lo);
    }
    *out = buf;
    *out_len = bytes;
    return true;
}

static bool hex_to_bytes_fixed(const char *hex, uint8_t *out, size_t want_len) {
    size_t n = strlen(hex);
    if (n != 2 * want_len) return false;
    for (size_t i = 0; i < want_len; ++i) {
        int hi = hexval(hex[2 * i]);
        int lo = hexval(hex[2 * i + 1]);
        if (hi < 0 || lo < 0) return false;
        out[i] = (uint8_t)((hi << 4) | lo);
    }
    return true;
}

static void bytes_to_hex(const uint8_t *in, size_t len, char *out) {
    static const char *hex = "0123456789abcdef";
    for (size_t i = 0; i < len; ++i) {
        out[2 * i]     = hex[in[i] >> 4];
        out[2 * i + 1] = hex[in[i] & 0x0F];
    }
    out[2 * len] = '\0';
}

/* Read a line from stdin into buf, strip newline. Returns false on EOF. */
static bool read_line(char *buf, size_t cap) {
    if (!fgets(buf, (int)cap, stdin)) return false;
    size_t n = strlen(buf);
    if (n && buf[n - 1] == '\n') buf[n - 1] = '\0';
    return true;
}

int main(void) {
    aes_init_sboxes();

    printf("AES-128 ECB interactive\n");
    printf("Choose mode\n");
    printf("1 for single 16-byte block test\n");
    printf("2 for multi-block ECB with PKCS#7 padding\n");
    printf("> ");
    char choice[8];
    if (!read_line(choice, sizeof(choice))) return 0;

    if (choice[0] == '1') {
        /* Single block test */
        const char *tv_pt  = "3243f6a8885a308d313198a2e0370734";
        const char *tv_key = "2b7e151628aed2a6abf7158809cf4f3c";
        const char *tv_ct  = "3925841d02dc09fbdc118597196a0b32";

        char pt_hex[256], key_hex[256];
        printf("Enter plaintext hex (32 hex chars)\n");
        printf("> ");
        if (!read_line(pt_hex, sizeof(pt_hex))) return 0;
        strip_spaces(pt_hex);

        printf("Enter key hex (32 hex chars)\n");
        printf("> ");
        if (!read_line(key_hex, sizeof(key_hex))) return 0;
        strip_spaces(key_hex);

        uint8_t pt[16], key[16], rk[176], ct[16];
        if (!hex_to_bytes_fixed(pt_hex, pt, 16)) {
            fprintf(stderr, "Error: plaintext must be exactly 32 hex characters.\n");
            return 1;
        }
        if (!hex_to_bytes_fixed(key_hex, key, 16)) {
            fprintf(stderr, "Error: key must be exactly 32 hex characters.\n");
            return 1;
        }

        aes128_expand_key(key, rk);
        aes128_encrypt_block(pt, ct, rk);

        char ct_hex[33];
        bytes_to_hex(ct, 16, ct_hex);
        printf("Computed ciphertext hex\n%s\n", ct_hex);

        bool inputs_match_tv = (strcmp(pt_hex, tv_pt) == 0) && (strcmp(key_hex, tv_key) == 0);
        if (inputs_match_tv) {
            if (strcmp(ct_hex, tv_ct) == 0) {
                printf("Verification against standard test vector: PASS\n");
            } else {
                printf("Verification against standard test vector: FAIL\n");
                printf("Expected\n%s\n", tv_ct);
            }

            /* Optional sanity decrypt */
            uint8_t dec[16];
            aes128_decrypt_block(ct, dec, rk);
            if (memcmp(dec, pt, 16) == 0) {
                printf("Decryption recovered the plaintext.\n");
            } else {
                printf("Decryption did not recover the plaintext.\n");
            }
        } else {
            printf("Note: your inputs are not the standard test vector\n");
            printf("Standard PT\n%s\n", tv_pt);
            printf("Standard KEY\n%s\n", tv_key);
            printf("Standard CT\n%s\n", tv_ct);
        }

    } else if (choice[0] == '2') {
        /* Multi-block ECB with padding */
        char pt_hex_in[1 << 14];   /* 16 KB of hex input capacity */
        char key_hex[256];

        printf("Enter plaintext hex of any length. No spaces. Even number of hex chars required.\n");
        printf("> ");
        if (!read_line(pt_hex_in, sizeof(pt_hex_in))) return 0;
        strip_spaces(pt_hex_in);

        printf("Enter key hex (32 hex chars)\n");
        printf("> ");
        if (!read_line(key_hex, sizeof(key_hex))) return 0;
        strip_spaces(key_hex);

        uint8_t key[16], rk[176];
        if (!hex_to_bytes_fixed(key_hex, key, 16)) {
            fprintf(stderr, "Error: key must be exactly 32 hex characters.\n");
            return 1;
        }
        aes128_expand_key(key, rk);

        uint8_t *pt = NULL;
        size_t pt_len = 0;
        if (!hex_to_bytes_strict(pt_hex_in, &pt, &pt_len)) {
            fprintf(stderr, "Error: plaintext hex must have even length and valid hex chars.\n");
            return 1;
        }

        size_t ct_len = 0;
        uint8_t *ct = aes128_ecb_encrypt(pt, pt_len, rk, &ct_len);
        free(pt);
        if (!ct) {
            fprintf(stderr, "ECB encrypt failed.\n");
            return 1;
        }

        char *ct_hex = (char *)malloc(ct_len * 2 + 1);
        if (!ct_hex) { free(ct); return 1; }
        bytes_to_hex(ct, ct_len, ct_hex);
        free(ct);

        printf("ECB ciphertext hex\n%s\n", ct_hex);
        free(ct_hex);

    } else {
        printf("Unknown choice. Exiting.\n");
    }

    return 0;
}

// test: Plaintext (hex): 3243f6a8885a308d313198a2e0370734

// Key (hex): 2b7e151628aed2a6abf7158809cf4f3c

// Verify that your ciphertext matches the standard AES test vector:

// Ciphertext (hex): 3925841d02dc09fbdc118597196a0b32