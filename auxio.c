/*
 toy_hqc_auxio.c
 Educational implementation:
  - Auxiliary decoder: inner=bit-repetition majority, outer=Reed-Solomon over GF(2^4) (n=15, k=11)
  - File I/O utilities: save/load keys and ciphertexts in binary
  - Meant for local experiments / learning only (NOT secure)

 Compile:
   gcc -O2 -std=c99 toy_hqc_auxio.c -o toy_hqc_auxio -lm

 Usage:
   ./toy_hqc_auxio

 Notes:
  - This is a simplified concatenated decoder: inner code is repetition over blocks of size REP,
    outer RS operates on GF(16) symbols.
  - Parameters chosen small for clarity.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>

/* ----------------------------
   Parameters (small toy values)
   ---------------------------- */
#define R           256   /* ring degree for polynomial parts (reused conceptually) */
#define REP         3     /* repetition factor for inner code (must be odd for majority) */

#define RS_m        4     /* GF(2^m) with m=4 => GF(16) */
#define RS_q        (1<<RS_m) /* 16 */
#define RS_n        15    /* code length */
#define RS_k        11    /* message length (so t=(n-k)/2=2 errors correctable) */
#define RS_t        ((RS_n - RS_k)/2)

/* ----------------------------
   Simple binary vector representation used for toy HQC
   ---------------------------- */
#define HQC_N  (R)         /* for simplicity */
typedef struct {
    uint8_t v[HQC_N];
} bv_t;

/* ----------------------------
   Simple file formats (binary):
   - key_hqc.bin: [N bytes x][N bytes y][N bytes h][N bytes s]
   - ct_hqc.bin:  [N bytes u][N bytes v]  (for toy)
   ---------------------------- */

/* ----------------------------
   Random helper
   ---------------------------- */
static void seed_random(){
    int fd = open("/dev/urandom", O_RDONLY);
    if(fd>=0){
        unsigned int s;
        if(read(fd,&s,sizeof(s))==sizeof(s)){ srand(s); close(fd); return;}
        close(fd);
    }
    srand((unsigned)time(NULL) ^ (unsigned)getpid());
}

/* ----------------------------
   Basic bv helpers
   ---------------------------- */
static void bv_zero(bv_t *b){ memset(b->v,0,sizeof(b->v)); }
static void bv_random(bv_t *b){
    for(int i=0;i<HQC_N;i++) b->v[i] = rand() & 1;
}
static void bv_sample_weight(bv_t *b,int w){
    bv_zero(b);
    int chosen[HQC_N]; memset(chosen,0,sizeof(chosen));
    int got=0;
    while(got<w){
        int r = rand() % HQC_N;
        if(!chosen[r]){ chosen[r]=1; b->v[r]=1; got++; }
    }
}
static int bv_weight(const bv_t *b){
    int s=0;
    for(int i=0;i<HQC_N;i++) s += b->v[i];
    return s;
}
static void bv_xor(bv_t *r,const bv_t *a,const bv_t *b){
    for(int i=0;i<HQC_N;i++) r->v[i] = a->v[i] ^ b->v[i];
}
static void bv_print_head(const char *label,const bv_t *b,int head){
    printf("%s (wt=%d) first %d bits: ", label, bv_weight(b), head);
    for(int i=0;i<head;i++) printf("%d", b->v[i]);
    printf("...\n");
}

/* ----------------------------
   File I/O utilities
   ---------------------------- */
static int save_hqc_key(const char *fname, const bv_t *x, const bv_t *y, const bv_t *h, const bv_t *s){
    FILE *f = fopen(fname,"wb");
    if(!f) return -1;
    fwrite(x->v,1,HQC_N,f);
    fwrite(y->v,1,HQC_N,f);
    fwrite(h->v,1,HQC_N,f);
    fwrite(s->v,1,HQC_N,f);
    fclose(f);
    return 0;
}
static int load_hqc_key(const char *fname, bv_t *x, bv_t *y, bv_t *h, bv_t *s){
    FILE *f = fopen(fname,"rb");
    if(!f) return -1;
    if(fread(x->v,1,HQC_N,f)!=HQC_N) { fclose(f); return -1; }
    if(fread(y->v,1,HQC_N,f)!=HQC_N) { fclose(f); return -1; }
    if(fread(h->v,1,HQC_N,f)!=HQC_N) { fclose(f); return -1; }
    if(fread(s->v,1,HQC_N,f)!=HQC_N) { fclose(f); return -1; }
    fclose(f);
    return 0;
}
static int save_hqc_ct(const char *fname, const bv_t *u, const bv_t *v){
    FILE *f = fopen(fname,"wb");
    if(!f) return -1;
    fwrite(u->v,1,HQC_N,f);
    fwrite(v->v,1,HQC_N,f);
    fclose(f);
    return 0;
}
static int load_hqc_ct(const char *fname, bv_t *u, bv_t *v){
    FILE *f = fopen(fname,"rb");
    if(!f) return -1;
    if(fread(u->v,1,HQC_N,f)!=HQC_N){ fclose(f); return -1; }
    if(fread(v->v,1,HQC_N,f)!=HQC_N){ fclose(f); return -1; }
    fclose(f); return 0;
}

/* ----------------------------
   Simple toy HQC encode / decode (to produce data for auxiliary decoder)
   We will encode a small message tmsg (length RS_k symbols) into a bit-vector of length HQC_N
   via concatenation:
     - Outer RS encodes message into RS_n symbols (GF(16))
     - Each GF(16) symbol is mapped to REP bits by a trivial mapping (e.g., low 4 bits -> REP copies)
   For simplicity in toy: encode RS symbols into 4 bits, then replicate each bit REP times.
   So total bits approx RS_n * 4 * REP. We fit into HQC_N (choose params accordingly).
   ---------------------------- */

/* For toy mapping we will use:
   - RS symbol (0..15) -> 4 bits
   - For each of 4 bit positions we produce REP copies in sequence
   So per symbol we output 4*REP bits.
   Ensure RS_n * 4 * REP <= HQC_N
*/
#if (RS_n * 4 * REP) > HQC_N
#error "HQC_N too small for chosen RS_n and REP"
#endif

/* ----------------------------
   GF(16) arithmetic (m=4)
   primitive poly: x^4 + x + 1 => 0b1_0011 = 0x13 (x^4 + x + 1)
   generator alpha: primitive element with alpha^4 = alpha + 1
   We'll build log/antilog tables for GF(16)
   ---------------------------- */
static uint8_t gf_exp[RS_q*2]; /* enough space */
static uint8_t gf_log[RS_q];
static void gf_init(){
    /* primitive poly: x^4 + x + 1 = 0b1_0011 = 0x13 */
    uint8_t prim = 0x13;
    uint8_t x = 1;
    for(int i=0;i<RS_q-1;i++){
        gf_exp[i] = x;
        gf_log[x] = i;
        x = x << 1;
        if(x & RS_q) x ^= prim; /* reduce if degree >= m */
    }
    /* extend */
    for(int i=RS_q-1;i<RS_q*2;i++) gf_exp[i] = gf_exp[i - (RS_q-1)];
    gf_log[0] = 0xFF; /* undefined */
}

/* gf add is xor */
static inline uint8_t gf_add(uint8_t a,uint8_t b){ return a ^ b; }
static inline uint8_t gf_sub(uint8_t a,uint8_t b){ return a ^ b; }
static uint8_t gf_mul(uint8_t a,uint8_t b){
    if(a==0 || b==0) return 0;
    int la = gf_log[a], lb = gf_log[b];
    return gf_exp[la + lb];
}
static uint8_t gf_inv(uint8_t a){
    if(a==0) return 0; /* undefined, but handle */
    int la = gf_log[a];
    return gf_exp[(RS_q-1) - la];
}
static uint8_t gf_pow(uint8_t a,int e){
    if(e==0) return 1;
    if(a==0) return 0;
    int la = gf_log[a];
    int idx = (la * e) % (RS_q-1);
    return gf_exp[idx];
}

/* ----------------------------
   Reed-Solomon (n=15,k=11) simple implementations
   - generator: alpha^0,... alpha^{n-1} as evaluation points (primitive)
   - encoding: evaluate message polynomial at alpha^i (systematic not implemented)
   - decoding: syndrome -> Berlekamp-Massey -> Chien -> Forney
   ---------------------------- */

/* polynomial arithmetic over GF(16) represented as arrays (degree descending) */
static void rs_poly_add(uint8_t *r, const uint8_t *a, const uint8_t *b, int len){
    for(int i=0;i<len;i++) r[i] = a[i] ^ b[i];
}

/* evaluate polynomial (coeffs in degree-descending? We'll use degree-ascending for simplicity)
   poly: coeffs[0..deg] (coeffs[0] = constant)
*/
static uint8_t rs_poly_eval(const uint8_t *coeffs,int deg,uint8_t x){
    uint8_t y = 0;
    uint8_t pow = 1;
    for(int i=0;i<=deg;i++){
        if(coeffs[i]) y ^= gf_mul(coeffs[i], pow);
        pow = gf_mul(pow, x);
    }
    return y;
}

/* RS encode: message msg[0..k-1] as coefficients (degree < k), produce codeword c[0..n-1] evaluated at alpha^i
   Here evaluation points are alpha^i for i=0..n-1
*/
static void rs_encode(const uint8_t *msg,uint8_t *cw){
    // msg length RS_k (coeffs), cw length RS_n
    for(int i=0;i<RS_n;i++){
        uint8_t pt = gf_exp[i]; // alpha^i
        cw[i] = rs_poly_eval(msg, RS_k-1, pt);
    }
}

/* compute syndromes S_j = cw(alpha^j) for j=1..2t
   syndromes stored S[0..2t-1]
*/
static void rs_compute_syndromes(const uint8_t *cw,uint8_t *S){
    for(int j=1;j<=2*RS_t;j++){
        uint8_t s = 0;
        for(int i=0;i<RS_n;i++){
            uint8_t pt = gf_exp[(i*j) % (RS_q -1)]; // alpha^{i*j}
            uint8_t term = gf_mul(cw[i], pt);
            s ^= term;
        }
        S[j-1] = s;
    }
}

/* Berlekamp-Massey to find error locator polynomial sigma of degree <= t
   Input S[0..2t-1], output sigma coeffs (degree <= t), and degree L
   Implement standard BM over GF(16)
*/
static void berlekamp_massey(const uint8_t *S,int Slen,uint8_t *sigma,int *Lout){
    // initialize
    uint8_t C[RS_t+1]; uint8_t B[RS_t+1];
    memset(C,0,sizeof(C)); memset(B,0,sizeof(B));
    C[0] = 1; B[0] = 1;
    int L = 0, m = 1;
    uint8_t b = 1;
    for(int n=0;n<Slen;n++){
        // compute discrepancy d
        uint8_t d = 0;
        for(int i=0;i<=L;i++){
            if(C[i] && S[n - i]) d ^= gf_mul(C[i], S[n - i]);
        }
        if(d==0){
            m++;
        } else {
            uint8_t coef = gf_mul(d, gf_inv(b));
            uint8_t T[RS_t+1]; memset(T,0,sizeof(T));
            // T = C - coef * x^m * B
            for(int i=0;i<=RS_t;i++) T[i] = C[i];
            for(int i=0;i<=RS_t;i++){
                if(B[i]){
                    int idx = i + m;
                    if(idx <= RS_t){
                        T[idx] ^= gf_mul(coef, B[i]);
                    }
                }
            }
            if(2*L <= n){
                // B <- C; b <- d; L <- n+1-L; m=1; C<-T
                memcpy(B, C, sizeof(C));
                b = d;
                int newL = n + 1 - L;
                m = 1;
                memcpy(C, T, sizeof(C));
                L = newL;
            } else {
                memcpy(C, T, sizeof(C));
                m++;
            }
        }
    }
    // copy sigma
    for(int i=0;i<=RS_t;i++) sigma[i] = 0;
    for(int i=0;i<=L;i++) sigma[i] = C[i];
    *Lout = L;
}

/* Chien search: find roots of sigma (locator) => error positions
   sigma degree L, test alpha^{-i} maybe; here we find positions 0..n-1 where sigma(alpha^{-i})==0
*/
static int chien_search(const uint8_t *sigma,int L,int *err_pos){
    int found = 0;
    // For position i (0..n-1) corresponds to evaluation at alpha^{-i}
    for(int i=0;i<RS_n;i++){
        uint8_t val = 0;
        uint8_t inv_pt = gf_exp[(RS_q-1 - i) % (RS_q-1)]; // alpha^{-i}
        uint8_t pow = 1;
        for(int j=0;j<=L;j++){
            if(sigma[j]) val ^= gf_mul(sigma[j], pow);
            pow = gf_mul(pow, inv_pt);
        }
        if(val == 0){
            err_pos[found++] = i;
            if(found > RS_t+2) break;
        }
    }
    return found;
}

/* Forney algorithm to compute error values (for non-binary RS).
   Here we do simplified: compute error value at position i as E = Omega(alpha^{-i}) / (sigma'(alpha^{-i})) .
   But we need Omega (error evaluator). For brevity and since RS_t small we will skip computing magnitudes and assume binary-like errors:
   For our toy we only aim to locate error positions (assume values are 1 in mapped bits), since inner code handles bit-level mapping.
*/
static void forney_dummy(const uint8_t *cw,int num_err,int *err_pos,uint8_t *err_vals){
    for(int i=0;i<num_err;i++) err_vals[i]=1; // assume value 1 (toy)
}

/* RS decode high-level: attempt to find error positions. Return number of errors found, positions in err_pos.
   If no errors or correctable, returns number found. If uncorrectable, returns -1.
*/
static int rs_decode(const uint8_t *cw,int *err_pos){
    uint8_t S[2*RS_t];
    rs_compute_syndromes(cw, S);
    int allzero = 1; for(int i=0;i<2*RS_t;i++) if(S[i]) { allzero=0; break; }
    if(allzero) return 0; // no errors
    uint8_t sigma[RS_t+1]; int L;
    berlekamp_massey(S, 2*RS_t, sigma, &L);
    if(L > RS_t) return -1; // too many errors
    int found = chien_search(sigma, L, err_pos);
    if(found == 0) return -1; // failed
    // we do not compute magnitudes (toy)
    return found;
}

/* ----------------------------
   Concatenation helpers:
   - map RS codeword (symbols 0..15) -> bits: for symbol s, produce 4 bits (lsb..msb), each repeated REP times
   - decode: inner repetition majority to get 4 bits per symbol, reconstruct symbol, feed to RS decode to find symbol errors.
   ---------------------------- */

/* pack RS symbols -> bit vector of length RS_n * 4 * REP (<= HQC_N) */
static void rs_symbols_to_bits(const uint8_t *sym, uint8_t *bits_out){
    int idx = 0;
    for(int i=0;i<RS_n;i++){
        uint8_t v = sym[i] & 0xF;
        for(int bit=0;bit<4;bit++){
            int bv = (v >> bit) & 1;
            for(int r=0;r<REP;r++){
                bits_out[idx++] = bv;
            }
        }
    }
}

/* inner decode: majority on each group of REP bits for each RS symbol/bit */
static void inner_decode_majority(const uint8_t *bits_in, uint8_t *sym_out){
    int idx = 0;
    for(int i=0;i<RS_n;i++){
        uint8_t symbol = 0;
        for(int bit=0;bit<4;bit++){
            int on = 0;
            for(int r=0;r<REP;r++){
                if(bits_in[idx++]) on++;
            }
            int maj = (on > (REP/2)) ? 1 : 0;
            symbol |= (maj << bit);
        }
        sym_out[i] = symbol & 0xF;
    }
}

/* Reconstruct bit-vector 'payload' (length RS_n*4*REP) from a bv_t (HQC_N) by taking first segment */
static void bv_to_payload(const bv_t *bv, uint8_t *payload){
    int total = RS_n * 4 * REP;
    for(int i=0;i<total;i++) payload[i] = bv->v[i];
}

/* Re-encode after error correction: apply symbol corrections (we assume error values 1 for simplicity)
   This function modifies bits_out accordingly: flip bits for corrected RS symbols.
*/
static void apply_rs_corrections(uint8_t *bits_out, int *err_pos, int num_err){
    for(int e=0;e<num_err;e++){
        int symidx = err_pos[e];
        // flip all bits corresponding to this symbol (naive)
        int base = symidx * 4 * REP;
        for(int j=0;j<4*REP;j++) bits_out[base + j] ^= 1;
    }
}

/* Convert bits payload back into bv_t (store into first bits) */
static void payload_to_bv(const uint8_t *payload, bv_t *bv){
    int total = RS_n * 4 * REP;
    bv_zero(bv);
    for(int i=0;i<total;i++) bv->v[i] = payload[i];
}

/* High-level auxiliary decode function:
   Input: received bit-vector w (bv_t), which encodes RS_n symbols with repetition.
   Steps:
    - Extract payload bits (first RS_n*4*REP bits)
    - inner majority -> symbols
    - RS decode to locate symbol errors (positions)
    - if errors found <= RS_t, flip bits corresponding to those symbols (apply corrections)
    - inner majority again to produce final symbols -> output recovered message symbols (first RS_k)
*/
static int auxiliary_decode(const bv_t *w, uint8_t *out_msg_symbols){
    uint8_t payload[RS_n * 4 * REP];
    bv_to_payload(w, payload);
    uint8_t sym[RS_n];
    inner_decode_majority(payload, sym);
    // attempt RS decode
    int err_pos[RS_t*2 + 4];
    int found = rs_decode(sym, err_pos);
    if(found < 0){
        return -1; // failure to decode
    }
    if(found > 0){
        // apply corrections on payload bits
        apply_rs_corrections(payload, err_pos, found);
        // redo inner majority
        inner_decode_majority(payload, sym);
    }
    // output first k symbols as message
    for(int i=0;i<RS_k;i++) out_msg_symbols[i] = sym[i];
    return 0;
}

/* ----------------------------
   Simple demo tying together:
   - Generate random RS message (RS_k symbols)
   - RS encode -> symbols cw[RS_n]
   - Map to bits -> embed into bv_t tmsg (first bits)
   - Simulate channel: flip some bits (introduce errors)
   - Run auxiliary_decode to recover message symbols
   - Show result
   ---------------------------- */

static void demo_auxiliary_decoder(){
    printf("=== Demo: auxiliary decoder (inner=rep, outer=RS GF(16) n=15 k=11) ===\n");
    // init GF
    gf_init();
    // sample random message symbols
    uint8_t msg[RS_k];
    for(int i=0;i<RS_k;i++) msg[i] = rand() & 0xF;
    printf("Original msg symbols: ");
    for(int i=0;i<RS_k;i++) printf("%x ", msg[i]);
    printf("\n");
    uint8_t cw[RS_n];
    rs_encode(msg, cw);
    printf("Encoded symbols (cw): ");
    for(int i=0;i<RS_n;i++) printf("%x ", cw[i]);
    printf("\n");
    // map to bits
    uint8_t payload[RS_n * 4 * REP];
    rs_symbols_to_bits(cw, payload);
    // place into bv_t
    bv_t w; bv_zero(&w);
    payload_to_bv(payload, &w);
    // introduce random bit flips (simulate noise)
    int flips = 8; // flip some bits across payload
    for(int i=0;i<flips;i++){
        int pos = rand() % (RS_n * 4 * REP);
        w.v[pos] ^= 1;
    }
    // show some
    bv_print_head("Received w (head)", &w, 64);
    // auxiliary decode
    uint8_t rec_msg[RS_k];
    int res = auxiliary_decode(&w, rec_msg);
    if(res == 0){
        printf("Aux decode succeeded. Recovered symbols: ");
        for(int i=0;i<RS_k;i++) printf("%x ", rec_msg[i]);
        printf("\n");
    } else {
        printf("Aux decode failed.\n");
    }
}

/* ----------------------------
   Demo: file I/O usage for keys and ciphertexts
   ---------------------------- */
static void demo_fileio(){
    printf("\n=== Demo: file I/O for HQC toy keys and ciphertexts ===\n");
    seed_random();
    bv_t x,y,h,s;
    bv_sample_weight(&x, 8);
    bv_sample_weight(&y, 8);
    bv_random(&h);
    // s = x XOR (h AND y) as in toy earlier
    bv_t hy; bv_zero(&hy);
    for(int i=0;i<HQC_N;i++) hy.v[i] = h.v[i] & y.v[i];
    bv_xor(&s, &x, &hy);
    // save
    if(save_hqc_key("key_hqc.bin",&x,&y,&h,&s)==0) printf("Saved key_hqc.bin\n");
    else printf("Failed saving key\n");
    // load into new variables
    bv_t x2,y2,h2,s2;
    if(load_hqc_key("key_hqc.bin",&x2,&y2,&h2,&s2)==0) printf("Loaded key_hqc.bin\n");
    else printf("Failed load key\n");
    bv_print_head("Loaded x", &x2, 32);
    bv_print_head("Loaded y", &y2, 32);
    // create simple ciphertext: u = random, v = random
    bv_t u,v;
    bv_random(&u); bv_random(&v);
    if(save_hqc_ct("ct_hqc.bin",&u,&v)==0) printf("Saved ct_hqc.bin\n");
    else printf("Failed save ct\n");
    bv_t uu,vv;
    if(load_hqc_ct("ct_hqc.bin",&uu,&vv)==0) printf("Loaded ct_hqc.bin\n");
    else printf("Failed load ct\n");
    bv_print_head("Loaded ct u", &uu, 32);
    bv_print_head("Loaded ct v", &vv, 32);
}

/* ----------------------------
   main
   ---------------------------- */
int main(void){
    seed_random();
    demo_auxiliary_decoder();
    demo_fileio();
    printf("\nAll done.\n");
    return 0;
}
