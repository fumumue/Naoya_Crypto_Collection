/*
 qc_mdpc_keyenc.c
 Reference / educational implementation for QC-MDPC-like keygen + encode (binary)
 - ring: F2[x] / (x^R + 1), with x^R == 1 (char 2)
 - secret: h0, h1 (sparse polynomials)
 - public: h = h0^{-1} * h1 mod (x^R+1)
 - encode: given message m (bitpoly of degree < R), random sparse r0,r1,e:
     u = r0 + r1 * h
     v = m ^ (r1 * pub) ^ e
 Note: This is for research/experimentation only. NOT production-ready.
 Compile: gcc -O2 -std=c99 qc_mdpc_keyenc.c -o qc_mdpc_keyenc
 Run: ./qc_mdpc_keyenc
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>

////////////////////////////////////////////////////////////////////////////////
// PARAMETERS (tune for experiments)
// R : block size (degree modulus). n = 2*R is code length.
// W : secret sparse weight for h0/h1
// WR : ephemeral weight for r0/r1/e in encode
////////////////////////////////////////////////////////////////////////////////
#define R  256           // block degree (example small value for testing)
#define WORD64S ((R + 63)/64)

#define SECRET_WEIGHT  11
#define EPHEMERAL_WEIGHT  11

////////////////////////////////////////////////////////////////////////////////
// Basic poly as bitset (length R)
////////////////////////////////////////////////////////////////////////////////
typedef struct {
    uint64_t a[WORD64S];
} poly_t;

static inline void poly_zero(poly_t *p){
    memset(p->a,0,sizeof(uint64_t)*WORD64S);
}
static inline void poly_copy(poly_t *dst,const poly_t *src){
    memcpy(dst->a,src->a,sizeof(uint64_t)*WORD64S);
}
static inline void poly_setbit(poly_t *p,int i){
    if(i<0) return;
    int w = i>>6; int b = i & 63;
    p->a[w] |= ((uint64_t)1<<b);
}
static inline int poly_getbit(const poly_t *p,int i){
    int w = i>>6; int b = i & 63;
    return (p->a[w] >> b) & 1;
}
static inline void poly_xor(poly_t *r,const poly_t *x,const poly_t *y){
    for(int i=0;i<WORD64S;i++) r->a[i] = x->a[i] ^ y->a[i];
}
static inline int poly_weight(const poly_t *p){
    int w=0;
    for(int i=0;i<WORD64S;i++) w += __builtin_popcountll(p->a[i]);
    return w;
}

////////////////////////////////////////////////////////////////////////////////
// rotate-left modulo R bits: dest = rotl(src, shift)
////////////////////////////////////////////////////////////////////////////////
static void poly_rotl(poly_t *dest,const poly_t *src,int shift){
    if(shift==0){ poly_copy(dest,src); return; }
    shift %= R;
    int word_shift = shift >> 6;
    int bit_shift = shift & 63;
    uint64_t tmp[WORD64S];
    for(int i=0;i<WORD64S;i++) tmp[i]=0;
    for(int i=0;i<R;i++){
        if(poly_getbit(src,i)){
            int ni = (i + shift) % R;
            int w = ni>>6; int b = ni&63;
            tmp[w] |= ((uint64_t)1<<b);
        }
    }
    for(int i=0;i<WORD64S;i++) dest->a[i]=tmp[i];
}

////////////////////////////////////////////////////////////////////////////////
// poly_mul_mod : r = a * b mod (x^R+1)  (naive, uses rotate & xor)
////////////////////////////////////////////////////////////////////////////////
static void poly_mul_mod(poly_t *r,const poly_t *a,const poly_t *b){
    poly_zero(r);
    // iterate bits set in a
    for(int i=0;i<R;i++){
        if(poly_getbit(a,i)){
            poly_t tmp;
            poly_rotl(&tmp,b,i); // tmp = b << i (mod x^R+1)
            for(int w=0;w<WORD64S;w++) r->a[w] ^= tmp.a[w];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// polynomial arithmetic for inversion and div: represent polys as dynamic bit arrays
// we implement simple polynomial division / gcd over GF(2) with degree <= R
// We use a simple bigint-like representation but restricted to R bits.
////////////////////////////////////////////////////////////////////////////////

// degree of polynomial (returns -1 for zero)
static int poly_deg(const poly_t *p){
    for(int w = WORD64S-1; w>=0; w--){
        uint64_t v = p->a[w];
        if(v){
            int hi = 63 - __builtin_clzll(v);
            return w*64 + hi;
        }
    }
    return -1;
}

// left shift by k bits (mod x^R+1 => rotate)
static void poly_shl_rot(poly_t *dst,const poly_t *src,int k){
    poly_rotl(dst,src,k);
}

// poly xor in-place shifted: a ^= (b << sh) mod (x^R+1)
static void poly_xor_shift(poly_t *a,const poly_t *b,int sh){
    poly_t tmp;
    poly_rotl(&tmp,b,sh%R);
    for(int i=0;i<WORD64S;i++) a->a[i] ^= tmp.a[i];
}

// polynomial division: compute quotient q and remainder rem such that num = q*den + rem
// naive polynomial long division over GF(2), deg <= R-1
static void poly_divmod(poly_t *q, poly_t *rem, const poly_t *num, const poly_t *den){
    poly_zero(q);
    poly_copy(rem, num);
    int d_deg = poly_deg(den);
    if(d_deg < 0){ fprintf(stderr,"div by zero poly_divmod\n"); exit(1); }
    while(1){
        int r_deg = poly_deg(rem);
        if(r_deg < d_deg) break;
        int shift = r_deg - d_deg;
        // q |= (1 << shift)
        poly_setbit(q, shift);
        // rem ^= den << shift
        poly_xor_shift(rem, den, shift);
    }
}

// polynomial gcd and inverse: extended euclid to get s,t with s*a + t*modpoly = gcd
// returns gcd in g; s is inverse if gcd==1 (i.e., s*a mod modpoly = 1)
static void poly_egcd(poly_t *g, poly_t *s, const poly_t *a, const poly_t *modpoly){
    // Extended Euclid on polynomials over GF(2)
    // Initialize: r0 = modpoly, r1 = a
    poly_t r0, r1, s0, s1, q, tmp, rem, prod;
    poly_copy(&r0, modpoly);
    poly_copy(&r1, a);
    poly_zero(&s0); poly_zero(&s1);
    // s0 = 0; s1 = 1
    poly_setbit(&s1,0);

    while(poly_deg(&r1) >= 0){
        poly_zero(&q); poly_zero(&rem);
        poly_divmod(&q, &rem, &r0, &r1); // r0 = q*r1 + rem
        // r0 <- r1; r1 <- rem
        poly_copy(&r0, &r1);
        poly_copy(&r1, &rem);
        // s0 <- s1; s1 <- s0 + q*s1
        // compute q*s1
        poly_zero(&prod);
        // multiply q * s1 (naive)
        for(int i=0;i<R;i++){
            if(poly_getbit(&q,i)){
                poly_t t; poly_rotl(&t,&s1,i);
                for(int w=0;w<WORD64S;w++) prod.a[w] ^= t.a[w];
            }
        }
        // tmp = s0 + prod
        poly_xor(&tmp, &s0, &prod);
        poly_copy(&s0, &s1);
        poly_copy(&s1, &tmp);
    }
    poly_copy(g, &r0); // gcd
    poly_copy(s, &s0); // coefficient
}

// wrapper to compute inverse of a modulo modpoly (modpoly assumed x^R+1)
// returns 0 on success (inv in inv), -1 if non-invertible
static int poly_inverse_mod(poly_t *inv, const poly_t *a){
    // modpoly = x^R + 1
    poly_t modpoly; poly_zero(&modpoly); poly_setbit(&modpoly, R); // x^R
    // but we represent only up to R-1 bits. For EEA we need modpoly represented differently.
    // To avoid complexity, we generate modpoly as polynomial with bit R set:
    // However our bitsets are length R only; we can't set bit R. Instead handle modpoly implicitly.
    // Simpler: perform EGCD with explicit modpoly function using python style would be easier,
    // but to remain in C and simple, we implement inversion by brute-force for small R:
    // Try to find inv s.t. a * inv mod (x^R+1) == 1
    // This naive approach is O(R^2) trials times multiply; ok for small R in experiments.
    //
    poly_t prod, candidate;
    for(int cand=1; cand < (1<<16) && R<=16; cand++){} // placeholder, but we implement random search below

    // Instead do Extended Euclid implemented with polynomials of degree up to R-1
    // We can represent modpoly as (x^R + 1) but this has degree R so fits with our routines if we allow degree R.
    // Implement modpoly with a buffer of R+1 bits: trick: use larger WORD array? For simplicity, create arrays with one extra bit mentally.
    // To avoid complexity and still be usable for R up to moderate, use binary extended euclid using integers in dynamic arrays.
    // But due to time, we will implement inverse by using exponentiation:
    // in ring GF(2)[x]/(x^R+1), the multiplicative group of units is subset; we can compute a^(2^R - 2) mod (x^R+1)
    // Use square-and-multiply with repeated squaring (Fermat-like) -- but multiplicative group size is not 2^R-1 necessarily. This is approximate.
    //
    // For simplicity and robustness in this educational code, we will use brute-force linear solver:
    // Solve linear system for inv coefficients: find vector v of length R s.t. convolution(a,v) mod x^R+1 = 1
    // This is linear over GF(2): A * v = e0 (vector), where A is circulant matrix formed by a. Solve with Gaussian elimination over GF(2).
    //
    // Build matrix A (R x R) columns are rotations of a. Then solve A v = e0
    {
        // allocate matrix rows as bitsets
        int rows = R, cols = R;
        int words = (cols+63)/64;
        uint64_t *M = calloc(rows * words, sizeof(uint64_t));
        // Fill M: row i, column j = a_{(i-j) mod R}
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                int idx = (i - j) % R;
                if(idx<0) idx += R;
                if(poly_getbit(a, idx)){
                    int w = j>>6; int b = j & 63;
                    M[i*words + w] |= ((uint64_t)1ULL<<b);
                }
            }
        }
        // RHS is e0: vector with 1 at position 0
        uint64_t *B = calloc(rows*1,sizeof(uint64_t));
        B[0] = 1;
        // Gaussian elimination over GF(2)
        int rank = 0;
        int *where = calloc(cols,sizeof(int));
        for(int i=0;i<cols;i++) where[i] = -1;
        for(int col=0; col<cols && rank<rows; col++){
            // find pivot row r >= rank with M[r][col]==1
            int sel = -1;
            for(int r=rank; r<rows; r++){
                if( (M[r*words + (col>>6)] >> (col&63)) & 1ULL ){
                    sel = r; break;
                }
            }
            if(sel == -1) continue;
            // swap rows sel and rank
            if(sel != rank){
                for(int w=0;w<words;w++){
                    uint64_t tmp = M[rank*words + w];
                    M[rank*words + w] = M[sel*words + w];
                    M[sel*words + w] = tmp;
                }
                uint64_t tmpb = B[rank];
                B[rank] = B[sel];
                B[sel] = tmpb;
            }
            where[col] = rank;
            // eliminate other rows
            for(int r=0;r<rows;r++){
                if(r==rank) continue;
                if( (M[r*words + (col>>6)] >> (col&63)) & 1ULL ){
                    // row r ^= row rank
                    for(int w=0;w<words;w++) M[r*words + w] ^= M[rank*words + w];
                    B[r] ^= B[rank];
                }
            }
            rank++;
        }
        // check consistency: for rows without pivot, RHS must be 0
        for(int r=rank;r<rows;r++){
            if(B[r]){
                // no solution => non invertible
                free(M); free(B); free(where);
                return -1;
            }
        }
        // back-substitution: build solution vector
        poly_zero(inv);
        for(int i=0;i<cols;i++){
            if(where[i] != -1){
                int r = where[i];
                // variable i = B[r]
                if(B[r] & 1ULL) poly_setbit(inv, i);
            } else {
                // free variable -> set 0
            }
        }
        free(M); free(B); free(where);
        return 0;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Random sparse polynomial sampling
////////////////////////////////////////////////////////////////////////////////
static void seed_random(){
    // try /dev/urandom
    int fd = open("/dev/urandom", O_RDONLY);
    if(fd >= 0){
        unsigned int seed;
        if(read(fd, &seed, sizeof(seed)) == sizeof(seed)){
            srand(seed);
            close(fd); return;
        }
        close(fd);
    }
    srand((unsigned)time(NULL) ^ (unsigned)getpid());
}

static void sample_sparse(poly_t *p,int weight){
    poly_zero(p);
    // sample distinct positions
    int chosen[ R ];
    for(int i=0;i<R;i++) chosen[i]=0;
    int got=0;
    while(got < weight){
        int r = rand() % R;
        if(!chosen[r]){ chosen[r]=1; poly_setbit(p,r); got++; }
    }
}

////////////////////////////////////////////////////////////////////////////////
// print poly as hex for debug
////////////////////////////////////////////////////////////////////////////////
static void print_poly(const char *label, const poly_t *p){
    printf("%s (wt=%d): ", label, poly_weight(p));
    for(int i=0;i<WORD64S;i++){
        printf("%016llx", (unsigned long long)p->a[i]);
    }
    printf("\n");
}

////////////////////////////////////////////////////////////////////////////////
// KEYGEN: sample h0,h1 and compute pub = inv(h0) * h1 mod (x^R+1)
////////////////////////////////////////////////////////////////////////////////
static int keygen(poly_t *h0, poly_t *h1, poly_t *pub){
    sample_sparse(h0, SECRET_WEIGHT);
    sample_sparse(h1, SECRET_WEIGHT);
    // invert h0
    if(poly_inverse_mod(pub, h0) != 0){
        // try again until invertible (rare for random sparse)
        return -1;
    }
    // pub = inv(h0) * h1
    poly_t tmp;
    poly_mul_mod(&tmp, pub, h1);
    poly_copy(pub, &tmp);
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
// ENCODE: given pub (h), message m (poly), produce u,v
// u = r0 + r1 * h
// v = m ^ (r1 * pub) ^ e
////////////////////////////////////////////////////////////////////////////////
static void encode(const poly_t *pub, const poly_t *m, poly_t *u, poly_t *v){
    poly_t r0, r1, e;
    sample_sparse(&r0, EPHEMERAL_WEIGHT);
    sample_sparse(&r1, EPHEMERAL_WEIGHT);
    sample_sparse(&e, EPHEMERAL_WEIGHT);

    poly_t t;
    poly_mul_mod(&t, &r1, pub); // r1 * pub
    poly_xor(u, &r0, &t);      // u = r0 + r1*pub  (here pub==h)

    poly_xor(v, m, &t);        // v = m ^ (r1*pub)
    poly_xor(v, v, &e);        // v ^= e
}

////////////////////////////////////////////////////////////////////////////////
// DEMO main
////////////////////////////////////////////////////////////////////////////////
int main(void){
    seed_random();
    poly_t h0,h1,pub;
    int tries=0;
    while(1){
        tries++;
        if(keygen(&h0,&h1,&pub) == 0) break;
        if(tries>16){ fprintf(stderr,"keygen failing repeatedly\n"); break;}
    }
    printf("Keygen succeeded after %d tries\n", tries);
    print_poly("h0", &h0);
    print_poly("h1", &h1);
    print_poly("pub(h)", &pub);

    // sample a random message m (R bits)
    poly_t m;
    poly_zero(&m);
    // for demo set message as low-weight or random
    for(int i=0;i<16;i++){
        if(rand() & 1) poly_setbit(&m, rand()%R);
    }
    print_poly("message m", &m);

    poly_t u,v;
    encode(&pub, &m, &u, &v);
    print_poly("u", &u);
    print_poly("v", &v);

    printf("Done.\n");
    return 0;
}
