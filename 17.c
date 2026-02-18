#include <stdio.h>
#include <stdint.h>

//#define DEG 256
//#define N 17

// ---- GF(17) 演算 ----
#define Q N  // 素体の法

int gf_add(int a, int b) {
    int r = (a + b) % Q;
    if (r < 0) r += Q;
    return r;
}

int gf_mul(int a, int b) {
    long long r = (long long)a * b % Q;
    if (r < 0) r += Q;
    return (int)r;
}

int gf_pow(int a, int e) {
    long long res = 1, base = a % Q;
    while (e > 0) {
        if (e & 1) res = (res * base) % Q;
        base = (base * base) % Q;
        e >>= 1;
    }
    return (int)res;
}

// 逆元 (a^(Q-2) mod Q) — フェルマーの小定理
int gf_inv(int a) {
    if (a == 0) {
        fprintf(stderr, "division by zero in GF(17)!\n");
        return 0;
    }
    return gf_pow(a, Q - 2);
}

int gf_div(int a, int b) {
    return gf_mul(a, gf_inv(b));
}

// ---- ラグランジュ補間 ----
void lagrange_interpolate(const int *xs, const int *ys, int k, int *result) {
    for (int i = 0; i < k; i++) result[i] = 0;

    for (int i = 0; i < k; i++) {
        // 基底多項式 l_i(x) の分子
        int numer[64] = {0}; // k が小さい前提
        numer[0] = 1;
        int numer_deg = 0;
        int denom = 1;

        for (int j = 0; j < k; j++) {
            if (j == i) continue;

            int new_numer[64] = {0};
            for (int d = 0; d <= numer_deg; d++) {
                int c = numer[d];
                if (c == 0) continue;

                // (x - xj) = x + (-xj) だが mod Q
                new_numer[d]   = gf_add(new_numer[d], gf_mul(c, (Q - xs[j]) % Q));
                new_numer[d+1] = gf_add(new_numer[d+1], c);
            }
            for (int d = 0; d <= numer_deg+1; d++) numer[d] = new_numer[d];
            numer_deg++;

            denom = gf_mul(denom, gf_add(xs[i], (Q - xs[j]) % Q)); // xi - xj
        }

        int scale = gf_div(ys[i], denom);

        for (int d = 0; d <= numer_deg; d++) {
            if (numer[d] == 0) continue;
            result[d] = gf_add(result[d], gf_mul(numer[d], scale));
        }
    }
}

/*
int gf_sub(int a, int b) { return (a - b + Q) % Q; }

#define DEG 8  // 多項式の最大次数

// バイナリ多項式の掛け算（mod 2）
void poly_mul(const uint8_t *a, int deg_a,
              const uint8_t *b, int deg_b,
              uint8_t *res) {
    memset(res, 0, (deg_a + deg_b + 1) * sizeof(uint8_t));
    for (int i = 0; i <= deg_a; i++) {
        if (a[i]) {
            for (int j = 0; j <= deg_b; j++) {
                res[i+j] ^= b[j];  // XOR が mod 2 加算
            }
        }
    }
}

// バイナリ多項式の割り算（商を求める）
void poly_div(const uint8_t *num, int deg_num,
              const uint8_t *den, int deg_den,
              uint8_t *quot) {
    uint8_t rem[DEG*2+1];
    memcpy(rem, num, (deg_num+1)*sizeof(uint8_t));
    memset(quot, 0, (deg_num - deg_den + 1) * sizeof(uint8_t));

    for (int i = deg_num; i >= deg_den; i--) {
        if (rem[i]) {
            quot[i - deg_den] = 1;
            for (int j = 0; j <= deg_den; j++) {
                rem[i - deg_den + j] ^= den[j];
            }
        }
    }
}

// 多項式の表示
void poly_print(const uint8_t *p, int deg) {
    for (int i = deg; i >= 0; i--) {
        if (p[i]) printf("1");
        else printf("0");
    }
    printf("\n");
}



// 逆元（フェルマーの小定理 a^(Q-2) mod Q）
int gf_inv_q(int a) {
    int res = 1, base = a % Q, e = Q - 2;
    while (e > 0) {
        if (e & 1) res = gf_mul(res, base);
        base = gf_mul(base, base);
        e >>= 1;
    }
    return res;
}

// 多項式の掛け算 mod Q
void poly_mul_q(const int *a, int deg_a, const int *b, int deg_b, int *res) {
    for (int i = 0; i <= deg_a + deg_b; i++) res[i] = 0;
    for (int i = 0; i <= deg_a; i++)
        for (int j = 0; j <= deg_b; j++)
            res[i + j] = gf_add(res[i + j], gf_mul(a[i], b[j]));
}

// 多項式の割り算（商）mod Q
void poly_div_q(const int *num, int deg_num, const int *den, int deg_den, int *quot) {
    int rem[DEG*2+1];
    for (int i = 0; i <= deg_num; i++) rem[i] = num[i];
    for (int i = 0; i <= deg_num - deg_den; i++) quot[i] = 0;

    for (int i = deg_num; i >= deg_den; i--) {
        if (rem[i] != 0) {
            int factor = gf_div(rem[i], den[deg_den]);
            quot[i - deg_den] = factor;
            for (int j = 0; j <= deg_den; j++)
                rem[i - deg_den + j] = gf_sub(rem[i - deg_den + j], gf_mul(factor, den[j]));
        }
    }
}
*/

/*
int main() {
    // l(x) = x^3 + 3x + 2
    int l[DEG+1] = {2,3,0,1};
    int deg_l = 3;

    // m(x) = 4x^2 + 5
    int m[DEG+1] = {5,0,4};
    int deg_m = 2;

    // スケーリング C_a = 7
    int lprime[DEG+1];
    for (int i = 0; i <= deg_l; i++) lprime[i] = gf_mul(7, l[i]);
    int deg_lprime = deg_l;

    // l'm
    int lm[DEG*2+1];
    poly_mul_q(lprime, deg_lprime, m, deg_m, lm);
    int deg_lm = deg_lprime + deg_m;

    printf("l' * m = ");
    poly_print(lm, deg_lm);

    // 割って m を復元
    int recovered[DEG+1];
    poly_div_q(lm, deg_lm, lprime, deg_lprime, recovered);

    printf("recovered m = ");
    poly_print(recovered, deg_m);

    return 0;
}
*/

/*
int main() {
    // l(x) = x^3 + x + 1
    uint8_t l[DEG+1] = {1,1,0,1};  // 1 + x + x^3
    int deg_l = 3;

    // m(x) = x^2 + 1
    uint8_t m[DEG+1] = {1,0,1};
    int deg_m = 2;

    // スケーリング C_a = 1 (バイナリでは意味ない)
    uint8_t lprime[DEG+1];
    memcpy(lprime, l, (deg_l+1)*sizeof(uint8_t));
    int deg_lprime = deg_l;

    // l'm
    uint8_t lm[DEG*2+1];
    poly_mul(lprime, deg_lprime, m, deg_m, lm);
    int deg_lm = deg_lprime + deg_m;

    printf("l' * m = ");
    poly_print(lm, deg_lm);

    // 割って m を復元
    uint8_t recovered[DEG+1];
    poly_div(lm, deg_lm, lprime, deg_lprime, recovered);

    printf("recovered m = ");
    poly_print(recovered, deg_m);

    return 0;
}


int mltn(int n,int x){
    int t=1;
    for(int i=0;i<n;i++){
    t*=x;
    t%=N;
    }
    return t;
}


// 多項式の代入値
long long int
trace(vec f,int deg, long long int x)
{
    long long int u = 0;
    vec v = (f);

    for (int i = 0; i < deg+1; i++)
    {
        if (v.x[i] > 0){
            u = (u + (v.x[i] * mltn(i, x))) % N;
            //printf("u%d %d\n",u,i);
        }
    }

    return u;
}




// ---- main ----
int main(void) {
    // 元の多項式 f(x) = 3 + 5x + 2x^2
    vec f = {3, 5, 2};
    int deg = 2;

    // 評価点
    int xs[3] = {1, 2, 4};
    int ys[3];

    for (int i = 0; i < 3; i++) {
        ys[i] = trace(f,deg, xs[i]);
        printf("f(%d) = %d\n", xs[i], ys[i]);
    }

    // 補間
    int rec[3];
    lagrange_interpolate(xs, ys, 3, rec);

    printf("\nrecovered polynomial:\n");
    for (int i = 0; i <= deg; i++) {
        printf("x^%d: %d\n", i, rec[i]);
    }

    return 0;
}
*/