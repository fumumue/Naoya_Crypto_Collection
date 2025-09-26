#include <stdio.h>

// ---- GF(17) 演算 ----
#define Q 17  // 素体の法

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
