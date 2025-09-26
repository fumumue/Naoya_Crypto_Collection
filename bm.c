#include <stdio.h>
#include <string.h>

#define N 15     // 符号長 (例)
#define K 9      // 次元
#define Q 17     // 体 GF(17)
#define MAXN 64  // 多項式次数の上限

// 有限体 mod Q 演算
int modq(int x) {
    x %= Q;
    if (x < 0) x += Q;
    return x;
}
int inv(int a) { // 逆元
    for (int i = 1; i < Q; i++) if ((a * i) % Q == 1) return i;
    return -1;
}

// Berlekamp–Massey アルゴリズム
int berlekamp_massey(int S[], int L, int Lambda[]) {
    int C[MAXN], B[MAXN];
    memset(C, 0, sizeof(C));
    memset(B, 0, sizeof(B));
    C[0] = 1; B[0] = 1;

    int degC = 0, degB = 0;
    int b = 1, m = 1;

    for (int n = 0; n < L; n++) {
        // discrepancy 計算
        int d = S[n];
        for (int i = 1; i <= degC; i++) {
            d = modq(d + C[i] * S[n - i]);
        }
        d = modq(d);

        if (d != 0) {
            int T[MAXN];
            memcpy(T, C, sizeof(C));
            int coef = modq(d * inv(b));
            for (int i = 0; i <= degB; i++) {
                C[i + m] = modq(C[i + m] - coef * B[i]);
            }
            if (2 * degC <= n) {
                memcpy(B, T, sizeof(B));
                degB = degC;
                b = d;
                m = 1;
                degC = n + 1 - degC;
            } else {
                m++;
            }
        } else {
            m++;
        }
    }

    // C[] が誤り位置多項式 Λ(x)
    int degLambda = 0;
    for (int i = MAXN - 1; i >= 0; i--) {
        if (C[i] != 0) { degLambda = i; break; }
    }
    for (int i = 0; i <= degLambda; i++) Lambda[i] = modq(C[i]);

    return degLambda;
}

// Λ(x) の根を探す（評価で brute force）
int find_roots(int Lambda[], int degLambda, int roots[]) {
    int count = 0;
    for (int x = 0; x < Q; x++) {
        int val = 0;
        for (int i = degLambda; i >= 0; i--) {
            val = modq(val * x + Lambda[i]);
        }
        if (val == 0) roots[count++] = x;
    }
    return count;
}

int main() {
    // 例: シンドローム列 (仮に長さ 12 個 = 2*(n-k) )
    int S[32] = {3,5,2,7,8,4,1,0,0,0,0,0};
    int Lambda[MAXN] = {0};

    int deg = berlekamp_massey(S, 12, Lambda);

    printf("誤り位置多項式 Λ(x) = ");
    for (int i = deg; i >= 0; i--) printf("%d ", Lambda[i]);
    printf("\n");

    int roots[32];
    int r = find_roots(Lambda, deg, roots);
    printf("誤り位置候補 = ");
    for (int i = 0; i < r; i++) printf("%d ", roots[i]);
    printf("\n");

    return 0;
}
