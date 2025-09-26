#include <stdio.h>
#include <stdlib.h>

#define P 257
#define N 8   // 符号長（例：8点だけでデモ）

// =============================
// 基本演算
// =============================
int modp(int a) {
    int r = a % P;
    return (r < 0) ? r + P : r;
}

int powmod(int a, int e) {
    long long res = 1, base = a % P;
    while (e) {
        if (e & 1) res = (res * base) % P;
        base = (base * base) % P;
        e >>= 1;
    }
    return (int)res;
}

int inv(int a) { return powmod(a, P-2); }

// Legendre記号
int legendre(int a) {
    if (a == 0) return 0;
    int t = powmod(a, (P-1)/2);
    return (t == P-1) ? -1 : t;
}

// Tonelli–Shanks 法で平方根を計算
int modsqrt(int a) {
    if (a == 0) return 0;
    if (legendre(a) != 1) return -1;
    if (P % 4 == 3) return powmod(a, (P+1)/4);

    int q = P-1, s = 0;
    while ((q & 1) == 0) { q >>= 1; s++; }

    int z = 2;
    while (legendre(z) != -1) z++;

    int m = s;
    long long c = powmod(z, q);
    long long t = powmod(a, q);
    long long r = powmod(a, (q+1)/2);

    while (t != 1) {
        int i=1;
        long long t2 = (t*t)%P;
        while (t2 != 1) {
            t2 = (t2*t2)%P;
            i++;
            if (i==m) return -1;
        }
        long long b = powmod(c, 1<<(m-i-1));
        m = i;
        c = (b*b)%P;
        r = (r*b)%P;
        t = (t*c)%P;
    }
    return (int)r;
}

// =============================
// 二次方程式 a y^2 + b y + c = 0
// =============================
int solve_quadratic(int a, int b, int c, int roots[2]) {
    if (a == 0) {
        if (b == 0) return 0;
        roots[0] = modp(-c * inv(b));
        return 1;
    }
    int disc = modp(b*b - 4*a*c);
    int sqrtD = modsqrt(disc);
    if (sqrtD < 0) return 0;
    int inv2a = inv(2*a);
    roots[0] = modp((-b + sqrtD) * inv2a);
    roots[1] = modp((-b - sqrtD) * inv2a);
    if (roots[0] == roots[1]) return 1;
    return 2;
}

// =============================
// 候補符号語を出力
// f(x) が線形: f(x) = a1*x + a0
// =============================
void print_codeword_linear(int a1, int a0) {
    printf("候補符号語: ");
    for (int i=0;i<N;i++) {
        int alpha = i+1; // 評価点
        int val = modp(a1*alpha + a0);
        printf("%d ", val);
    }
    printf("\n");
}

// =============================
// デモ
// =============================
int main() {
    // 例: Q(x,y) = y^2 + (2x+3)y + (x+1)
    // → A2=1, A1=2x+3, A0=x+1

    int roots[2];
    for (int x=1;x<=3;x++) { // 評価点 x=1,2,3 でデモ
        int A2=1;
        int A1=modp(2*x+3);
        int A0=modp(x+1);
        int n = solve_quadratic(A2,A1,A0,roots);

        printf("x=%d の候補 y: ", x);
        for (int i=0;i<n;i++) printf("%d ", roots[i]);
        printf("\n");
    }

    // 例: 単純な線形候補 f(x) = 5x+7 の符号語
    print_codeword_linear(5,7);

    return 0;
}
