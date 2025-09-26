#include <stdio.h>
#include <stdlib.h>

#define P 257
#define MAX_DEG 64

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

int inv(int a) {
    return powmod(a, P-2); // Fermat
}

// Legendre記号
int legendre(int a) {
    if (a == 0) return 0;
    int t = powmod(a, (P-1)/2);
    return (t == P-1) ? -1 : t;
}

// Tonelli–Shanks (平方根計算)
int modsqrt(int a) {
    if (a == 0) return 0;
    if (legendre(a) != 1) return -1; // 解なし

    if (P % 4 == 3) {
        return powmod(a, (P+1)/4);
    }

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
// 二次方程式: a y^2 + b y + c = 0
// =============================
int solve_quadratic(int a, int b, int c, int roots[2]) {
    if (a == 0) {
        if (b == 0) return 0;
        roots[0] = modp(-c * inv(b));
        return 1;
    }
    int disc = modp(b*b - 4*a*c);
    int sqrtD = modsqrt(disc);
    if (sqrtD < 0) return 0; // 解なし
    int inv2a = inv(2*a);
    roots[0] = modp((-b + sqrtD) * inv2a);
    roots[1] = modp((-b - sqrtD) * inv2a);
    if (roots[0] == roots[1]) return 1;
    return 2;
}

// =============================
// 三次方程式: a y^3 + b y^2 + c y + d = 0
// =============================
int solve_cubic(int a, int b, int c, int d, int roots[3]) {
    if (a == 0) return solve_quadratic(b, c, d, roots);

    // 変換: y = x - b/(3a) → x^3 + px + q
    int inva = inv(a);
    long long bb = (long long)b*inva % P;
    long long cc = (long long)c*inva % P;
    long long dd = (long long)d*inva % P;

    long long p = modp(cc - (bb*bb % P)/3);
    long long q = modp((2*bb*bb%P*bb%P)/27 - bb*cc/3 + dd);

    // 判別式 Δ = (q/2)^2 + (p/3)^3
    long long Q = modp(q*q/4);
    long long P3 = modp(p*p*p/27);
    long long D = modp(Q + P3);

    int count=0;

    if (D == 0) {
        // 重解あり
        int u = modp(-q/2);
        int r = modsqrt(u);
        if (r>=0) {
            roots[count++] = modp(r - bb/3);
            roots[count++] = modp(-r - bb/3);
        }
    } else {
        // 一般カルダノ: u^3 + v^3 = -q, uv = -p/3
        // 実装を簡単にするため、総当たり探索で解を探す
        for (int y=0; y<P; y++) {
            int val = modp(((long long)a*y%P*y%P*y)%P + b*y%P*y%P + c*y%P + d);
            if (val==0) roots[count++] = y;
            if (count==3) break;
        }
    }
    return count;
}

// =============================
// デモ
// =============================
int main() {
    // 例: y^2 + 2y + 3 = 0
    int roots[3];
    int n = solve_quadratic(1,2,3,roots);
    printf("二次方程式の解: ");
    for (int i=0;i<n;i++) printf("%d ", roots[i]);
    printf("\n");

    // 例: y^3 + y + 1 = 0
    int m = solve_cubic(1,0,1,1,roots);
    printf("三次方程式の解: ");
    for (int i=0;i<m;i++) printf("%d ", roots[i]);
    printf("\n");

    return 0;
}
