#include <stdio.h>
#include <stdlib.h>

#define P 7       // 体 GF(7) を使う
#define MAX_DEG 16

// 多項式 f(x) 構造体（係数は mod P）
typedef struct {
    int deg;
    int coef[MAX_DEG]; // coef[i] = x^i の係数
} Poly;

// modを正に調整
int modp(int a) {
    int r = a % P;
    return (r < 0) ? r + P : r;
}

// 正規化（最高次項が0なら次数を下げる）
void normalize(Poly *f) {
    while (f->deg > 0 && f->coef[f->deg] == 0) f->deg--;
}

// コピー
void copyPoly(Poly *dst, Poly *src) {
    dst->deg = src->deg;
    for (int i = 0; i <= src->deg; i++)
        dst->coef[i] = src->coef[i];
}

// 多項式の加算
void polyAdd(Poly *res, Poly *a, Poly *b) {
    res->deg = (a->deg > b->deg ? a->deg : b->deg);
    for (int i = 0; i <= res->deg; i++) {
        int ca = (i <= a->deg ? a->coef[i] : 0);
        int cb = (i <= b->deg ? b->coef[i] : 0);
        res->coef[i] = modp(ca + cb);
    }
    normalize(res);
}

// 多項式の減算
void polySub(Poly *res, Poly *a, Poly *b) {
    res->deg = (a->deg > b->deg ? a->deg : b->deg);
    for (int i = 0; i <= res->deg; i++) {
        int ca = (i <= a->deg ? a->coef[i] : 0);
        int cb = (i <= b->deg ? b->coef[i] : 0);
        res->coef[i] = modp(ca - cb);
    }
    normalize(res);
}

// 多項式の乗算
void polyMul(Poly *res, Poly *a, Poly *b) {
    for (int i = 0; i < MAX_DEG; i++) res->coef[i] = 0;
    res->deg = a->deg + b->deg;
    for (int i = 0; i <= a->deg; i++) {
        for (int j = 0; j <= b->deg; j++) {
            res->coef[i+j] = modp(res->coef[i+j] + a->coef[i] * b->coef[j]);
        }
    }
    normalize(res);
}

// 多項式の割り算 f = f mod g
void polyMod(Poly *f, Poly *g) {
    while (f->deg >= g->deg && f->deg > 0) {
        int coef = modp(f->coef[f->deg]); // gはモニックと仮定
        int shift = f->deg - g->deg;
        for (int i = 0; i <= g->deg; i++) {
            f->coef[i + shift] = modp(f->coef[i + shift] - coef * g->coef[i]);
        }
        normalize(f);
    }
}

// GCD(f,g)
void polyGCD(Poly *f, Poly *g, Poly *res) {
    Poly a, b, r;
    copyPoly(&a, f);
    copyPoly(&b, g);
    while (b.deg > 0 || b.coef[0] != 0) {
        copyPoly(&r, &a);
        polyMod(&r, &b);
        copyPoly(&a, &b);
        copyPoly(&b, &r);
    }
    copyPoly(res, &a);
}

// 多項式の表示
void printPoly(Poly *f, char var) {
    for (int i = f->deg; i >= 0; i--) {
        if (f->coef[i] != 0) {
            printf("%d*%c^%d ", f->coef[i], var, i);
            if (i > 0) printf("+ ");
        }
    }
    if (f->deg == 0 && f->coef[0] == 0) printf("0");
    printf("\n");
}

// ======= デモ: Q(x,y) の因数から候補 f(x) を抽出 =======
int main() {
    // 例: Q(x,y) = (y - (x+1)) * (y - (2x+3)) over GF(7)
    // → 候補 f1(x) = x+1, f2(x) = 2x+3

    // f1(x) = x+1
    Poly f1 = {1, {1,1}};
    // f2(x) = 2x+3
    Poly f2 = {1, {3,2}};

    printf("候補多項式（リスト復号の出力）:\n");
    printf("f1(x) = "); printPoly(&f1, 'x');
    printf("f2(x) = "); printPoly(&f2, 'x');

    // 実際には Q(x,y) の因数分解や GCD を使ってこれらを取り出す
    // 今回は簡単化して、あらかじめ構築した形を出力するデモにしている
    return 0;
}
