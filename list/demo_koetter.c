#include <stdio.h>
#include <stdlib.h>

#define P 7           // 体 GF(7)
#define MAX_DEG 32    // 最大次数
#define MAX_BASIS 16  // 最大基底数

// 2変数多項式 Q(x,y) = sum coef[i][j] * x^i y^j
typedef struct {
    int coef[MAX_DEG][MAX_DEG];
    int degX, degY;
} Poly2D;

// mod演算
int modp(int a) {
    int r = a % P;
    return (r < 0) ? r + P : r;
}

// 初期化
void initPoly2D(Poly2D *Q) {
    for (int i=0; i<MAX_DEG; i++)
        for (int j=0; j<MAX_DEG; j++)
            Q->coef[i][j] = 0;
    Q->degX = Q->degY = 0;
}

// コピー
void copyPoly2D(Poly2D *dst, Poly2D *src) {
    for (int i=0; i<MAX_DEG; i++)
        for (int j=0; j<MAX_DEG; j++)
            dst->coef[i][j] = src->coef[i][j];
    dst->degX = src->degX;
    dst->degY = src->degY;
}

// 表示
void printPoly2D(Poly2D *Q) {
    printf("Q(x,y) = ");
    int empty = 1;
    for (int i=0; i<=Q->degX; i++) {
        for (int j=0; j<=Q->degY; j++) {
            if (Q->coef[i][j] != 0) {
                printf("%d*x^%d*y^%d + ", Q->coef[i][j], i, j);
                empty = 0;
            }
        }
    }
    if (empty) printf("0");
    printf("\n");
}

// xシフト: Q(x,y) ← (x - a)Q(x,y)
void shiftX(Poly2D *Q, int a) {
    Poly2D newQ;
    initPoly2D(&newQ);
    for (int i=0; i<=Q->degX; i++) {
        for (int j=0; j<=Q->degY; j++) {
            int c = Q->coef[i][j];
            if (c != 0) {
                newQ.coef[i+1][j] = modp(newQ.coef[i+1][j] + c);
                newQ.coef[i][j]   = modp(newQ.coef[i][j] - a*c);
            }
        }
    }
    newQ.degX = Q->degX+1;
    newQ.degY = Q->degY;
    *Q = newQ;
}

// yシフト: Q(x,y) ← (y - b)Q(x,y)
void shiftY(Poly2D *Q, int b) {
    Poly2D newQ;
    initPoly2D(&newQ);
    for (int i=0; i<=Q->degX; i++) {
        for (int j=0; j<=Q->degY; j++) {
            int c = Q->coef[i][j];
            if (c != 0) {
                newQ.coef[i][j+1] = modp(newQ.coef[i][j+1] + c);
                newQ.coef[i][j]   = modp(newQ.coef[i][j] - b*c);
            }
        }
    }
    newQ.degX = Q->degX;
    newQ.degY = Q->degY+1;
    *Q = newQ;
}

// 多点補間 Koetter アルゴリズム
void koetterInterpolation(int n, int alpha[], int beta[], int m, int L, Poly2D *Q) {
    Poly2D basis[MAX_BASIS];

    // 初期基底 B = {1, y, y^2, ..., y^L}
    for (int j=0; j<=L; j++) {
        initPoly2D(&basis[j]);
        basis[j].coef[0][j] = 1;
        basis[j].degY = j;
    }

    // 各点 (alpha, beta) について多重度条件を課す
    for (int k=0; k<n; k++) {
        int a = alpha[k];
        int b = beta[k];
        for (int mult=0; mult<m; mult++) {
            for (int j=0; j<=L; j++) {
                // 評価して 0 でなければ補正（簡易版）
                int val = basis[j].coef[0][j]; // 本来は偏導関数評価
                if (val != 0) {
                    shiftX(&basis[j], a);
                    shiftY(&basis[j], b);
                }
            }
        }
    }

    // 最小次数の多項式を選ぶ
    int best = 0;
    int bestDeg = 1e9;
    for (int j=0; j<=L; j++) {
        int deg = basis[j].degX + (L+1)*basis[j].degY;
        if (deg < bestDeg) {
            bestDeg = deg;
            best = j;
        }
    }
    copyPoly2D(Q, &basis[best]);
}

int main() {
    // 入力点
    int n = 3;
    int alpha[3] = {1, 2, 3};
    int beta[3]  = {2, 4, 1};
    int m = 2; // 多重度
    int L = 2; // 最大y次数

    Poly2D Q;
    koetterInterpolation(n, alpha, beta, m, L, &Q);

    printPoly2D(&Q);
    return 0;
}
