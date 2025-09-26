#include <stdio.h>
#include <stdlib.h>

#define P 257        // GF(257)
#define MAX_DEG 64   // 次数上限
#define MAX_BASIS 16 // 基底上限

// 2変数多項式 Q(x,y) = sum coef[i][j] x^i y^j
typedef struct {
    int coef[MAX_DEG][MAX_DEG];
    int degX, degY;
} Poly2D;

// mod 演算
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

// (x - a)Q に置換（簡易版）
void shiftX(Poly2D *Q, int a) {
    Poly2D newQ; initPoly2D(&newQ);
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

// (y - b)Q に置換（簡易版）
void shiftY(Poly2D *Q, int b) {
    Poly2D newQ; initPoly2D(&newQ);
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

// 多点補間（簡略）
void koetterInterpolation(int n, int alpha[], int beta[], int m, int L, Poly2D *Q) {
    Poly2D basis[MAX_BASIS];

    // 初期基底 {1, y, y^2, ..., y^L}
    for (int j=0; j<=L; j++) {
        initPoly2D(&basis[j]);
        basis[j].coef[0][j] = 1;
        basis[j].degY = j;
    }

    // 各点に対して多重度条件を課す（ここでは単純なシフト）
    for (int k=0; k<n; k++) {
        int a = alpha[k];
        int b = beta[k];
        for (int mult=0; mult<m; mult++) {
            for (int j=0; j<=L; j++) {
                if (basis[j].coef[0][j] != 0) {
                    shiftX(&basis[j], a);
                    shiftY(&basis[j], b);
                }
            }
        }
    }

    // 最小重み次数の多項式を選ぶ
    int best = 0, bestDeg = 1e9;
    for (int j=0; j<=L; j++) {
        int deg = basis[j].degX + (L+1)*basis[j].degY;
        if (deg < bestDeg) {
            bestDeg = deg;
            best = j;
        }
    }
    copyPoly2D(Q, &basis[best]);
}

// 簡易因数分解: y を解いて候補 f(x) を抽出
void factorCandidates(Poly2D *Q) {
    // 仮定: Q(x,y) = a(x) + b(x)*y （y の次数は1）
    // → 候補 f(x) = -a(x)/b(x)
    printf("候補多項式の探索...\n");
    int degY = Q->degY;

    if (degY == 0) {
        printf("Q は y を含まないので候補なし\n");
        return;
    }
    if (degY > 1) {
        printf("y の次数が %d > 1 → 実装外（拡張要）\n", degY);
        return;
    }

    // b(x), a(x) を抽出
    int a[MAX_DEG] = {0};
    int b[MAX_DEG] = {0};
    for (int i=0; i<=Q->degX; i++) {
        a[i] = Q->coef[i][0];
        b[i] = Q->coef[i][1];
    }

    printf("候補 f(x) = -a(x)/b(x)\n");
    printf("a(x) = ");
    for (int i=0; i<=Q->degX; i++) if (a[i]) printf("%d*x^%d + ", a[i], i);
    printf("\n");
    printf("b(x) = ");
    for (int i=0; i<=Q->degX; i++) if (b[i]) printf("%d*x^%d + ", b[i], i);
    printf("\n");
}

int main() {
    int alpha[5] = {1, 2, 3, 4, 5};
    int beta[5]  = {2, 4, 1, 7, 3};
    int n = 5;
    int m = 2;  // 多重度
    int L = 1;  // y の次数上限を 1 に制限

    Poly2D Q;
    koetterInterpolation(n, alpha, beta, m, L, &Q);

    printPoly2D(&Q);

    // 候補抽出
    factorCandidates(&Q);

    return 0;
}
