#include <stdio.h>
#include <stdlib.h>

#define P 7       // GF(7)
#define MAX_DEG 16
#define MAX_POINTS 32

// 2変数多項式 Q(x,y) の簡易表現: Q(x,y) = sum_{i,j} coef[i][j] * x^i y^j
typedef struct {
    int coef[MAX_DEG][MAX_DEG]; // coef[i][j] = x^i y^j の係数
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

// Q(x,y) の表示
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

// Koetter補間アルゴリズム（最小版）
void koetterInterpolation(int n, int alpha[], int y[], Poly2D *Q) {
    // 初期基底: Q(x,y) = y
    initPoly2D(Q);
    Q->coef[0][1] = 1; // y
    Q->degY = 1;

    // 各点に対して条件を追加していく
    for (int k=0; k<n; k++) {
        int a = alpha[k];
        int b = y[k];

        // Q(a,b) を評価
        int val = 0;
        for (int i=0; i<=Q->degX; i++) {
            for (int j=0; j<=Q->degY; j++) {
                int c = Q->coef[i][j];
                if (c != 0) {
                    // a^i b^j
                    int powa=1, powb=1;
                    for (int ii=0; ii<i; ii++) powa = modp(powa*a);
                    for (int jj=0; jj<j; jj++) powb = modp(powb*b);
                    val = modp(val + c * powa * powb);
                }
            }
        }

        // もし Q(a,b)=0 でなければ、補正を加える
        if (val != 0) {
            // 簡易版: Q(x,y) ← (x - a)*Q(x,y)
            Poly2D newQ;
            initPoly2D(&newQ);
            for (int i=0; i<=Q->degX; i++) {
                for (int j=0; j<=Q->degY; j++) {
                    int c = Q->coef[i][j];
                    if (c != 0) {
                        newQ.coef[i+1][j] = modp(newQ.coef[i+1][j] + c);   // x*Q
                        newQ.coef[i][j]   = modp(newQ.coef[i][j] - a*c);  // -a*Q
                    }
                }
            }
            newQ.degX = Q->degX+1;
            newQ.degY = Q->degY;
            *Q = newQ;
        }
    }
}

int main() {
    // 入力点 (alpha, y)
    int n = 3;
    int alpha[3] = {1, 2, 3};
    int y[3]     = {2, 4, 1};

    Poly2D Q;
    koetterInterpolation(n, alpha, y, &Q);

    printPoly2D(&Q);
    return 0;
}
