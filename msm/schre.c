// schur_rank.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define MAXK 128   // 最大行数 G の行数
#define MAXN 512   // 最大列数（符号長）
#define MAXS (MAXK*(MAXK+1)/2) // Schur 行数上限

// モジュロ q (素数) を外部で設定して使う
static int Q = 7;

// 正規化 mod Q
static inline int modq(int x) {
    x %= Q;
    if (x < 0) x += Q;
    return x;
}
static inline int gf_mul(int a, int b) {
    return (int)((1LL * (a % Q) * (b % Q)) % Q);
}

// G (k x n) の行 i と行 j の要素ごとの積を S_row に書き込む
// G: flattened row-major, G[r * n + c]
// S_row must have length >= n
void schur_row_product(const int *G, int k, int n, int i, int j, int *S_row) {
    const int *Gi = G + i * n;
    const int *Gj = G + j * n;
    for (int c = 0; c < n; ++c) {
        S_row[c] = gf_mul(Gi[c], Gj[c]);
    }
}

// Gauss (行基本変形) で行列の rank を求める（行数 r, 列数 n）
// 行列は row-major ints, mutated in place
int gf_row_rank(int *A, int rows, int cols) {
    int rank = 0;
    for (int col = 0, row = 0; col < cols && row < rows; ++col) {
        // pivot を見つける（行 row..rows-1 で A[r][col] != 0）
        int sel = -1;
        for (int r = row; r < rows; ++r) {
            if (A[r * cols + col] % Q != 0) { sel = r; break; }
        }
        if (sel == -1) continue;
        // swap sel と row
        if (sel != row) {
            for (int c = col; c < cols; ++c) {
                int tmp = A[sel * cols + c];
                A[sel * cols + c] = A[row * cols + c];
                A[row * cols + c] = tmp;
            }
        }
        // normalize pivot to 1 by multiplying row by inv
        int pivot = modq(A[row * cols + col]);
        // compute inverse of pivot mod Q (Q is small prime)
        int inv = 0;
        for (int t = 1; t < Q; ++t) if ((pivot * t) % Q == 1) { inv = t; break; }
        assert(inv != 0);
        for (int c = col; c < cols; ++c) {
            A[row * cols + c] = modq(A[row * cols + c] * inv);
        }
        // eliminate other rows
        for (int r = 0; r < rows; ++r) {
            if (r == row) continue;
            int factor = A[r * cols + col];
            if (factor == 0) continue;
            for (int c = col; c < cols; ++c) {
                int val = modq(A[r * cols + c] - factor * A[row * cols + c]);
                A[r * cols + c] = val;
            }
        }
        ++row;
        ++rank;
    }
    return rank;
}

// G (k x n) -> S (srows x n) を作る。S の行は g_i ◦ g_j for 0<=i<=j<k
// 戻り値: srows (number of rows actually written). caller が Sbuf を確保すること。
// Sbuf は十分な大きさ (MAXS * n) を想定
int build_schur_matrix(const int *G, int k, int n, int *Sbuf) {
    int idx = 0;
    for (int i = 0; i < k; ++i) {
        for (int j = i; j < k; ++j) {
            schur_row_product(G, k, n, i, j, Sbuf + idx * n);
            ++idx;
        }
    }
    return idx;
}

// 高レベル関数: G->Sを作り rank を返す
int schur_span_rank(const int *G, int k, int n, int q) {
    Q = q; // set modulus
    static int S[MAXS * MAXN]; // static buffer (注意: サイズ制約に注意)
    if (k > MAXK || n > MAXN) {
        fprintf(stderr, "k or n too large (MAXK=%d MAXN=%d)\n", MAXK, MAXN);
        exit(1);
    }
    int srows = build_schur_matrix(G, k, n, S);
    // make a copy because gf_row_rank mutates matrix
    int *M = (int *)malloc(sizeof(int) * srows * n);
    if (!M) { perror("malloc"); exit(1); }
    for (int i = 0; i < srows * n; ++i) M[i] = modq(S[i]);
    int rank = gf_row_rank(M, srows, n);
    free(M);
    return rank;
}

// 小さな動作サンプル
int main(void) {
    // サンプル：G (k x n)
    // 例 q=7, k=3, n=7 (ここは実験用に小さめ)
    int q = 7;
    int k = 3;
    int n = 7;
    // G 行列を適当に埋める（例）
    int G[ MAXK * MAXN ];
    memset(G, 0, sizeof(G));
    // 例: 3 行を代入
    int sample[3][7] = {
        {1,1,1,1,1,1,1},
        {1,3,2,6,4,5,0},
        {1,2,4,1,2,4,0}
    };
    for (int i = 0; i < k; ++i) for (int j = 0; j < n; ++j) G[i*n + j] = sample[i][j] % q;

    int rank = schur_span_rank(G, k, n, q);
    printf("Schur span rank = %d (over GF(%d))\n", rank, q);
    return 0;
}
