#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define MAXM 161   // 最大行数
#define MAXN 161   // 最大列数

// 掃き出し法で A を RREF にする (GF(2))
int gauss_elim(uint8_t A[MAXM][MAXN], int m, int n, int pivots[MAXN]) {
    int r = 0; // 現在の行
    for (int c = 0; c < n && r < m; c++) {
        // ピボット探す
        int pivot = -1;
        for (int i = r; i < m; i++) {
            if (A[i][c]) { pivot = i; break; }
        }
        if (pivot == -1) continue; // その列は全部0 → ピボット無し

        // 行をスワップ
        if (pivot != r) {
            for (int j = 0; j < n; j++) {
                uint8_t tmp = A[r][j];
                A[r][j] = A[pivot][j];
                A[pivot][j] = tmp;
            }
        }

        // ピボットを記録
        pivots[c] = r;

        // 他の行を消去 (XOR)
        for (int i = 0; i < m; i++) {
            if (i != r && A[i][c]) {
                for (int j = c; j < n; j++) {
                    A[i][j] ^= A[r][j];
                }
            }
        }

        r++;
    }
    
    for(int i=0;i<16;i++){
        for(int j=0;j<49;j++)
        printf("%d,",A[i][j]);
    printf("\n");
    }
    //exit(1);
    
    
    return r; // ランク
}

// カーネルを計算
int kernel(uint8_t A[MAXM][MAXN], int m, int n, uint8_t basis[MAXN][MAXN]) {
    int pivots[MAXN];
    for (int i = 0; i < n; i++) pivots[i] = -1;

    // 掃き出し法
    uint8_t B[MAXM][MAXN];
    memcpy(B, A, sizeof(uint8_t)*MAXM*MAXN);
    int rank = gauss_elim(B, m, n, pivots);

    int num_free = 0;

    // 自由変数ごとに基底ベクトルを構成
    for (int j = 0; j < n; j++) {
        if (pivots[j] == -1) {
            // 自由変数 j を 1 にして解を作る
            for (int k = 0; k < n; k++) basis[num_free][k] = 0;
            basis[num_free][j] = 1;

            // ピボット列を逆算
            for (int c = 0; c < n; c++) {
                int r = pivots[c];
                if (r != -1) {
                    basis[num_free][c] = B[r][j]; // j 列の係数で決まる
                }
            }
            num_free++;
        }
    }
    return num_free; // カーネルの次元
}

#define K_prod 48
#define N_prod 161
#define H_prod (161-48)

// G: k×n の生成行列
// H: (n-k)×n の検査行列を出力
int make_parity_check(uint8_t G[K_prod][N_prod], int k, int n, uint8_t H[H_prod][N_prod]) {
    // G の行空間の直交補を計算
    // つまり G * x^T = 0 を満たすベクトル集合を求める
    uint8_t basis[MAXN][MAXN];
    int dim = kernel(G, k, n, basis);
    printf("dim=%d\n",dim);
    //exit(1);
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < n; j++) {
            H[i][j] = basis[i][j];
        }
    }
    return dim; // H の行数
}

/*
// ====== 動作テスト ======
int main(void) {
    // 例: A = [[1 0 1 1],
    //          [0 1 1 0]]
    int m = 2, n = 4;
    uint8_t A[MAXM][MAXN] = {0};
    A[0][0] = 1; A[0][2] = 1; A[0][3] = 1;
    A[1][1] = 1; A[1][2] = 1;

    uint8_t basis[MAXN][MAXN];
    int dim = kernel(A, m, n, basis);

    printf("カーネルの次元: %d\n", dim);
    for (int i = 0; i < dim; i++) {
        printf("basis %d: ", i);
        for (int j = 0; j < n; j++) printf("%d", basis[i][j]);
        printf("\n");
    }
    return 0;
}
*/