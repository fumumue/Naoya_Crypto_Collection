#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#define N 32  // 行列サイズ（32x32を例に）

// 1行を32ビット整数でまとめる
typedef union {
    uint32_t u32;   // まとめて処理
    uint8_t b[4];   // バイト単位でアクセスも可能
} Row;

// GF(2)上で正則かどうか判定
int is_invertible(Row *A, int n) {
    int i, j, k;
    for (i = 0; i < n; i++) {
        // ピボットを探す
        int pivot = -1;
        for (j = i; j < n; j++) {
            if ((A[j].u32 >> (n - 1 - i)) & 1) { // i列目をチェック
                pivot = j;
                break;
            }
        }
        if (pivot == -1) return 0; // ピボットなし → 正則でない

        // 行を入れ替え
        if (pivot != i) {
            Row tmp = A[i];
            A[i] = A[pivot];
            A[pivot] = tmp;
        }

        // ピボット行で他の行を消去
        for (j = 0; j < n; j++) {
            if (j != i && ((A[j].u32 >> (n - 1 - i)) & 1)) {
                A[j].u32 ^= A[i].u32; // まとめて消去
            }
        }
    }

    return 1; // すべてのピボットが見つかった → 正則
}

// テスト
int main() {
    Row *A = malloc(N * sizeof(Row));
    unsigned int u=0;
    srand(clock());
    // 適当な32x32バイナリ行列を初期化（例として単位行列）
    for (int i = 0; i < N; i++) {
        for(int j=0;j<32;j++){
        u<<=1;
        u^=rand()%2;
        }

        A[i].u32 = u; //random(); //1u << (N - 1 - i);
    }

    if (is_invertible(A, N)) {
        printf("行列はGF(2)上で正則です\n");
    } else {
        printf("行列はGF(2)上で正則ではありません\n");
    }

    free(A);
    return 0;
}
