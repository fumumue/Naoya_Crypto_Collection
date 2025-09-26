// schur_test.c
// Schur product distinguisher test for codes over GF(P)
// Compile: gcc schur_test.c -o schur_test
// Usage: 生成行列 G を用意し、関数 schur_rank_test を呼ぶ

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define P 257  // 素体の法（例: GF(257)）
static inline int modp(long long x){ int r = (int)(x % P); if(r < 0) r += P; return r; }

// --- 線形代数: ガウス消去法 mod P で行列の階数を返す ---
int rank_modp(int *M, int rows, int cols){
    int r = 0;
    for(int c = 0; c < cols && r < rows; c++){
        // pivot を探す
        int pivot = -1;
        for(int i = r; i < rows; i++){
            if(M[i*cols + c] % P != 0){ pivot = i; break; }
        }
        if(pivot == -1) continue;

        // pivot 行を swap
        if(pivot != r){
            for(int j = 0; j < cols; j++){
                int tmp = M[r*cols + j];
                M[r*cols + j] = M[pivot*cols + j];
                M[pivot*cols + j] = tmp;
            }
        }

        // pivot の逆元
        int a = M[r*cols + c] % P; if(a < 0) a += P;
        // Fermat 逆元
        long long base = a, res = 1;
        int exp = P-2;
        while(exp){
            if(exp & 1) res = (res * base) % P;
            base = (base * base) % P;
            exp >>= 1;
        }
        int inv = (int)res;

        // pivot 行を正規化
        for(int j = c; j < cols; j++)
            M[r*cols + j] = modp((long long)M[r*cols + j] * inv);

        // 他の行を消去
        for(int i = 0; i < rows; i++){
            if(i == r) continue;
            int factor = M[i*cols + c];
            if(factor == 0) continue;
            for(int j = c; j < cols; j++){
                M[i*cols + j] = modp((long long)M[i*cols + j] - (long long)factor * M[r*cols + j]);
            }
        }
        r++;
    }
    return r;
}

// --- Schur product test ---
// G: k x n (row-major) 生成行列
// 戻り値: Schur 積符号の rank
int schur_rank_test(int *G, int k, int n){
    int maxv = k*(k+1)/2; // vstack of all i<=j
    int *S = malloc(sizeof(int)*maxv*n);
    int row = 0;
    for(int i = 0; i < k; i++){
        for(int j = i; j < k; j++){
            for(int c = 0; c < n; c++){
                S[row*n + c] = modp((long long)G[i*n + c] * G[j*n + c]);
            }
            row++;
        }
    }
    int rank = rank_modp(S, row, n);
    free(S);
    return rank;
}

// --- テスト用 main ---
// ランダム G を作って実行する例
int main(void){
    srand((unsigned)time(NULL));
    int k = 32, n = 64;
    int *G = malloc(sizeof(int)*k*n);
    // ランダム生成行列
    for(int i=0;i<k*n;i++) G[i] = rand() % P;

    printf("G =\n");
    for(int i=0;i<k;i++){
        for(int j=0;j<n;j++) printf("%3d ", G[i*n+j]);
        printf("\n");
    }

    int r = schur_rank_test(G, k, n);
    printf("Schur rank = %d (out of max %d)\n", r, n);

    free(G);
    return 0;
}
