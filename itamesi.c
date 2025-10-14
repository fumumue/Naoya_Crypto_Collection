// product_decoder.c
// compile: gcc -O3 -std=c11 -o product_decoder product_decoder.c
// Note: fill in your LUT functions (golay_lut_decode, hamming_lut_decode or use direct LUT arrays)

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "ker.c"
#define N_ROWS 23   // n2 (rows, length of columns = n_rows)
#define N_COLS 7    // n1 (cols, length of rows = n_cols)


/* Hamming(7,4) 用：シンドローム(3bit) -> error bit index mask (1<<pos) or 0 if none.
   返り値は 7-bit mask (LSB = bit0) 表示で、複数ビット訂正ないなら単一ビットだけ。 */
uint8_t hamming_lut_lookup(int syndrome3) {
    // ここをあなたの table に置き換えてください
    // syndrome (3bit整数: 0..7) → error vector (7bit整数)
static const uint8_t humm_lut[8] = {
    0b0000000, // 0: no error
    0b0000001, // 1: error in bit0
    0b0000010, // 2: error in bit1
    0b0000100, // 3: error in bit2
    0b0001000, // 4: error in bit3
    0b0010000, // 5: error in bit4
    0b0100000, // 6: error in bit5
    0b1000000  // 7: error in bit6
};


    //const int lut[8]={0,0b1000000,0b100000,0b1000,0b10000000,0b10000,0b100,0b10};
    // 例（標準ハミング列対応）:
    //static const uint8_t lut[8] = {0x00, 0x01, 0x02, 0x08, 0x04, 0x10, 0x20, 0x40};
    //return humm_lut[syndrome3 & 7];
    return humm_lut[syndrome3];
}


/* Golay(23,12) 用：シンドローム(11bit) -> 23-bit error mask を返す（uint32_t下位23ビットを使用）
   あなたが持っている LUT を呼び出せる形にして下さい。 */

   uint32_t golay_lut_lookup(uint16_t syndrome11) {
    // ここをあなたの golay LUT 呼び出しに差し替える
    // ダミー（常に 0 = no error）
    #include "sind.h"
    //(void)syndrome11;
    return syndrome[syndrome11];
}

/* ---------- 検査行列を用いたシンドローム計算（実装例） ----------
   もしあなたが行/列用の H 行列を持っていれば、それを使ってシンドロームを計算してください。
   下はハミング(7,4)用の例（3x7）と Golay(23,12)用の例を想定しています。
*/

/* Example: Hamming (3x7) 棋 */
static const uint8_t H_hamming[3][7] = {
    {1,0,1,0,1,0,1},
    {0,1,1,0,0,1,1},
    {0,0,0,1,1,1,1}
};
static const uint8_t G_hamming[4][7] = {
    {1,0,0,0, 0,1,1},
    {0,1,0,0, 1,0,1},
    {0,0,1,0, 1,1,0},
    {0,0,0,1, 1,1,1}
};


static const uint8_t H_golay[11][23]={
{1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,1 ,0 ,0 ,1 ,0},
{0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,1 ,0 ,0 ,1},
{0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,1 ,1 ,0},
{0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,1 ,1},
{0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,1 ,0 ,0 ,0 ,1 ,1 ,1 ,1},
{0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,1 ,1 ,1 ,0 ,1 ,0 ,1 ,0 ,1},
{0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,1 ,0 ,1 ,1 ,1 ,1 ,0 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,1 ,0 ,1 ,1 ,1 ,1 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,1 ,0 ,1 ,1 ,1 ,1 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,1 ,0 ,1 ,1 ,1 ,1},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,1 ,0 ,0 ,1 ,0 ,1}
};

static const uint8_t G_golay2[12][23]={
{1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0},
{0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0},
{0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0},
{0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0},
{0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1}
};
static const uint8_t G_golay[12][23]={
{1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1},
{0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,1 ,0 ,0 ,1},
{0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,1 ,0 ,0 ,1 ,0 ,1 ,0 ,1},
{0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,1 ,0 ,1 ,1},
{0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,1 ,1 ,0 ,1 ,1 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,1 ,1 ,0 ,1 ,1 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,0 ,0 ,1 ,1 ,0 ,1 ,1},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,1 ,0 ,1 ,1 ,1 ,1 ,0 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,1 ,0 ,1 ,1 ,1 ,1 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,1 ,0 ,1 ,1 ,1 ,1},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1 ,0},
{0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,0 ,1 ,0 ,1 ,1 ,1 ,0 ,0 ,0 ,1 ,1}
};

/* ---------- ユーティリティ: シンドローム計算 ---------- */

/* row: length N_COLS ビット -> ハミングのシンドローム(0..7) */
int compute_row_syndrome(uint8_t row_bits[N_COLS]) {
    int s0 = 0, s1 = 0, s2 = 0;
    for (int j = 0; j < N_COLS; j++) {
        s0 ^= (H_hamming[0][j] & row_bits[j]);
        s1 ^= (H_hamming[1][j] & row_bits[j]);
        s2 ^= (H_hamming[2][j] & row_bits[j]);
    }
    return (s0) | (s1<<1) | (s2<<2);
}

/* col: length N_ROWS ビット -> golay syndrome (11-bit) 
   ここはあなたの Golay 検査行列 H_golay を使って実装するか、
   もし LUT が「受信列そのもの」-> error を返すタイプならそちらを使う */
uint16_t compute_col_syndrome(uint8_t col_bits[N_ROWS]) {
    // ここではダミー：常に 0 を返す。実装するなら H_golay を用いる。
    int syn=0;
    for(int i=0;i<11;i++){
        int s=0;
        syn<<=1;
        for(int j=0;j<N_ROWS;j++)
        s^=col_bits[j]&H_golay[i][j];
        syn^=s;
    }

    return syn;
}

// 検算: G * H^T = 0 ?
void check(uint8_t G[48][161], uint8_t H[][161], int Hr) {
    for(int i=0;i<48;i++){
        for(int j=0;j<Hr;j++){
            int sum=0;
            for(int c=0;c<161;c++) sum ^= (G[i][c] & H[j][c]);
            if(sum!=0){
                printf("NG: G[%d]*H[%d]^T != 0\n",i,j);
                return;
            }
        }
    }
    printf("OK: G*H^T = 0 (mod 2)\n");
}

// Kronecker product of matrices A (r1 x c1) and B (r2 x c2).
void kronecker(const uint8_t A[12][23], int r1, int c1,
               uint8_t B[4][7], int r2, int c2,
               uint8_t *C) {
    // Result size: (r1*r2) x (c1*c2)
    for (int i = 0; i < r1; i++) {
        for (int j = 0; j < c1; j++) {
            for (int u = 0; u < r2; u++) {
                for (int v = 0; v < c2; v++) {
                    C[(i*r2+u)*(c1*c2) + (j*c2+v)] = A[i][j] * B[u][v];
                }
            }
        }
    }
}

// Product code generator matrix from G1 (r1 x c1), G2 (r2 x c2)
// Result C has (r1*r2) x (c1*c2)
void product_code(uint8_t G1[][23], int r1, int c1,
                  uint8_t G2[][7],  int r2, int c2,
                  uint8_t *C) {
    for (int u = 0; u < r2; u++) {           // 各行 g2 を選ぶ
        for (int i = 0; i < r1; i++) {       // G1 の各行
            int row = u*r1 + i;
            for (int j = 0; j < c1; j++) {
                for (int v = 0; v < c2; v++) {
                    // 行列ブロックに展開
                    C[row*(c1*c2) + j*c2 + v] = G1[i][j] & G2[u][v];
                }
            }
        }
    }
}

// Hamming(7,4) のシンドローム計算
int compute_syndrome7(uint8_t received) {
    // parity-check equations
    int s1 = ((received >> 0) ^ (received >> 2) ^ (received >> 4) ^ (received >> 6)) & 1;
    int s2 = ((received >> 1) ^ (received >> 2) ^ (received >> 5) ^ (received >> 6)) & 1;
    int s3 = ((received >> 3) ^ (received >> 4) ^ (received >> 5) ^ (received >> 6)) & 1;

    // 3ビットの整数にまとめる (s1がLSB)
    int s=(s3 << 2) | (s2 << 1) | s1;

    return s;
}

uint8_t decode_hamming74(uint8_t codeword) {
    // syndrome計算 (自分のH行列にあわせて実装する)
    int syn = compute_syndrome7(codeword); // 値は0..7
    
    // LUTからエラーパターンを取得
    uint8_t errvec = hamming_lut_lookup(syn);
    
    // codewordを補正
    return codeword ^ errvec;
}


// Hamming(7,4) 符号化
uint8_t encode_hamming74(uint8_t data) {
    // dataは下位4ビットが情報ビット (d1..d4)
    int d1 = (data >> 0) & 1;
    int d2 = (data >> 1) & 1;
    int d3 = (data >> 2) & 1;
    int d4 = (data >> 3) & 1;

    // パリティ計算
    int p1 = d1 ^ d2 ^ d4;
    int p2 = d1 ^ d3 ^ d4;
    int p3 = d2 ^ d3 ^ d4;

    // コード語のビット並び: p1,p2,d1,p3,d2,d3,d4
    uint8_t codeword = (p1 << 0) | (p2 << 1) | (d1 << 2) |
                       (p3 << 3) | (d2 << 4) | (d3 << 5) | (d4 << 6);

    return codeword;
}


uint32_t encode_golay(uint16_t data[12]){
int i,j;
uint32_t codeword=0;

for(i=0;i<23;i++){
    codeword<<=1;
    for(j=0;j<12;j++)
codeword^=data[j]&G_golay[j][i];
}

return codeword;
}


/* ---------- メイン復号ルーチン ---------- */

/* matrix: [N_ROWS][N_COLS] 0/1 bit matrix
   returns: number of iterations performed; decoded matrix in-place */
int decode_product_iterative(uint8_t matrix[N_ROWS][N_COLS], int max_iters) {
    int it;
    for (it = 0; it < max_iters; it++) {
        int total_flips = 0;
        printf("\n\n");
        for(int ii=0;ii<23;ii++){
            for(int jj=0;jj<7;jj++)
            printf("%d,",matrix[ii][jj]);
        printf("\n");
        }
        printf("\n\n");

        // --- 行パス: 各行 (長さ N_COLS) を Hamming で復号 ---
        uint8_t row[N_COLS];
        for (int r = 0; r < N_ROWS; r++) {
            for (int c = 0; c < N_COLS; c++) row[c] = matrix[r][c];
            int syn = compute_row_syndrome(row);          // 0..7
            if(syn>0)
                printf("@");
            uint8_t errmask = hamming_lut_lookup(syn);   // 7-bit mask
            if (errmask) {
                printf("rrr=%b\n",errmask);
                // apply mask to row and count flips
                for (int c = 0; c < N_COLS; c++) {
                    if (errmask & (1u << c)) {
                        matrix[r][c] ^= 1;
                        total_flips++;
                    }
                }
            }
        }

        // --- 列パス: 各列 (長さ N_ROWS) を Golay で復号 ---
        uint8_t col[N_ROWS];
        int first=1;
        for (int c = 0; c < N_COLS; c++) {
            for (int r = 0; r < N_ROWS; r++) col[r] = matrix[r][c];
            uint16_t syn = compute_col_syndrome(col);    // 11-bit
            if(syn>0)
                printf("#");
            uint32_t errmask23 = golay_lut_lookup(syn); // lower 23 bits
            if (errmask23) {
                if(first==1)
                printf("eee=%b\n",errmask23);
                first=0;
                for (int r = 0; r < N_ROWS; r++) {
                    if (errmask23 & (1u << r)) {
                        matrix[r][c] ^= 1;
                        total_flips++;
                    }
                }
            }
        }

        if (total_flips == 0) {
            printf("\n\n");
            for(int jj=0;jj<23;jj++){
                for(int ii=0;ii<7;ii++){
                    printf("%d,",matrix[jj][ii]);
                }
            printf("\n");
            }
            // 収束
            return it + 1;
        }
    }
    return max_iters;
}


// ====== テスト用関数 ======

// 符号語を作る: m(長さk) * G(k×n) = c(長さn)
void encode_message(uint8_t *m, int k, uint8_t *G, int n, uint8_t *c) {
    memset(c, 0, n);
    for (int j = 0; j < n; j++) {
        int sum = 0;
        for (int i = 0; i < k; i++) {
            sum ^= (m[i] & G[i*n + j]);
        }
        c[j] = sum;
    }
}

// H * c^T を計算（長さr）
void syndrome(uint8_t *c, int n, uint8_t *H, int r, uint8_t *s) {
    for (int i = 0; i < r; i++) {
        int sum = 0;
        for (int j = 0; j < n; j++) {
            sum ^= (c[j] & H[i*n + j]);
        }
        s[i] = sum;
    }
}

// 1ビット誤りを訂正できるか（デモ用）
void flip_bit(uint8_t *c, int pos) {
    c[pos] ^= 1;
}

// ====== テストランナー ======
void unit_test(uint8_t *G, int k, int n, uint8_t *H, int r) {
    printf("=== Unit Test Start ===\n");

    // 1. メッセージを準備（例: 000...1）
    uint8_t m[256] = {0};
    m[0] = 1; // 簡単のため最初のビットを1に

    // 2. 符号語に変換
    uint8_t c[512];
    encode_message(m, k, G, n, c);

    printf("Encoded codeword (first 32 bits): ");
    for (int i = 0; i < 7*23 && i < n; i++) printf("%d", c[i]);
    printf("...\n");

    // 3. 検査（シンドローム）
    uint8_t s[512];
    syndrome(c, n, H, r, s);

    int ok = 1;
    for (int i = 0; i < r; i++) if (s[i] != 0) { ok = 0; printf("Ah!\n"); break; }
    printf("Check codeword: %s\n", ok ? "PASS (H*c^T=0)" : "FAIL");

    // 4. 1ビット誤りを導入
    flip_bit(c, 0);
    printf("Introduced 1-bit error at pos 0\n");

    syndrome(c, n, H, r, s);
    int nonzero = 0;
    for (int i = 0; i < r; i++) if (s[i] != 0) { nonzero = 1; break; }
    printf("Syndrome after error: %s\n", nonzero ? "nonzero (OK)" : "zero (NG)");

    printf("=== Unit Test End ===\n\n");
}

// example: map a flat codeword (length 161) into matrix[23][7]
void map_codeword_to_matrix(uint8_t *code, uint8_t matrix[N_ROWS][N_COLS]) {
    for (int i = 0; i < N_ROWS; i++)
        for (int j = 0; j < N_COLS; j++)
            matrix[i][j] = code[i * N_COLS + j];
}

/* ---------- Example: small driver (testing) ---------- */
int main(void) {
    // example: allocate matrix and fill with some received bits
    uint8_t matrix[N_ROWS][N_COLS];
    memset(matrix, 0, sizeof(matrix));
    uint8_t C[12*4][23*7]={0};
    uint8_t H[23*7-12*4][23*7]={0};

    // TODO: fill matrix with your received (noisy) bits
    int d=encode_hamming74(0b11);
    d^=1;
    int e=(decode_hamming74(d)>>3);
    printf("m=%b\n",e);
    //pro(G_golay,H_golay);
    kronecker(G_golay,12,23,G_hamming,4,7,&C[0][0]);
    //product_code(G_golay,12,23,G_hamming,4,7,&C);
    uint8_t piv[161];
    //gauss_elim(C,48,161,piv);
    for(int i=0;i<12*4;i++)
    {
        //for(int j=0;j<23*7;j++)
        printf("%d,",C[i][i]);
        printf("\n");
    }
    printf("\n\n");
   //exit(1);
    
    int Hr;
    //exit(1);
    //make_H(C,H,&Hr);
    make_parity_check(C,48,161,H);
    check(C,H,113);
    //exit(1);

    //unit_test((uint8_t*)C, 48, 161, (uint8_t*)H, 113);
    //exit(1);

    
    for(int i=0;i<8;i++)
    C[0][i*7+i]^=1;
    for(int i=0;i<7;i++)
    C[0][i*23+i]^=1;
    for(int j=0;j<23;j++){
        for(int i=0;i<7;i++){
        matrix[j][i]=C[0][j*7+i];
        printf("%d,",matrix[j][i]);
        }
        printf("\n");
    }
    int iters = decode_product_iterative(matrix, 10);
    printf("done iterations: %d\n", iters);

    // show decoded matrix
    for (int c = 0; c < N_COLS; c++){
        for (int r = 0; r < N_ROWS; r++) {
            printf("%d", matrix[r][c]);
        }
        printf("\n");
    }
    return 0;
}
