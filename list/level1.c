/* sudan_m1_fallback.c
 *
 * - GF(257) 上で重複度=1 補間 (Q = A(x) + B(x) y) を試す。
 * - 補間成功 -> A,B から f(α) = -A(α)/B(α) を K 点で補間 -> f(x) を検算 -> 出力
 * - 補間失敗 -> フォールバックとして受信語の先頭 K 点で Lagrange 補間し f(x) を出力（必ず 1 個返す）
 *
 * コンパイル: gcc -O2 -o sudan_m1_fallback sudan_m1_fallback.c
 * 実行: ./sudan_m1_fallback
 *
 * 注意:
 *  - フォールバックは必ず 1 個返しますが、正しさは保証しません。
 *  - 本格的に信頼できるリスト復号が欲しい場合は degA/degB 増加か多重度 m>1 の実装が必要。
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define P 257

/* === ユーザ調整パラメータ（必要に応じて変更） === */
#define N 8       /* 符号長（受信点数、テストは小さめに） */
#define K 4       /* 復元多項式次数上限+1 (deg_f <= K-1) */
#define degA 3    /* A(x) の最大次数 */
#define degB 3    /* B(x) の最大次数 */

/* === 有限体基本関数 === */
int modp(int a){ int r = a % P; return (r < 0) ? r + P : r; }
int powmod(int a,int e){
    long long res = 1, base = modp(a);
    while(e){
        if(e & 1) res = (res * base) % P;
        base = (base * base) % P;
        e >>= 1;
    }
    return (int)res;
}
int inv(int a){ return powmod(modp(a), P-2); }

/* evaluate polynomial coeff[0..deg] at x */
int poly_eval_coeff(const int *coeff, int deg, int x){
    long long res = 0;
    long long xp = 1;
    for(int i=0;i<=deg;i++){
        res = (res + xp * (coeff[i] % P + P)) % P;
        xp = (xp * x) % P;
    }
    return (int)modp((int)res);
}

/* === ガウス消去・零空間基底作成（mat rows x cols） ===
   mat_in は rows x cols、要素は mod P の値。
   出力: rank、null_dim、null_basis (max_basis x cols)
*/
int gauss_nullspace(int rows, int cols, int mat_in[][cols], int *rank,
                    int max_basis, int null_basis[][cols], int *null_dim) {
    int mat[rows][cols];
    for(int i=0;i<rows;i++) for(int j=0;j<cols;j++) mat[i][j] = modp(mat_in[i][j]);

    int where[cols]; for(int j=0;j<cols;j++) where[j] = -1;
    int r = 0, c = 0;
    while(r < rows && c < cols){
        int sel = -1;
        for(int i=r;i<rows;i++) if(mat[i][c] != 0){ sel = i; break; }
        if(sel == -1){ c++; continue; }
        if(sel != r){
            for(int j=c;j<cols;j++){ int tmp = mat[r][j]; mat[r][j] = mat[sel][j]; mat[sel][j] = tmp; }
        }
        int invp = inv(mat[r][c]);
        for(int j=c;j<cols;j++) mat[r][j] = modp(mat[r][j] * invp);
        for(int i=0;i<rows;i++){
            if(i==r) continue;
            if(mat[i][c] != 0){
                int factor = mat[i][c];
                for(int j=c;j<cols;j++){
                    mat[i][j] = modp(mat[i][j] - factor * mat[r][j]);
                }
            }
        }
        where[c] = r;
        r++; c++;
    }
    *rank = r;

    /* free variables are columns with where[j] == -1 */
    int free_idx[cols], fcnt = 0;
    for(int j=0;j<cols;j++) if(where[j] == -1) free_idx[fcnt++] = j;
    *null_dim = fcnt;
    if(fcnt == 0) return 0;

    int bcount = 0;
    for(int t=0; t<fcnt && bcount < max_basis; t++){
        int fv = free_idx[t];
        for(int j=0;j<cols;j++) null_basis[bcount][j] = 0;
        null_basis[bcount][fv] = 1;
        for(int j=0;j<cols;j++){
            if(where[j] != -1){
                int row = where[j];
                int coeff = mat[row][fv]; /* entry in row-reduced matrix */
                null_basis[bcount][j] = modp(-coeff);
            }
        }
        bcount++;
    }
    return 0;
}

/* Lagrange interpolation from t points (xs[],ys[]) producing poly out[0..t-1] */
void lagrange_interpolate(int *xs, int *ys, int t, int *out){
    for(int i=0;i<t;i++) out[i] = 0;
    for(int j=0;j<t;j++){
        int numer[t]; for(int i=0;i<t;i++) numer[i]=0;
        numer[0] = 1; int numer_deg = 0;
        int denom = 1;
        for(int m=0;m<t;m++){
            if(m==j) continue;
            int next[t]; for(int i=0;i<t;i++) next[i]=0;
            for(int a=0;a<=numer_deg;a++){
                next[a] = modp(next[a] - numer[a] * xs[m]);
                next[a+1] = modp(next[a+1] + numer[a]);
            }
            numer_deg++;
            for(int u=0; u<=numer_deg; u++) numer[u] = modp(numer[u]);
            denom = modp(denom * (xs[j] - xs[m]));
        }
        int scale = modp( ys[j] * inv(denom) );
        for(int u=0; u<=numer_deg; u++){
            out[u] = modp(out[u] + scale * numer[u]);
        }
    }
}

/* attempt Sudan m=1 interpolation -> produce one candidate f of degree <= Kwant-1
   return 1 on success (fcoef_out filled), 0 on failure */
int attempt_sudan_m1_one(int alpha[], int r[], int n, int Kwant, int *fcoef_out){
    int vars = (degA + 1) + (degB + 1);
    int mat[N][ (degA+1) + (degB+1) ];
    memset(mat, 0, sizeof(mat));
    for(int i=0;i<n;i++){
        int x = alpha[i];
        long long xp = 1;
        for(int j=0;j<=degA;j++){ mat[i][j] = modp((int)xp); xp = (xp * x) % P; }
        xp = 1;
        for(int j=0;j<=degB;j++){ mat[i][(degA+1) + j] = modp((int)((xp * r[i]) % P)); xp = (xp * x) % P; }
    }

    int max_basis = vars;
    int null_basis[max_basis][vars];
    for(int i=0;i<max_basis;i++) for(int j=0;j<vars;j++) null_basis[i][j]=0;
    int null_dim=0, rank=0;
    gauss_nullspace(n, vars, mat, &rank, max_basis, null_basis, &null_dim);
    if(null_dim == 0) return 0; /* 補間失敗: 非自明解なし */

    /* take the first null basis vector as a concrete solution */
    int unknown[vars];
    for(int j=0;j<vars;j++) unknown[j] = null_basis[0][j];
    int allzero = 1;
    for(int j=0;j<vars;j++) if(unknown[j] != 0){ allzero = 0; break; }
    if(allzero){
        /* try other basis vectors */
        int found = 0;
        for(int b=1;b<null_dim && b<max_basis;b++){
            int az = 1;
            for(int j=0;j<vars;j++) if(null_basis[b][j] != 0){ az = 0; break; }
            if(!az){ for(int j=0;j<vars;j++) unknown[j] = null_basis[b][j]; found=1; break; }
        }
        if(!found) return 0;
    }

    int Acoef[degA+1], Bcoef[degB+1];
    for(int j=0;j<=degA;j++) Acoef[j] = unknown[j];
    for(int j=0;j<=degB;j++) Bcoef[j] = unknown[(degA+1)+j];

    /* find highest-degree nonzero in B and scale to make it 1 */
    int lead = -1;
    for(int j=degB;j>=0;j--) if(Bcoef[j] != 0){ lead = j; break; }
    if(lead == -1) return 0; /* B == 0 -> invalid */
    int s = inv(Bcoef[lead]);
    for(int j=0;j<=degA;j++) Acoef[j] = modp(Acoef[j] * s);
    for(int j=0;j<=degB;j++) Bcoef[j] = modp(Bcoef[j] * s);

    /* collect Kwant sample points where B(alpha) != 0 */
    int xs[K], ys[K], got = 0;
    for(int i=0;i<n && got < Kwant; i++){
        int x = alpha[i];
        int Bv = poly_eval_coeff(Bcoef, degB, x);
        if(Bv != 0){
            int Av = poly_eval_coeff(Acoef, degA, x);
            ys[got] = modp( - (long long)Av * inv(Bv) );
            xs[got] = x;
            got++;
        }
    }
    if(got < Kwant) return 0; /* insufficient sample points */

    int fpoly[K];
    for(int i=0;i<K;i++) fpoly[i] = 0;
    lagrange_interpolate(xs, ys, Kwant, fpoly);

    /* verify Q(alpha, f(alpha)) == 0 for all alpha */
    for(int i=0;i<n;i++){
        int x = alpha[i];
        int Av = poly_eval_coeff(Acoef, degA, x);
        int Bv = poly_eval_coeff(Bcoef, degB, x);
        long long xp = 1, sum = 0;
        for(int t=0;t<Kwant;t++){ sum = (sum + xp * fpoly[t]) % P; xp = (xp * x) % P; }
        int fval = modp((int)sum);
        int Qv = modp( Av + (long long)Bv * fval );
        if(Qv != 0) return 0; /* 検算失敗 */
    }

    for(int t=0;t<Kwant;t++) fcoef_out[t] = modp(fpoly[t]);
    return 1;
}

/* fallback: take first K points (indices 0..K-1) and Lagrange interpolate -> always returns one poly */
void fallback_interpolate(int alpha[], int r[], int n, int Kwant, int *fcoef_out){
    int xs[Kwant], ys[Kwant];
    int got = 0;
    for(int i=0;i<n && got < Kwant; i++){
        xs[got] = alpha[i];
        ys[got] = r[i];
        got++;
    }
    lagrange_interpolate(xs, ys, Kwant, fcoef_out);
}

/* print polynomial coeffs */
void print_poly(int *coef, int deg){
    int empty = 1;
    for(int i=0;i<=deg;i++){
        if(coef[i] != 0){ printf("%d*x^%d ", coef[i], i); empty = 0; }
    }
    if(empty) printf("0");
    printf("\n");
}

/* === main: demo with example r[]; modify r[] as needed === */
int main(void){
    int alpha[N], r[N];
    for(int i=0;i<N;i++) alpha[i] = i+1;

    /* Example true message f(x) = 3 + 2x + 5x^2 */
    for(int i=0;i<N;i++){
        int x = alpha[i];
        r[i] = modp(3 + 2*x + 5*x*x);
    }
    /* Inject some errors for demo */
    r[2] = modp(r[2] + 7);
    r[5] = modp(r[5] + 9);

    int fcoef[K];
    int ok = attempt_sudan_m1_one(alpha, r, N, K, fcoef);
    if(ok){
        printf("[Sudan m=1] Found candidate f(x): ");
        print_poly(fcoef, K-1);
        return 0;
    }

    printf("Sudan m=1 補間失敗または検算失敗 → フォールバックで候補生成します。\n");
    fallback_interpolate(alpha, r, N, K, fcoef);
    printf("[Fallback] Candidate f(x) (interpolated from first K points): ");
    print_poly(fcoef, K-1);
    return 0;
}
