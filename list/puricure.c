#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define P 257
#define N 8       // 符号長（テスト用に小さめ）
#define K 4       // 情報長（多項式次数+1）

// 有限体 mod P
int modp(int a){ int r=a%P; return r<0?r+P:r; }
int powmod(int a,int e){
    long long res=1, base=modp(a);
    while(e){ if(e&1) res=(res*base)%P; base=(base*base)%P; e>>=1; }
    return (int)res;
}
int inv(int a){ return powmod(modp(a), P-2); }

// ---------- 符号語から補間多項式 Q(x,y) を構築（簡易版） ----------
typedef struct {
    int deg;
    int coef[K]; // 多項式 f(x) の係数
} Poly;

// Lagrange 補間で r=(r1,...,rn) を通る f(x) を復元
void interpolate_poly(int alpha[], int r[], int n, Poly *f){
    memset(f->coef,0,sizeof(f->coef));
    f->deg = K-1;
    for(int j=0;j<K;j++){
        int numer[K]; memset(numer,0,sizeof(numer)); numer[0]=1; int numer_deg=0;
        int denom=1;
        for(int m=0;m<K;m++) if(m!=j){
            // multiply by (x - alpha[m])
            int next[K]; memset(next,0,sizeof(next));
            for(int a=0;a<=numer_deg;a++){
                next[a] = modp(next[a] - numer[a]*alpha[m]);
                next[a+1] = modp(next[a+1] + numer[a]);
            }
            numer_deg++;
            for(int t=0;t<=numer_deg;t++) numer[t]=modp(next[t]);
            denom = modp(denom * (alpha[j]-alpha[m]));
        }
        int scale = modp(r[j]*inv(denom));
        for(int t=0;t<=numer_deg;t++) f->coef[t]=modp(f->coef[t]+scale*numer[t]);
    }
}

// 符号語を入力して補間多項式を計算し候補を表示
void decode_and_print(int alpha[], int r[], int n){
    Poly f;
    interpolate_poly(alpha,r,K,&f);
    printf("候補 f(x) の係数: ");
    for(int i=0;i<=f.deg;i++) printf("%d*x^%d ", f.coef[i], i);
    printf("\n");
}

int main(){
    // 符号語（テスト用、実際には受信語 r を入れる）
    int alpha[N], r[N];
    for(int i=0;i<N;i++){ alpha[i]=i+1; }
    // 本来は符号語=情報多項式を評価したもの。ここでは例として f(x)=3+2x+5x^2
    for(int i=0;i<N;i++){
        int x=alpha[i];
        r[i] = modp(3 + 2*x + 5*x*x);
    }

    decode_and_print(alpha,r,N);
    return 0;
}
