#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define P 257        // GF(257)
#define N 8          // 符号長 (例)
#define K 4          // 情報長
#define M 3          // 最大多重度

// 有限体 mod P
int modp(int a){ int r=a%P; return r<0?r+P:r; }
int powmod(int a,int e){
    long long res=1, base=modp(a);
    while(e){ if(e&1) res=(res*base)%P; base=(base*base)%P; e>>=1; }
    return (int)res;
}
int inv(int a){ return powmod(modp(a), P-2); }

// 多項式構造体
typedef struct {
    int deg;
    int coef[K]; // deg <= K-1 を仮定
} Poly;

// 評価 f(alpha)
int poly_eval(Poly *f, int x){
    long long res=0, powx=1;
    for(int i=0;i<=f->deg;i++){
        res += (long long)f->coef[i]*powx;
        powx = (powx*x)%P;
    }
    return modp((int)res);
}

// 候補生成: 次数 <K のすべての多項式を列挙（単純版）
void brute_force_candidates(int alpha[], int r[], int n){
    Poly f;
    f.deg = K-1;
    printf("候補符号語一覧:\n");
    int count=0;

    // 単純に各係数を小さい範囲で探索（例: 0..P-1 だと大きすぎるので縮小）
    // 実際の因数分解を置き換えるためにサンプルとして探索幅を制限
    for(int a0=0;a0<10;a0++){
      for(int a1=0;a1<10;a1++){
        for(int a2=0;a2<10;a2++){
          for(int a3=0;a3<10;a3++){
            f.coef[0]=a0; f.coef[1]=a1; f.coef[2]=a2; f.coef[3]=a3;

            // 一致数をカウント
            int matches=0;
            for(int i=0;i<n;i++){
              int val=poly_eval(&f,alpha[i]);
              if(val==r[i]) matches++;
            }
            // 多重度Mに基づき閾値を設定（簡略化）
            if(matches>=K-1){
              printf("f(x) = ");
              for(int d=0;d<K;d++) printf("%d*x^%d ", f.coef[d], d);
              printf("\n");
              count++;
            }
          }
        }
      }
    }
    printf("合計 %d 候補\n", count);
}

int main(){
    // 評価点
    int alpha[N], r[N];
    for(int i=0;i<N;i++) alpha[i]=i+1;

    // 元の多項式: f(x)=3+2x+5x^2
    for(int i=0;i<N;i++){
        int x=alpha[i];
        r[i]=modp(3 + 2*x + 5*x*x);
    }
    // 少しエラーを注入
    r[2]=modp(r[2]+7);
    r[5]=modp(r[5]+9);

    brute_force_candidates(alpha,r,N);
    return 0;
}
