#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//符号のパラーメータの指定。通常[N,K,T]として、
//Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
//を表す。ここではDは符号長にしている。
#define N 8191  // set small prime ex. p=2053
#define M N // puncture code
#define K (1024) // degree of polynomial
#define E (12)   // bit size of prime
#define DEG 2048 // set (K * E) < N
#define T (K / 2) // weight of error vector


unsigned int mat[N][N] = {0};

unsigned int ma[N*2][N] = {0};
int vb[N][N] = {0};
int gt[N][N] = {0};

