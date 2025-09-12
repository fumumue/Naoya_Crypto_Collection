#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//符号のパラーメータの指定。通常[N,K,T]として、
//Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
//を表す。ここではDは符号長にしている。
#define N 257 // set small prime ex. p=2053
#define M N // puncture code
#define K (4) // degree of polynomial
#define E (3)   // bit size of prime
#define DEG N*2 // set (K * E) < N
#define T (K / 2) // weight of error vector
#define Q K*2

unsigned int mat[N*2][N*E] = {0};
unsigned int ma[N*2][N*E] = {0};
int vb[N * 2][N*2] = {0};
int gt[N * 2][N * 2] = {0};

