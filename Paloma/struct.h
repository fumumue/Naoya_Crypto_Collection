
/* -*- mode: C; coding:utf-8 -*- */
#include <assert.h>


//monomial
typedef struct
{
  unsigned short n; //単項式の次数
  unsigned short a; //単項式の係数
} oterm;

//polynomial
typedef struct
{
  oterm t[3333]; //単項式の配列として多項式を表現する
} OP;

typedef struct 
{
  unsigned short x[3333]; //配列の添字を次数に、配列の値を係数に持つ多項式の表現
} vec;

typedef struct {
  OP q;
  OP r;
} rem;

typedef struct
{
  unsigned short v[N];
  int f;
} MT;

//extra gcd
typedef struct
{
  OP u; //inverse of polynomial?
  OP v; //error locater
 OP d; //gcd
} EX;

typedef union
{ //test(SIMN)
  unsigned long long int u[256];
  unsigned short s[1024];
} SU;

typedef union
{
  __uint128_t z[K/16];
  unsigned long long int u[K/8];
  unsigned int t[K/4];
  unsigned short x[K/2];
  unsigned char d[K];
} arrayul;

typedef union
{
  unsigned long long int u[N];
  unsigned char d[N * 8];
} arrayull;

typedef union a4
{
  unsigned char ar[4];
  unsigned int n;
} array;

typedef struct a8
{
  unsigned char ar[8];
} array8;

typedef union
{
  unsigned int h[16];
  unsigned char c[64];
} array16;

typedef union aN
{
  unsigned int d[64];
  unsigned long long int u[32];
  unsigned char ar[256];
  //
} arrayn;

typedef struct pub
{
  unsigned char a[K];
  unsigned char b[K];
} set;


#define I8T char
#define U8C(v) (v##U)

#define U8V(v) ((unsigned char)(v)&U8C(0xFF))
#define ROTL8(v, n) \
  (U8V((v) << (n)) | ((v) >> (8 - (n))))

#define R(x, n) (((x) << (n)) | ((x) >> (32 - (n))))

/*
unsigned int rotate_left(unsigned int x, int n)
{
  assert(0 < n && n < 32);
  return (x << n) | (x >> (32 - n));
}
*/

typedef struct {
  unsigned short x[N][N];
  OP f;
  int row; //行
  int col; //列
} MTX;

