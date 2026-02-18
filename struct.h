
/* -*- mode: C; coding:utf-8 -*- */

//monomial
typedef struct
{
   int n; //単項式の次数
   int a; //単項式の係数
} oterm;

//polynomial
typedef struct
{
  oterm t[DEG]; //単項式の配列として多項式を表現する
} OP;

typedef union
{
    int x[DEG]; //配列の添字を次数に、配列の値を係数に持つ多項式の表現
    short a[DEG*2];
    char c[DEG*4];
} vec;

typedef struct {
  vec q;
  vec r;
} rem;


#define I8T char
#define U8C(v) (v##U)

#define U8V(v) (( char)(v)&U8C(0xFF))
#define ROTL8(v, n) \
  (U8V((v) << (n)) | ((v) >> (8 - (n))))

#define R(x, n) (((x) << (n)) | ((x) >> (32 - (n))))


 int rotate_left( int x, int n)
{
  assert(0 < n && n < 32);
  return (x << n) | (x >> (32 - n));
}


typedef struct {
   int x[10][10];
   vec z[N];
  OP f;
  int row; //行
  int col; //列
  int flg;
} MTX;


typedef struct {
   float x[N][N];
  OP f;
  int row; //行
  int col; //列
  int flg;
} MTA;


typedef union {
    unsigned long long int x[K/2];
    unsigned d[K];
    unsigned short a[K*2];
    unsigned char c[K*4];
} uni;

typedef struct
{
    vec f;
    vec g;
    vec h;
} ymo;

