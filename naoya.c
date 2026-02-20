#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <math.h>
#include <assert.h>
#include <x86intrin.h> // SIMD命令を使用するためのヘッダファイル

#include "golay.c"
#include "hqc_golay.c"
#include "inv_mat.c"
#include "17.c"


#define SEPARABLE 0
#define MATRIX_SIZE K
#define SHM_KEY 128

#define Prim 3

int g[K + 1] = {0};

// ランダム多項式の生成
static void
ginit(void)
{
    int j, count = 0, k = 0;
    int gg[K + 1] = {0};

    printf("in ginit\n");

    g[K] = 1;          // xor128();
    g[0] = rand() % N; // or N
    k = rand() % (K - 1);
    if (k > 0)
    {
        while (count < k)
        {
            printf("in whule\n");
            j = rand() % (K);
            if (j < K && j > 0 && g[j] == 0)
            {
                g[j] = rand() % N; // or N;
                count++;
            }
        }
    }

    for (j = 0; j < K + 1; j++)
        gg[j] = g[K - j];

    memcpy(g, gg, sizeof(g));
}


// 有限体の元の平方を計算する
int isqrt(int u)
{
  int i, j, k;

  for (i = 0; i < N; i++)
  {
    if ((i*i)%N == u)
      return i;
  }

  printf("来ちゃいけないところに来ました\n");
  exit(1);
}

//内積
int intermul(vec a,vec b){
int i;
int c=0;

for(i=0;i<K;i++)
c+=(a.x[i]*b.x[i])%N;

return c%N;
}

// ベクトルのノルム
int norm(vec a){
int i,n=0;

n+=a.x[i]*a.x[i]%N;
n=isqrt(n);

return n;
}
 int oinv2(int a, int R)
{
     int i;

    if (a == 0)
        return 0;
    //return inv(a,R);


    if (a < 0)
    {
        //printf("a=%d", a);
        a = R + a%R;
        //printf("-a=%d\n", a);
        // exit(1);
    }
     if (a == 1)
         return 1;
    for (i = 1; i < R; i++)
    {
        if ((i * a) % R == 1)
            return i;
    }

    printf("no return2\n");
    exit(1);
}

// 多項式を表示する(default)
void printpol(vec a)
{
    int i, n,flg=0;

    n = deg(a);

    // printf ("baka\n");
    //  assert(("baka\n", n >= 0));

    for (i = n; i > -1; i--)
    {
        if (a.x[i] != 0)
        {
            printf("%d*", a.x[i]);
            // if (i > 0)
            printf("x^%d", i);
            if(i>0)
            printf("+");
        }
    }
    //  printf("\n");

    return;
}

// 多項式を表示する(default)
void printpoln(vec a)
{
    int i, n,flg=0;

    n = deg(a);

    // printf ("baka\n");
    //  assert(("baka\n", n >= 0));

    for (i = n; i > -1; i--)
    {
        if (a.x[i] != 0)
        {
            printf("%d*", a.x[i]);
            // if (i > 0)
            printf("x^%d", i);
            if(i>0)
            printf("+");
        }
    }
      printf("\n");

    return;
}



vec cof( int R, vec f)
{
    int i, k;
    vec b = {0}, h = {0};

    printf("R=%d\n", R);
    // exit(1);
    b = f; // o2v(f);
    k = deg(b);
    printpol(b);
    printf(" =b debugi\n");

    for (i = 0; i < k + 1; i++)
    {
        while(b.x[i]<0)
        b.x[i]+=R;
        h.x[i] = (b.x[i]) % R;
    }
    // g = v2o(h);
    printpol(h);
    printf(" =h in cof\n");
    return h;
}

vec kof( int c, vec f)
{
    int i, k;
    vec b = {0}, h = {0};

    printf("c=%d\n", c);
    // exit(1);
    b = f; // o2v(f);
    k = deg(b);
    printpol(b);
    printf(" =b debugi\n");
    for (i = 0; i < k + 1; i++)
    {
        h.x[i] = (c * b.x[i]) % N;
    }
    // g = v2o(h);
    printpol(h);
    printf(" =h in oinv2\n");
    return h;
}

vec kof2( int c, vec f)
{
    int i, k;
    vec b = {0}, h = {0};

    c = inv(c, N);
    printf("c=%d\n", c);
    // exit(1);
    b = f; // o2v(f);
    k = deg(b);
    printpol(b);
    printf(" =b debugi\n");
    for (i = 0; i < k + 1; i++)
    {
        h.x[i] = (c * b.x[i]) % N;
    }
    // g = v2o(h);
    printpol(h);
    printf(" =h in oinv2\n");
    return h;
}

vec vadd2(vec a, vec b,int R)
{
    int i;
    vec c = {0};

    // printf("deg=%d %d\n",deg(a),deg(b));

    for (i = 0; i < DEG; i++)
        c.x[i] = (a.x[i] + b.x[i]) % R;

    return c;
}

vec lsft(vec a)
{
    vec b = {0};
    int o = deg(a);

    for (int i = 0; i < o + 1; i++)
    {
        b.x[i + 1] = a.x[i];
    }
    // b.x[K*2]=0;

    return b;
}

vec rsft(vec a)
{
    vec b = {0};
    int o = deg(a);

    for (int i = 0; i < o + 1; i++)
        b.x[i] = a.x[i + 1];
    // b.x[0]=0;

    return b;
}


vec convolution( vec a, vec b, int n ) {
    int k, j;
    vec r={0};


    for( k = 0; k < K+1; k++ ){
        for( j = 0; j < K+1; j++ )
            r.x[(k+j) % K] += ( a.x[k] * b.x[j] % n );
    }

    for( k = 0; k < K+1; k++ ){
        r.x[k] = ( r.x[k] % n );
    }

    return r;
}


 // RS-Code generater
 void van(int kk)
 {
     int i, j;

     printf("van der\n");

     for (i = 0; i < kk; i++)
     {
         //mat[i][0] = vb[0][i] = 1; //mltn(i,2);
         //mat[i][1]=vb[1][i]=mltn(i,2);
         //printf("%d,", vb[0][i]);
     }
     printf("\n");

     // #pragma omp parallel for private(i, j)
    for (j = 0; j < K; j++)
     {
        for (i = 0; i < N; i++)
        {
             vb[j][i] = mltn(j, mltn((i)%N,Prim));
            if(i==N-1 && j>0)
            vb[j][i]=0;
             //printf("a%d ",(j+1)%N);
             printf("g%d,", vb[j][i]);
             mat[j][i] = vb[j][i];
         }
         printf("\n");
    }
         for(i=0;i<kk;i++){
         for(j=0;j<N;j++)
         printf("%d,",mat[i][j]);
        printf("\n");
         }
         printf("\n");

         //exit(1);
 }


 // RS-List Decoder
 void list(vec v)
 {
     int i, j;

     printf("van der\n");
     printpoln(v);

     // #pragma omp parallel for private(i, j)
    for (j = 0; j < N; j++)
     {
        for (i = 0; i < N+1; i++)
        {
             vb[i][j] = mltn(i, (j+1)%N);
             //if(i==N-1 && j>0)
             //vb[j][i]=0;
            printf("f%d,", vb[i][j]);
             mat[j][i] = vb[i][j];
         }
         printf("\n");
    }
    for(j=0;j<T+1;j++){
        for(i=0;i<N;i++){
        mat[i][N-T+j]=v.x[i]*mltn(j,mat[i][1])%N;
        }
    }


    printf("mat=\n");
    for(i=0;i<N;i++){
        for(j=0;j<N+1;j++)
        printf("b%d,",mat[i][j]);
       printf("\n");
        }
        printf("\n");

        printf("vb=\n");
        for(i=0;i<N;i++){
            for(j=0;j<N+1;j++)
            printf("b%d,",vb[i][j]);
           printf("\n");
            }
            printf("\n");
//exit(1);
}



 // RS-Code generater
 vec monde(vec x)
 {
    int i, j;
    MTX A={0};
    vec y={0};

     printf("geer\n");

    for(i=0;i<N;i++)
    printf("%d,",x.x[i]);
    printf("\n");
    // exit(1);
    // #pragma omp parallel for private(i, j)
    for (i = 0; i < N-T; i++)
    {
        for (j = 1; j < N+1; j++)
        {
            mat[j-1][i] = mltn(i, j);
            // mat[i-1][j-1] = vb[i-1][j];
         }
         printf("\n");
     }
    for (j = 1; j < N+1; j++)
        {
            for (i = 0; i < K; i++)
            {
                    //A.x[j-1][i+N-T]
             y.x[j-1]+= mltn(i,j)*x.x[i]%N;
             // mat[i-1][j-1] = vb[i-1][j];
             printf("A%d,%d %d %d\n",mltn(i,j),x.x[j-1],i,j);
          }
          y.x[j-1]%=N;
          printf("y=%d\n",y.x[j-1]);
      }

     //for(i=0;i<N;i++)
     {
        for(j=0;j<N+1;j++)
        printf("h%d,", y.x[j]);
        printf("\n");
     }
     //exit(1);

     return y;
 }


void ogt(int kk)
{
    int i, j;

    // #pragma omp parallel for private(i, j)
    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < kk - i; j++)
        {
            gt[i][j + i] = g[j];
        }
    }
    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < kk; j++)
            printf("h%d,", gt[i][j]);
        printf("\n");
    }
    // exit(1);
}


// 配列の値を係数として多項式に設定する
vec setpol( int f[], int n)
{
    OP g;
    vec v = {0};
    int i;

    for (i = 0; i < n; i++)
    {
        v.x[n - 1 - i] = f[i];
    }


    return v;
}

vec mkpol2(int s){
    vec a={0};
    int i;
    for(i=0;i<s;i++)
    a.x[i]=rand()%N;
a.x[s]=1;

return a;
}

vec mkpol3(int s,int R){
    vec a={0};
    int i;
    for(i=0;i<s;i++)
    a.x[i]=rand()%R;
a.x[s]=1;

return a;
}


vec mkpol()
{
    int i, j, k, flg, ii = 0;
    vec w = {0};

    do
    {
        // fail = 0;
        j = 0;
        k = 0;
        flg = 0;
        // l = 0;
        memset(g, 0, sizeof(g));
        // memset(ta, 0, sizeof(ta));
        memset(w.x, 0, sizeof(w));
        ginit();
        ii++;
        if (ii > 100)
        {
            printf("erro=%d\n", ii);
            exit(1);
        }

        for (i = 0; i < K; i++)
        {
            if (g[K - 1] > 0)
                flg = 1;
            if (i % 2 == 1 && g[i] > 0 && i < K)
                k++;
        }

        // 偶数項だけにならないようにする
        if ((k > 0 && flg == 0) || (k > 1 && flg == 1))
        // if(k>0)
        {
            w = setpol(g, K + 1);
            j = 1;
            // if(isquad(w)==-1)
            // exit(1);
        }
        // exit(1);

    } while (j == 0);

    printpol((w));
    printf(" ==g\n");
    // exit(1);

    return w;
}



void printsage(vec a)
{
    int i, j=deg(a),flg=0;
    oterm b;

    printf("poly=");
    for (i = 0; i < DEG; i++)
    {
        if (a.x[i] != 0)
        {
            printf("%d*", a.x[i]);
            // if (i > 0)
            printf("x^%d", i);
            if (i < j)
            printf("+");
            flg=1;
            //printf("+%d*X**%d", a.x[i], i); // for GF(2^m)
        }
    }
}

// 多項式の代入値
long long int
trace(vec f, long long int x)
{
    long long int u = 0;
    vec v = (f);
    int d = deg((v)) + 1;

    for (int i = 0; i < d; i++)
    {
        if (v.x[i] > 0){
            u = (u + (v.x[i] * mltn(i, x))) % N;
            //printf("u%d %d\n",u,i);
        }
    }

    return u;
}


// 多項式の代入値
unsigned int
trace2(vec f, unsigned int x)
{
    unsigned int u = 0;
    vec v = (f);
    int d = deg((v)) + 1;

    for (int i = 0; i < d; i++)
    {
        if (v.x[i] > 0){
            u = (u + (v.x[i] * mltn(i, x)));
            //printf("u%d %d\n",u,i);
        }
    }

    return u;
}


// aに何をかけたらbになるか
 int
equ2( int a,  int b,int R)
{
    int i;
    // for(int i=0;i<N;i++)
    if (b == 0)
        return 0;
    if (a == 1)
        return b;
    if(a==b)
        return 1;

    return (inv(a, R) * b) % R;

}

// 多項式を単行式で割る
oterm vLTdiv2(vec f, oterm t,int R)
{
    oterm tt = {0}, s = {
                        0};

    tt = vLT(f);
    if (tt.n < t.n)
    {
        s.n = 0;
        s.a = 0;
    }
    else if (tt.n == t.n)
    {
        s.n = 0;
        s.a = equ2(t.a, tt.a,R);
    }
    else if (tt.n > t.n)
    {
        s.n = tt.n - t.n;
        s.a = equ2(t.a, tt.a,R);
        printf("ss-%d %d %d\n",s.a,t.a,tt.a);
    }
    else if (t.n == 0 && t.a > 0)
    {
        s.a = (tt.a * inv(t.a, R)) % R;
        s.n = tt.n;
    }

    return s;
}

// 多項式を項ずつ掛ける
vec vterml2(vec f, oterm t,int R)
{
    // f = conv(f);
    // ssert(op_verify(f));
    int i;
    vec h = {0};

    // f=conv(f);
    // k = deg (o2v(f));

    for (i = 0; i < DEG; i++)
    {
        // h.t[i].n = f.t[i].n + t.n;
        if (f.x[i] > 0)
            h.x[i + t.n] = (f.x[i] * t.a) % R;
    }

    // h = conv(h);
    //  assert(op_verify(h));
    return h;
}


// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
vec vsub2(vec a, vec b,int R)
{
    vec c = {0};
    // int i, j, k, l = 0;
    vec h = {0}, f2 = {0}, g2 = {0};

    for (int i = 0; i < DEG; i++)
    {
        if (a.x[i] >= b.x[i])
            c.x[i] = (a.x[i] - b.x[i]) % R;
        if (a.x[i] < b.x[i])
            c.x[i] = (R + a.x[i] - b.x[i]) % R;
    }

    return c;
}

int vm = 0;
// 多項式の剰余を取る
vec vmod(vec f, vec g)
{
    vec h = {0};
    oterm b = {0}, c = {0};

    if (deg(g) == 0)
        return g;
    vm++;
    // printf("vmod-bl=%d k=%d\n",deg(f),deg(g));
    if (vLT(f).n < vLT(g).n)
    {
        //    exit(1);
        return f;
    }

    b = vLT(g);

    // printpol(f);
    // printf(" ==f\n");
    while (1)
    {
        //b = vLT(g);
        // printf("@\n");
        c = vLTdiv(f, b);
        h = vterml(g, c);
        f = vsub(f, h);
        // printsage(g);
        if (deg((f)) == 0 || deg((h)) == 0)
        {
            break;
        }

        if (c.n == 0)
            break;
    }
    // printf("vmod-baka== %d %d\n",deg(f),deg(g));
    return f;
}

// 多項式の剰余を取る
vec vmod2(vec f, vec g,int R)
{
    vec h = {0};
    oterm b = {0}, c = {0};

    if (deg(g) == 0)
        return g;
    vm++;
    // printf("vmod-bl=%d k=%d\n",deg(f),deg(g));
    if (vLT(f).n < vLT(g).n)
    {
        //    exit(1);
        return f;
    }
    if(vLT(f).n==vLT(g).n && vLT(f).a==vLT(g).a)
    return vsub2(f,g,R);

    //b = vLT(g);

    // printpol(f);
    // printf(" ==f\n");
    while (1)
    {
        //printpol(g);
        //printf(" ==gI\n");
        b = vLT(g);
        printpol(f);
        printf(" ==IfI\n");
        // printf("@\n");
        c = vLTdiv2(f, b,R);
        h = vterml2(g, c,R);
        f = vsub2(f, h,R);
        printf("Ic.a=%d,%d b.a=%d %d\n",c.a,c.n,b.a,b.n);
        printpol(f);
        printf(" ==If2\n");
        printpol(g);
        printf(" ==Ig2\n");
        printpol(h);
        printf(" ==Ih\n");
        //if(vLT(h).a==0)
        //exit(1);
        // printsage(g);
        printf("%d %d lol\n",c.a,c.n);
        if (c.n == 0)
            break;
        if (deg((f)) == 0 || deg((h)) == 0)
        {
            break;
        }


    }
    // printf("vmod-baka== %d %d\n",deg(f),deg(g));
    return f;
}


// 多項式のべき乗
vec opow(vec f, int n)
{
    // int i;
    vec g = {0};

    g = f;

    for (int i = 1; i < n; i++)
        g = vmul(g, f,N);

    return g;
}

vec vpowmod(vec f, vec mod, int n)
{
    vec v = {0};
    vec ret = {0};

    v.x[0] = 1;
    ret = (v);
    while (n > 0)
    {
        if (n % 2 == 1)
            ret = vmod(vmul(ret, f,N), mod); // n の最下位bitが 1 ならば x^(2^i) をかける
        f = vmod(vmul(f, f,N), mod);
        n >>= 1; // n を1bit 左にずらす
    }
    return ret;
}

// gcd
vec ogcd(vec xx, vec yy)
{
    vec tt = {0}, tmp, h = {0};
    // ee.x[K] = 1;

    h.x[0] = 1;
    // h.x[0] = 0;
    if (deg((xx)) < deg((yy)))
    {
        tmp = xx;
        xx = yy;
        yy = tmp;
    }
    // tt = vmod(xx, yy);
    tt = vmod(xx, yy);
    while (deg(tt) > 0)
    {
        // printf("Oh!\n");
        xx = yy;
        yy = tt;
        if (deg(yy) > 0)
        {
            tt = vmod(xx, yy);
        }
        if (vLT(tt).a == 0)
            return yy;
    }
    if (vLT(yy).a == 0)
    {
        return tt;
    }
    else
    {
        return h;
    }
    //  return yy;
}

int diag(MTX a, int n)
{
    return (a.x[n][n] * a.x[n + 1][n + 1] - a.x[n][n + 1] * a.x[n + 1][n]) % N;
}

// resultant（シルベスター行列）
int resl(vec f, vec g)
{
    MTX a = {0};
    int dia[N] = {0};
    /*
    f.x[0]=16;
    f.x[1]=0;
    f.x[2]=4;
    f.x[3]=4;
    f.x[4]=1;
    g.x[0]=8;
    g.x[1]=9;
    g.x[2]=10;
    g.x[3]=9;
    printf("\n");
    */
    int n = deg(f), m = deg(g);
    if (n < m)
    {
        for (int i = 0; i < n + 1; i++)
        {
            for (int j = 0; j < m + 1; j++)
            {
                a.x[i + j][i] = f.x[n - j];
            }
        }
        for (int i = 0; i < n + m; i++)
        {
            for (int j = 0; j < n + m; j++)
            {
                a.x[j + i][i + m] = g.x[m - j];
            }
        }
    }
    if (n >= m)
    {
        for (int i = 0; i < m + 1; i++)
        {
            for (int j = 0; j < n + 1; j++)
            {
                a.x[i + j][i] = f.x[n - j];
            }
        }
        for (int i = 0; i < n + m + 1; i++)
        {
            for (int j = 0; j < n + m + 1; j++)
            {
                a.x[j + i][i + m] = g.x[m - j];
            }
        }
    }
    /*
    for(int i=0;i<n+m;i++){
        for(int j=0;j<m+n;j++)
        printf("%d,",a.x[i][j]);
        printf("\n");
    }
    printf("\n");
    */
    int tmp[N] = {0};
    int i, j, k, t;
    for (i = 0; i < m + n - 1; i++)
    {

        for (k = i; k < m + n - 1; k++)
        { // m+n
            // printf("%d ",k);
            t = a.x[k + 1][i];
            for (int j = i; j < n + m; j++)
            {
                tmp[j] = a.x[k + 1][j] - (a.x[i][j] * equ(a.x[i][i], a.x[k + 1][i])) % N;
                // printf("i=%d (j=%d k+1=%d) n=%d ks=%d %d %d t=%d =%d\n",i,j,k+1,a.x[k+1][j],(a.x[i][j]*equ(a.x[i][i],a.x[k+1][i]))%N,a.x[k][j],(a.x[i][j]),t,(N+tmp[j])%N);
            }
            // printf("\n");
            for (int j = 0; j < n + m; j++)
            {
                a.x[k + 1][j] = tmp[j];
                if (a.x[k + 1][j] < 0)
                    a.x[k + 1][j] = N + a.x[k + 1][j];
            }
            /*
            for(int u=0;u<n+m;u++){
                for(int v=0;v<n+m;v++)
                printf("%d ",a.x[u][v]);
                printf("\n");
            }
            printf(" %d %d %d\n",k,m+n,i);
            */
        }
        dia[i] = a.x[i][i];
    }

    int y = diag(a, n + m - 2);

    for (i = 0; i < m + n - 2; i++)
    {
        y = (y * dia[i]) % N;
        if (dia[i] == 0)
            return 0;
    }
    printf("y=%d\n", y);
    // exit(1);
    /*
    vec c=ogcd(f,g);
    if((deg(c)>0 && y>0)){ //} || (deg(c)==0 && y==0)){
    printsage(c);
    printf(" ==baka\n");
    printsage(f);
    printf(" ==f\n");
    printsage(g);
    printf(" ==g\n");
    exit(1);
    }
    */
    if (y > 0)
        return 0;
    if (y == 0)
        return -1;

    return 0;
}

int cnty = 0;
vec vpp(vec f, vec mod, int n)
{
    int i;
    vec s = {0};
    // t = f;
    s = f;
    printf("@\n");
    // 繰り返し２乗法
    for (i = 1; i < n; i++)
    {
        s = vmod(vmul(s, f,N), mod);
    }

    return s;
}

// GCD for decode
vec vgcd(vec xx, vec yy)
{
    vec tt;

    while (deg(yy) > 0)
    {
        tt = vmod(xx, yy);
        xx = yy;
        yy = tt;
    }
    if (yy.x[0] > 0)
        tt = kof2(yy.x[0], xx);
    printpol((yy));
    printf(" =========yy\n");
    printpol((tt));
    printf(" =========tt\n");

    return tt;
}

unsigned int oinv(int a, unsigned int n)
{
    unsigned int i;

    if (a == 0)
        return 0;
    if (a < 0)
    {
        //printf("a=%d", a);
        a = N + a%N;
        //printf("-a=%d\n", a);
        // exit(1);
    }
    // if (a == 1)
    //     return 1;
    for (i = 1; i < n; i++)
    {
        if ((i * a) % N == 1)
            return i;
    }
    printf("no return\n");
    exit(1);
}


// 行列の逆行列を計算する関数
MTX inverseMatrix(MTX A, MTX A_inv, int start_row, int end_row)
{
    int i, j, k;
    int temp;

    // 単位行列を初期化
    for (i = 0; i < K / 2; i++)
    {
        for (j = 0; j < K / 2 + 1; j++)
        {
            A_inv.x[i][j] = (i == j) ? 1 : 0;
        }
    }

    // ガウス・ジョルダン法による逆行列の計算
    for (k = start_row; k < end_row; k++)
    {
        temp = A.x[k][k];
        for (j = 0; j < K / 2 + 1; j++)
        {
            A.x[k][j] = A.x[k][j] * oinv(temp, N) % N;
            A_inv.x[k][j] = A_inv.x[k][j] * oinv(temp, N) % N;
        }
        for (i = start_row; i < end_row; i++)
        {
            if (i != k)
            {
                temp = A.x[i][k] % N;
                for (j = 0; j < K / 2 + 1; j++)
                {
                    A.x[i][j] -= (A.x[k][j] * temp) % N;
                    A_inv.x[i][j] -= (A_inv.x[k][j] * temp) % N;
                }
            }
        }
    }
    vec x = {0};
    for (i = 0; i < K / 2; i++)
    {
        if (N > A.x[i][K / 2])
        {
            x.x[K / 2 - i] = (N - A.x[i][K / 2]) % N;
        }
        else
        {
            x.x[K / 2 - i] = A.x[i][K / 2] % N;
        }
    }
    for (int i = 0; i < K / 2; i++)
    {
        printf("in inverse ");
        for (int j = 0; j < K / 2; j++)
        {
            if (A_inv.x[i][j] < 0)
                A_inv.x[i][j] = N + A_inv.x[i][j]%N;
            printf("%d ", A_inv.x[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    // exit(1);

    /*
        x.x[0] = 1;

        vec vv = {0};
        OP pol = {0};
        pol = setpol(x.x, K / 2 + 1);
        printpol(o2v(pol));
        printf(" ==key\n");
        int key=0;
        for (i = 0; i < N; i++)
        {
            // v.x[i] = 0;
            if (trace(pol, i) % N == 0)
            {
                printf("error position=%d\n", i);
                vv.x[key++] = i;
            }
        }
        for (i = 0; i < K / 2; i++)
        {
            for (j = 0; j < K / 2 + 1; j++)
                printf("%d,", A_inv.x[i][j]%N);
            printf("\n");
        }
        // exit(1);
      */
    return A_inv;
}


// #define NN 16
vec sol(MTX a, int start, int end)
{
    int p, d;
    int i, j, k;
    vec v = {0};

    for (i = start; i < end; i++)
    {
        p = a.x[i][i];

        for (j = 0; j < (K / 2 + 1); j++)
        {
            a.x[i][j] = (a.x[i][j] * inv(p, N)) % N;
        }

        for (j = 0; j < K / 2; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K / 2 + 1); k++)
                {
                    if (a.x[j][k] > (d * a.x[i][k]) % N)
                    {
                        a.x[j][k] -= (d * a.x[i][k]) % N;
                    }
                    else
                    {
                        a.x[j][k] = (N + (a.x[j][k] - (d * a.x[i][k]) % N)) % N;
                    }
                }
            }
        }
    }
    vec x = {0};
    for (i = start; i < end; i++)
    {
        if (N > a.x[i][K / 2])
        {
            x.x[K / 2 - i] = (N - a.x[i][K / 2]) % N;
        }
        else
        {
            x.x[K / 2 - i] = a.x[i][K / 2] % N;
        }
    }

    x.x[0] = 1;

    vec vv = {0};
    vec pol = {0};
    pol = setpol(x.x, K / 2 + 1);
    printpol((pol));
    printf(" ==key\n");
    int key = 0;
    for (i = 0; i < N; i++)
    {
        // v.x[i] = 0;
        if (trace(pol, i) % N == 0)
        {
            printf("error position=%d\n", i);
            vv.x[key++] = i;
        }
    }

    return vv;
}

// #define NN 16
vec sol3(MTX a, int start, int end)
{
    int p, d;
    int i, j, k;
    vec v = {0};

    for (i = start; i < end; i++)
    {
        p = a.x[i][i];

        for (j = 0; j < (K / 2 + 1); j++)
        {
            a.x[i][j] = (a.x[i][j] * inv(p, N)) % N;
        }

        for (j = 0; j < K / 2; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K / 2 + 1); k++)
                {
                    if (a.x[j][k] > (d * a.x[i][k]) % N)
                    {
                        a.x[j][k] -= (d * a.x[i][k]) % N;
                    }
                    else
                    {
                        a.x[j][k] = (N + (a.x[j][k] - (d * a.x[i][k]) % N)) % N;
                    }
                }
            }
        }
    }
    vec x = {0};
    for (i = start; i < end; i++)
    {
        if (N > a.x[i][K / 2])
        {
            x.x[K / 2 - i] = (N - a.x[i][K / 2]) % N;
        }
        else
        {
            x.x[K / 2 - i] = a.x[i][K / 2] % N;
        }
    }

    x.x[0] = 1;

    vec vv = {0};
    vec pol = {0};
    pol = setpol(x.x, K / 2 + 1);
    printpol((pol));
    printf(" ==key\n");
    int key = 0;
    for (i = 0; i < N; i++)
    {
        // v.x[i] = 0;
        if (trace(pol, i) % N == 0)
        {
            printf("error position=%d\n", i);
            vv.x[key++] = i;
        }
    }

    return vv;
}


// #define NN 16
vec sol2(MTX a, int start, int end)
{
    int p, d;
    int i, j, k;
    vec v = {0};

    for (i = start; i < end; i++)
    {
        p = a.x[i][i];

        for (j = 0; j < (K / 2 + 1); j++)
        {
            a.x[i][j] = (a.x[i][j] * inv(p, N)) % N;
        }

        for (j = 0; j < K / 2; j++)
        {
            if (i != j)
            {
                d = a.x[j][i];

                for (k = i; k < (K / 2 + 1); k++)
                {
                    if (a.x[j][k] > (d * a.x[i][k]) % N)
                    {
                        a.x[j][k] -= (d * a.x[i][k]) % N;
                    }
                    else
                    {
                        a.x[j][k] = (N + (a.x[j][k] - (d * a.x[i][k]) % N)) % N;
                    }
                }
            }
        }
    }
    vec x = {0};
    for (i = start; i < end; i++)
    {
        if (N > a.x[i][K / 2])
        {
            x.x[K / 2 - i] = (N - a.x[i][K / 2]) % N;
        }
        else
        {
            x.x[K / 2 - i] = a.x[i][K / 2] % N;
        }
    }

    x.x[0] = 1;

    vec vv = {0};
    vec pol = {0};
    pol = (setpol(x.x, K / 2 + 1));
    printpol((pol));
    printf(" ==key\n");
    int key = 0;
    for (i = 0; i < N; i++)
    {
        // v.x[i] = 0;
        if (trace(pol, i) % N == 0)
        {
            printf("error position=%d\n", i);
            vv.x[key++] = i;
        }
    }

    return vv;
}

// 多項式のべき乗余
vec opowmod(vec f, vec mod, int n)
{
    // int i, j = 0;
    vec g = f;
    printsage(mod);
    printf(" ma\n");
    // 繰り返し２乗法
    for (int i = 1; i < n; i++)
    {
        // f = vmul(f, f);
        g = vmul(g, f,N);
        if (deg(g) > deg(mod))
        {
            // printsage(g);
            // printf(" tadaima!\n");
            g = vmod(g, mod);
            // printsage(g);
            // printf(" tadaima2!\n");
        }
    }
    printsage(g);
    printf(" ==ge!\n");
    // exit(1);
    return g;
}

int is_equ(vec a, vec b)
{
    for (int i = 0; i < N * N; i++)
        if (a.x[i] != b.x[i])
            return -1;

    return 0;
}

// GF(2^m) then set m in this function.
int ben_or(vec f)
{
    int n; //, pid;

    vec s = {0}, u = {0}, r = {0};
    vec v = {0}; //, ff=o2v(f);
    // if GF(8192) is 2^m and m==13 or if GF(4096) and m==12 if GF(16384) is testing
    // int m = E;
    //  m=12 as a for GF(4096)=2^12 defined @ gloal.h or here,for example m=4 and GF(16)

    v.x[1] = 1;
    s = (v);
    // for (int i = 0; i < K / 2; i++)
    r = s;
    n = deg((f));

    if (vLT(f).n == 0)
    {
        printf("f==0\n");
        exit(1);
    }
    if (n == 0)
        return -1;

    // r(x)^{q^i} square pow mod
    for (int i = 0; i < K / 2; i++)
    {
        printf(":i=%d", i);
        // irreducible over GH(8192) 2^13
        // if(r.x[0]==65535)
        // return -1;
        // printsage(r);
        // printf(" --p\n");

        memset(r.x, 0, sizeof(r.x));
        v = vpowmod(v, f, N);
        r = v;
        // r.x[l]=1;

        u = vsub(r, (s));
        u = vmod(u, f);

        if (deg(u) > 0)
        {
            // printsage(u);
            // printf(" you\n");
            // printsage(f);
            printf(" me\n");
            u = ogcd(f, u);
            // int le=resl(f,u);
            // if(le==0 && deg(u)==0){
            //     printf("baka^^\n");
            // exit(1);
            // return -1;
            // }
            printf("you\n");
        }
        else
        {
            return -1;
        }
        if (deg(u) > 0) //  || vLT(u).a > 0)
        {
            // if(fequ(u,f)==1)
            {
                // flg[i]= -1;
                printf("ae\n");
                return -1;
            }
        }
    }

    return 0;
}

vec mkd(vec w, int kk)
{
    int i, j, k, l, ii = 0;

    unsigned int tr[N] = {0};
    unsigned int ta[N] = {0};
    vec v = {0}, pp = {0}, tt = {0};
    unsigned int po[K + 1] = {1, 0, 1, 0, 5};
    // vec w={0};
    vec r = {0};

aa:

    // printf("\n");
    memset(mat, 0, sizeof(mat));
    // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
    // 既約多項式しか使わない。

    l = 0;
    ii = 0;
    // irreducible gvecpa code (既役多項式が必要なら、ここのコメントを外すこと。)

    w = mkpol(K);
    //l = ben_or((w));
    l=0;
    while (l == -1)
        goto aa;
    printsage((w));
    printf("\n");
    //exit(1);
    //     printf("wwwwwww\n");
    //  exit(1);
    //  separable gvecpa code
    //  w = mkpol();
    r = (w);
    //  r=vmul(w,w);
    memset(ta, 0, sizeof(ta));
    // w = setpol(g, K + 1);
    printpol((r));
    printf(" =poly\n");
    // exit(1);

    // 多項式の値が0でないことを確認
    for (int i = 0; i < K; i++)
    {
        ta[i] = trace(w, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
            // fail = 1;
            // exit(1);
            goto aa;
        }
    }
    for (int i = 0; i < K; i++)
    {
        tr[i] = inv(ta[i], N);
        // printf("%d,", tr[i]);
    }
    memset(g, 0, sizeof(g));
    // g[0] = 1;

    // 多項式を固定したい場合コメントアウトする。
    printpol(r);
    printf("\n");
    printsage((r));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");
    memset(v.x, 0, sizeof(v.x));
    //  v=rev(w);
    van(kk);
    //  v=(w);
    ogt(kk);
    // exit(1);
    //  wait();

    // #pragma omp parallel for

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (int j = 0; j < K; j++)
    {
        for (int i = 0; i < M; i++)
        {
            ma[i][j] = (vb[j][i] * tr[i]) % N;
        }
    }

    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < M; j++)
        {
            for (int k = 0; k < K; k++)
            {
                mat[j][i] = (mat[j][i] + (gt[k][i] * ma[j][k])) % N;
            }
            printf("c%d,", mat[j][i]);
        }
        printf("\n");
    }

    /*
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < kk; j++)
            {
                mat[j][i] = vb[j][i];
            }
        }
    */
    // printf("\n");
    // exit(1);
    /*
    for( int j = 0; j < N; j++)
    {
        for( int i= 0; i < kk; i++)
            printf("%d,", mat[j][i]);
        printf("\n");
    }
    //exit(1);
    //wait();
*/

    return (w);
}

vec mkd2(vec w, int kk, int start, int end)
{
    int i, j, k, l, ii = 0;

     int tr[N] = {0};
     int ta[N] = {0};
    vec v = {0}, pp = {0}, tt = {0};
     int po[K + 1] = {1, 0, 1, 0, 5};
    // vec w={0};
    vec r = {0};

aa:

    // printf("\n");
    memset(mat, 0, sizeof(mat));
    // 既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
    // 既約多項式しか使わない。

    l = 0;
    ii = 0;
    // irreducible gvecpa code (既役多項式が必要なら、ここのコメントを外すこと。)

    w = mkpol();

    //l = ben_or((w));
    while (l == -1)
        goto aa;
    printsage((w));
    printf("\n");
    //exit(1);

    //     printf("wwwwwww\n");
    //  exit(1);
    //  separable gvecpa code
    //  w = mkpol();
    r = (w);
    //  r=vmul(w,w);
    memset(ta, 0, sizeof(ta));
    // w = setpol(g, K + 1);
    printpol((r));
    printf(" =poly\n");
    // exit(1);

    // 多項式の値が0でないことを確認
    for (int i = start; i < end; i++)
    {
        ta[i] = trace(w, i);
        if (ta[i] == 0)
        {
            printf("eval 0 @ %d\n", i);
            // fail = 1;
            // exit(1);
            goto aa;
        }
    }
    for (int i = start; i < end; i++)
    {
        tr[i] = inv(ta[i], N);
        // printf("%d,", tr[i]);
    }
    memset(g, 0, sizeof(g));
    // g[0] = 1;

    // 多項式を固定したい場合コメントアウトする。
    printpol(r);
    printf("\n");
    printsage((r));
    printf("\n");
    printf("sagemath で既約性を検査してください！\n");
    memset(v.x, 0, sizeof(v.x));
    //  v=rev(w);
    van(kk);
    //  v=(w);
    ogt(kk);
    // exit(1);
    //  wait();

    // #pragma omp parallel for

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (int j = start; j < end; j++)
    {
        for (int i = 0; i < M; i++)
        {
            ma[i][j] = (vb[j][i] * tr[i]) % N;
        }
    }

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < M; j++)
        {
            for (int k = 0; k < K; k++)
            {
                mat[j][i] = (mat[j][i] + (gt[k][i] * ma[j][k])) % N;
            }
            //printf("c%d,", mat[j][i]);
        }
        //printf("\n");
    }

    /*
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < kk; j++)
            {
                mat[j][i] = vb[j][i];
            }
        }
    */
    // printf("\n");
    // exit(1);
    /*
    for( int j = 0; j < N; j++)
    {
        for( int i= 0; i < kk; i++)
            printf("%d,", mat[j][i]);
        printf("\n");
    }
    //exit(1);
    //wait();
*/

    return (w);
}

// Goppa Code's Parity Check (Berlekamp type)
void vv(int kk)
{
    int i, j;
    vec r = mkpol();
     int tr[N];
     int ta[N] = {0};

    printf("van der\n");

    for (i = 0; i < kk; i++)
    {
        for (j = 0; j < N; j++)
        {
            vb[i][j] = mltn(i, j);
        }
        // printf("\n");
    }

    int l = -1;
    vec pp = {0}, tt = {0};

aa:
    // exit(1);
    r = mkpol();

    for (i = 0; i < N; i++)
    {
        ta[i] = trace(r, i);
        if (ta[i] == 0)
        {
            printf("trace 0 @ %d\n", i);
            // fail = 1;
            goto aa;
        }
    }

    for (i = 0; i < N; i++)
    {
        tr[i] = inv(ta[i], N);
        // printf("%d,", tr[i]);
    }

    printf("\nすげ、オレもうイキそ・・・\n");
    // keygen(g);
    // exit(1);

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < kk; j++)
        {
            mat[i][j] = (vb[j][i] * tr[i]) % N;
        }
    }
    for (i = 0; i < K; i++)
    {
        for (j = 0; j < N; j++)
            printf("c%d,", mat[j][i]);
        printf("\n");
    }
}

void mkerr(unsigned int *z1, int num)
{
    int j, l;

    j = 0;

    memset(z1, 0, sizeof(2 * M));

    while (j < num)
    {
        l = rand() % (N - 2)+1;
        // printf ("l=%d\n", l);
        if (0 == z1[l] && l > 0)
        {
            z1[l] = 1; //rand()%N;
            // printf("l=%d\n", l);
            if(z1[l]>0)
            j++;
        }
    }
}

vec synd( int zz[], int kk)
{
     int syn[K] = {0}, s = 0;
    int i, j;
    vec f = {0};

    printf("in synd2\n");

    for (i = 0; i < kk; i++)
    {
        syn[i] = 0;
        s = 0;
        // #pragma omp parallel num_threads(16)
        for (j = 0; j < N; j++)
        {
            s = (s + (zz[j] * mat[j][i])) % N;
        }
        syn[i] = s;
        // printf ("syn%d,", syn[i]);
    }
    // printf ("\n");

    f = setpol(syn, kk);
    printpol((f));
    printf(" syn============= %d\n", deg((f)));
    //  exit(1);

    return f;
}


int ink(int vx){

    if(vx==0)
    return 1;
    for(int i=0;i<N-1;i++){
    if((vx)==mltn(i,Prim)){
        return i;
        }
    }
    //exit(1);
    //return N-1;
}


// chen探索
vec chen(vec f)
{
    vec e = {0};
    int i, n, x = 0, count = 0;
     int z;

    n = deg((f));
    for (x = 0; x < M; x++)
    {
        z = 0;
        for (i = 0; i < n + 1; i++)
        {
            if (f.x[i] > 0)
                z += (mltn(i, x) * f.x[i]) % N;
        }
        if (z % N == 0)
        {
            e.x[count] = ink(x);
            count++;
            printf("change %d\n", ink(x));
        }
    }
    if(count<T)
    {
        printf("kaba\n");
        exit(1);
    }
    return e;
}


vec bms( int s[])
{
    int L = 0, m = -1, d[K] = {0}, k = 0, i, e;
    vec f = {0}, g = {0}, h, v;

    f.x[0] = g.x[0] = 1;

    while (k <= (2 * T - 1))
    {
        e = 0;
        for (i = 0; i < L; i++)
            e = (e+f.x[i]*s[k - i])%N;

        d[k] = (f.x[i]*s[k - i] + e)%N; // s[k] ^ e;
        if (d[k] > 0)
        {
            h = f;
            memset(v.x, 0, sizeof(v.x));
            v.x[k - m] = 1;

             int a;
            a = (m < 0) ? 1 : inv(d[m],N);
            f = vadd2(f, vmul(kof2((d[k]* a), g), v,N),N);
            if (L <= k / 2)
            {
                L = k + 1 - L;
                m = k;
                g = h;
            }
        }
        k++;
    }

    return f;
}


vec pmul(vec a, vec b)
{
    int i, j, k, l;
    vec c = {0};

    k = deg(a) + 1;
    l = deg(b) + 1;
    printf("k=%d,l=%d", k, l);
    // exit(1);
    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < l; j++)
            if (a.x[i] > 0)
            {
                c.x[i + j] = (c.x[i + j] + a.x[i] * b.x[j]) % N;
                // printf("%d=c ",c.x[i+j]);
            }
        // printf("\n");
    }
    /*
    printf("\n");
    printpol(v2o(c));
    printf(" ==c\n");
    printpol(v2o(a));
    printf(" ==a\n");
    printpol(v2o(b));
    printf(" ==b\n");
    // exit(1);
    */
    return c;
}


// simple_bm.c -- 安全な BM 実装（vec ベース）
// 前提: vec 型と vmul, vadd, vsub, deg, coeff access (v.x[i]) がある。
//       体の法は MOD (例: q), inv(a, MOD) が利用可能。
// 入力:  s[0..Slen-1]  (シンドローム, Slen は与える長さ)
// 出力:  f (誤り位置多項式) を vec として返す (次数は返り値で伝えても可)

int bm_vec(const unsigned int *s, int Slen, vec *out_f) {
    // C(x) と B(x) を vec で持つ。C は現在の Λ, B は直前のコピー
    vec C = {0}, B = {0}, Tt = {0};
    // 初期: C(x) = 1, B(x) = 1
    C.x[0] = 1;
    B.x[0] = 1;

    int L = 0;          // 現在の Λ の次数
    int m = 1;          // シフト量
    unsigned int b = 1; // 最後に使った不一致量（discrepancy）

    for (int n = 0; n < Slen; n++) {
        // discrepancy d = S[n] + sum_{i=1..L} C[i]*S[n-i]
        int d = s[n] % N; // N は体の法（例: q）
        for (int i = 1; i <= L; i++) {
            if (C.x[i] == 0) continue;
            int idx = n - i;
            if (idx < 0) break; // S[n-i] が無ければ 0 と仮定（あるいはSlenを確保しておく）
            d = (d + (C.x[i] * s[idx])) % N;
        }
        d = (d % N + N) % N;

        if (d == 0) {
            m++;
        } else {
            // T = C (backup)
            memcpy(&Tt, &C, sizeof(vec));
            // C = C - (d/b) * x^m * B
            int coef = (d * inv(b, N)) % N; // invはmod逆元
            // shift B by m: compute tmp = x^m * B
            vec tmp = {0};
            for (int i = 0; i <= deg(B); i++) {
                if (B.x[i] == 0) continue;
                tmp.x[i + m] = (B.x[i]) % N;
            }
            // scale tmp by coef and subtract from C: C = C - coef*tmp
            for (int i = 0; i < DEG; i++) {
                if (tmp.x[i] == 0) continue;
                // C.x[i] = (C.x[i] - coef*tmp.x[i]) mod N
                int v = (C.x[i] - coef * tmp.x[i]) % N;
                if (v < 0) v += N;
                C.x[i] = v;
            }

            if (2 * L <= n) {
                // B = T, b = d, L = n+1-L, m = 1
                memcpy(&B, &Tt, sizeof(vec));
                b = d % N;
                int newL = n + 1 - L;
                L = newL;
                m = 1;
            } else {
                m++;
            }
        }
    } // end for n

    // 結果: C が Λ(x) の係数を持つ（定数項 C.x[0] = 1）
    // 正常化: 係数を mod N で整える
    for (int i = 0; i < DEG; i++) if (C.x[i]) C.x[i] = (C.x[i] % N + N) % N;

    // 出力
    memcpy(out_f, &C, sizeof(vec));
    return L; // Λ の次数
}

ymo bm_itr(unsigned int s[])
{
    vec U1[2][2] = {0}, U2[2][2][2] = {0}, null = {0};
    int i, j, k;
    ymo t = {0};

    U2[0][0][0].x[0] = 1;       // f[0];
    U2[0][0][1].x[0] = 0;       // fai[0];
    U2[0][1][0].x[0] = 0;       // g[0];
    U2[0][1][1].x[0] = N - (1); // thi[0];
    int m = 0, d = 0, p = 2 * d - m - 1, myu = 0;
    printf("m=%d d=%d myu=%d p=%d\n", m, d, myu, p);
    for (m = 0; m < K; m++)
    {
        d = deg(U2[0][0][0]);
        p = 2 * d - m - 1;
        myu = 0;
        for (int i = 0; i <= d; i++)
            myu = (myu + U2[0][0][0].x[i] * s[i + (m - d)]) % N;

        printf("m=%d ad=%d myu=%d p=%d\n", m, d, myu, p);
        memset(U1, 0, sizeof(U1));
        if (myu == 0 || p >= 0)
        {
            U1[0][0].x[0] = 1;
            U1[0][1].x[p] = N - (myu);
            U1[1][0].x[0] = 0;
            U1[1][1].x[0] = 1;
            // exit(1);
        }
        else if (myu > 0 && p < 0)
        {
            if (p < 0)
            {
                p = (-1 * (p));
            }
            U1[0][0].x[p] = 1;
            U1[0][1].x[0] = N - (myu);
            U1[1][0].x[0] = inv(myu, N);
            U1[1][1].x[0] = 0;
        }
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                printpoln(U1[i][k]);
                for (int k = 0; k < 2; k++){
                    printf("ii! %d %d %d %d\n",i,k, deg(U1[i][k]), deg(U2[0][k][j]));
                    if(deg(U1[0][0])>N-K){
                        printf("blake\n");
                    for(int ii=N-K;ii<DEG;ii++)
                        U1[0][0].x[ii]=0;
                    }
                    U2[1][i][j] = (vadd((U2[1][i][j]), (vmul(U1[i][k], U2[0][k][j],N))));
                }   if(deg(U2[1][i][j])>DEG){
                    printf("breake\n");
                        //exit(1);
                }
            }
        }
        memcpy(U2[0], U2[1], sizeof(U2[0]));
        memset(U2[1], 0, sizeof(U2[1]));
    }
    t.f = U2[0][0][0];
    t.g = U2[0][1][0];
    t.h = U2[0][0][1];
    if (deg(t.f) == T)
    {
        printsage((t.f));
        printf(" ==chen00\n");
        return t;
    }
    else
    {
        t.f = U2[1][0][0];
        printsage((t.f));
        printf("baka\n");
        exit(1);
    }
}


ymo bm_itr2( int s[])
{
    vec U1[2][2] = {0}, U2[2][2][2] = {0}, null = {0};
    int i, j, k;
    ymo t = {0};

    U2[0][0][0].x[0] = 1;       // f[0];
    U2[0][0][1].x[0] = 0;       // fai[0];
    U2[0][1][0].x[0] = 0;       // g[0];
    U2[0][1][1].x[0] = N - (1); // thi[0];
    int m = 0, d = 0, p = (2 * d - m - 1)%N, myu = 0;
    printf("m=%d d=%d myu=%d p=%d\n", m, d, myu, p);
    for (m = 0; m < K; m++)
    {
        d = deg(U2[0][0][0]);
        p = (2 * d - m - 1)%N;
        myu = 0;
        for (int i = 0; i <= d; i++)
            myu = (myu + U2[0][0][0].x[i] * s[i + (m - d)]) % N;

        printf("m=%d ad=%d myu=%d p=%d\n", m, d, myu, p);
        memset(U1, 0, sizeof(U1));
        if (myu == 0 || p >= 0)
        {
            U1[0][0].x[0] = 1;
            U1[0][1].x[p] = N - (myu);
            U1[1][0].x[0] = 0;
            U1[1][1].x[0] = 1;
            // exit(1);
        }
        else if (myu > 0 && p < 0)
        {
            if (p < 0)
            {
                p = -1 * (p);
            }
            U1[0][0].x[p] = 1;
            U1[0][1].x[0] = N - (myu);
            U1[1][0].x[0] = inv(myu, N);
            U1[1][1].x[0] = 0;
        }
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 2; k++)
                    U2[1][i][j] = (vadd((U2[1][i][j]), (vmul(U1[i][k], U2[0][k][j],N))));
            }
        }
        memcpy(U2[0], U2[1], sizeof(U2[0]));
        memset(U2[1], 0, sizeof(U2[1]));
    }
    t.f = U2[0][0][0];
    t.g = U2[0][1][0];
    t.h = U2[0][0][1];
    if (deg(t.f) == T)
    {
        printsage((t.f));
        printf(" ==chen00\n");
        return t;
    }
    else
    {
        t.f = U2[1][0][0];
        printsage((t.f));
        printf("baka\n");
        exit(1);
    }
}

// 行列の掛け算関数
void matrix_multiply(int A[MATRIX_SIZE][MATRIX_SIZE], int B[MATRIX_SIZE][MATRIX_SIZE], int *C, int start_row, int end_row)
{
    for (int i = start_row; i < end_row; i++)
    {
        for (int j = 0; j < MATRIX_SIZE; j++)
        {
            int sum = 0.0;
            for (int k = 0; k < MATRIX_SIZE; k++)
            {
                sum += A[i][k] * B[k][j];
            }
            C[i * MATRIX_SIZE + j] = sum;
        }
    }
}

int ipow(int b, int n)
{
    int l = 1;

    if (n == 0)
        return 1;

    for (int i = 0; i < n; i++)
        l = b * l % N;

    return l;
}


vec recv(MTX t, vec v)
{
    vec x = {0};
    int i, j, k;
    for (i = 0; i < K/2; i++)
    {
        for (k = 0; k < K/2; k++)
            x.x[i] += (v.x[k] * t.x[i][k])%N;
        x.x[i] %= N;
        if(x.x[i]<0)
        x.x[i]=(x.x[i]+N)%N;
    }

    return x;
}

vec ev(vec x,vec v)
{
    int i, j, k;
    MTX mmk = {0};
    MTX inv_A = {0};
    vec tx = {0};

    for (i = 0; i < K / 2; i++)
    {
        for (int j = 0; j < K / 2; j++)
        {
            mmk.x[j][i] = mat[x.x[i]][j];
            printf("%d %df", mat[x.x[j]][i], x.x[j]);
        }
        printf("\n");
    }
    printf("\n(");
    for (i = 0; i < K / 2; i++)
        printf("%d ", x.x[i]);
    printf(")\n");
    // exit(1);
    //mmk.x[0][K / 2] = 2;
    //mmk.x[1][K / 2] = 5;
    for (i = 0; i < K / 2; i++)
    {
        mmk.x[i][K / 2] = v.x[i];
        for (int j = 0; j < K / 2; j++)
            printf("%d^ ", mmk.x[i][j]);
        printf("\n");
    }

    // tx.x[1]=v.x[1];
    for (i = 0; i < K / 2; i++)
        tx.x[i] = v.x[i];
        //tx.x[1] = 5; //v.x[i];
    // inv_A=ver(mmk);
    inv_A = inverseMatrix(mmk, inv_A, 0, K / 2);
    printf("inv %d %d %d\n", inv_A.x[0][0], inv_A.x[0][1], inv_A.x[0][2]);
    printf("inv %d %d %d\n", inv_A.x[1][0], inv_A.x[1][1], inv_A.x[1][2]);
    tx = recv(inv_A, tx);
    for (i = 0; i < K / 2; i++)
        printf("error value is %d\n", tx.x[i]);

        return tx;
}

//モニック多項式にする
vec coeff2(vec f,  int d,int R)
{
  int i, j, k;
  vec a, b;

  k = deg((f)) + 1;
  for (i = 0; i < k; i++)
    f.x[i] = (f.x[i]*oinv2(d,R))%R;

  return f;
}


//多項式の商を取る
vec vdiv2(vec f, vec g,int R)
{

  int i = 0, j, n, k;
  vec h = {0}, e = {0}, tt = {0};
  oterm a, b = {0}, c = {0};

  if (vLT(f).n == 0 && vLT(g).a == 0)
  {
    printf("baka^\n");
    //return f;
    exit(1);
  }
  if (vLT(g).a == 0)
  {
    exit(1);
  }
  if (vLT(g).n == 0 && vLT(g).a > 1)
    return coeff2(f, vLT(g).a,R);

  k = deg(g);
  b = vLT(g);
  if (b.a == 1 && b.n == 0)
    return f;
  if (b.a == 0 && b.n == 0)
  {
    printf("baka in vdiv\n");
    exit(1);
  }
  if (deg((f)) < deg((g)))
  {
    return f;
    //  a=LT(f);
  }

  printpol(f);
  printf(" ==f\n");
  printpol(g);
  printf(" ==g\n");
  i = 0;
  while (vLT(f).n > 0 && vLT(g).n > 0)
  {
    c = vLTdiv2(f, b,R);
    assert(c.n < DEG);
    tt.x[c.n] = c.a%R;
    //i++;

    h = vterml2(g, c,R);
    printpol(h);
    printf(" ==h\n");
    f = vsub2(f, h,R);
    printpol(f);
    printf(" ==f\n");
    if (deg((f)) == 0 || deg((g)) == 0)
    {
      //printf ("blake2\n");
      break;
    }

    if (c.n == 0)
      break;
  }

  // tt は逆順に入ってるので入れ替える
  return tt;
}

/**
 *  f = f << 1
 */
//pieko
vec shiftRotateL(vec f) {
    int i;
    int fn = f.x[deg(f) - 1];
    for (int i=deg(f)-1; i>0; i--)
        f.x[i] = f.x[i-1];
    f.x[0] = fn;

    return f;
}

/**
 *  f = f >> 1
 */
 //pieko
vec shiftRotateR(vec f) {
    int i;
    int f0 = f.x[0];
    for (int i=0; i < deg(f) ; i++)
        f.x[i] = f.x[i+1];
    f.x[deg(f) - 1] = f0;

    return f;
}

int mod(int x, int R){
    return ((x % R) + R) % R;
}


// invert of polynomial
vec vinv(vec a, vec n)
{
    vec d = n;
    vec x = {0};
    vec s = {0};
    vec t={0},r={0};
    vec tt=n;

    s.x[0]=1;

    while (deg(a)>0)
    {
        vec q = vdiv(d , a);
        r = vmod(d , a);
        d = a;
        a = r;
        t = vsub(x, vmul(q, s,N));
        x = s;
        s = t;
    }
    d = a;
    a = r;
    x = s;
    s = t;

    vec gcd = d; // $\gcd(a, n)$
    vec u=vmod(vadd(x , n),tt);
    u=vdiv(u, d);

    return  u;
}

// invert of polynomial mod R
vec vinv2(vec a, vec n,int R)
{
    vec d = n;
    vec x = {0};
    vec s = {0};
    vec t={0},r={0};
    vec tt=n;

    s.x[0]=1;

    while (deg(a)>0)
    {
        vec q = vdiv2(d , a,R);
        r = vmod2(d , a,R);
        d = a;
        a = r;
        t = vsub2(x, vmul(q, s,R),R);
        x = s;
        s = t;
    }
    d = a;
    a = r;
    x = s;
    s = t;

    vec gcd = d; // $\gcd(a, n)$
    vec u=vmod2(vadd2(x , n,R),tt,R);
    u=vdiv2(u, d,R);

    return  u;
}


vec trim(vec a,int R){
    int i;
    vec v={0};

    printpol(a);
    printf(" ==a\n");
    for(i=0;i<DEG;i++){
    {
    while(a.x[i]<0)
    a.x[i]+=R;
    }
    }
    for(i=0;i<DEG;i++){
    v.x[i]=a.x[i]%R;
    }

    return v;
}

vec pes(vec u,vec v){
vec g={0};

while(deg(v)>0){
g=u;
u=vdiv2(v,u,32);
v=g;
}
return u;
}

vec deli(vec a, vec b)
{

  vec v = {0};

  for (int i = 0; i < deg(b); i++)
    v.x[i] = a.x[i];

  return v;
}

vec vcoef(vec v)
{
  unsigned int n=0, k = deg(v);

  // if(v.x[0]==0)
  // return v;

  if (v.x[0] > 1)
    n = inv(v.x[0],N);
  for (int i = 0; i < k + 1; i++)
    v.x[i] = n*v.x[i]%N;

  return v;
}

// ニュートン法で逆元を求める (something buggy)
vec invpol(vec a)
{
  vec v = {0}, x = {0},g={0},f={0};
  int i;

  f.x[0]=2;
  x.x[2] = 1;
  //a = mon(oinv(trace(a,0),N),a);
  //g.x[0] = (trace(a,0));
  g.x[0] = 1;

  //if (a.x[0] > 1)
  //  a = vcoef(a);

  i = 1;
  while (i < K+1)
  {
    g=vsub(vmul(g,f,N),vmul(a,vmul(g,g,N),N));
    //v = vmul(vmul(v, v), a);
    if (i > 1)
      x = vmul(x, x,N);
    g = deli(g, x);
    i*=2;
    printpol(g);
    printf(" ('A`)\n");
  }
  if(vmul(a,g,N).x[0]!=1)
  {
    printf("baka inv\n");
    //exit(1);
    }
    //g=mon(equ(39,37),g);

  return g;
}


vec invpol2(vec a,vec I,int R)
{
  vec v = {0}, x = {0},g={0},f={0};
  int i;

  f.x[0]=2;
  x.x[2] = 1;
  //a = mon(oinv(trace(a,0),N),a);
  //g.x[0] = (trace(a,0));
  g.x[0] = 1;

  //if (a.x[0] > 1)
  //  a = vcoef(a);

  i = 1;
  while (i < K+1)
  {
    g=(vsub2(vmul(g,f,R),vmul(a,vmul(g,g,R),R),R));
    //v = vmul(vmul(v, v), a);
    if (i > 1)
      x = vmul(x, x,R);
    printpol(g);
    printf(" ----b\n");
    g = deli(g, x);
    //g=vmod2(g,I,32);
    i*=2;
    printpol(g);
    printf(" ('A`)\n");
    if(vLT(g).a==30)
    exit(1);
  }
  if(vmul(a,g,R).x[0]!=1)
  {
    printf("baka inv2\n");
    exit(1);
    }
    //g=mon(equ(39,37),g);

  return g;
}

vec ev3(vec x,vec v)
{
    int i, j, k;
    MTX mmk = {0};
    MTX inv_A = {0};
    vec tx = {0};

    for (i = 0; i < K / 2; i++)
    {
        for (int j = 0; j < K / 2; j++)
        {
            mmk.x[j][i] = mat[x.x[i]][j];
            printf("%d %df", mat[x.x[j]][i], x.x[j]);
        }
        printf("\n");
    }
    printf("\n(");
    for (i = 0; i < K / 2; i++)
        printf("%d ", x.x[i]);
    printf(")\n");
    // exit(1);
    //mmk.x[0][K / 2] = 2;
    //mmk.x[1][K / 2] = 5;
    for (i = 0; i < K / 2; i++)
    {
        mmk.x[i][K / 2] = v.x[i];
        for (int j = 0; j < K / 2; j++)
            printf("%d^ ", mmk.x[i][j]);
        printf("\n");
    }

    // tx.x[1]=v.x[1];
    for (i = 0; i < K / 2; i++)
        tx.x[i] = v.x[i];
        //tx.x[1] = 5; //v.x[i];
    // inv_A=ver(mmk);
    inv_A = inverseMatrix(mmk, inv_A, 0, K / 2);
    printf("inv %d %d %d\n", inv_A.x[0][0], inv_A.x[0][1], inv_A.x[0][2]);
    printf("inv %d %d %d\n", inv_A.x[1][0], inv_A.x[1][1], inv_A.x[1][2]);
    tx = recv(inv_A, tx);
    for (i = 0; i < K / 2; i++)
        printf("error value is %d\n", tx.x[i]);

        return tx;
}


typedef struct {
    vec g;
    vec h;
} fair;

fair soka(){
    int i;
    vec g[K]={0},h[N-K]={0},mod={0},gg={0},hh={0};
    fair kubi={0};

    gg.x[0]=1;
    hh.x[0]=1;
    for(i=1;i<N-K;i++){
    g[i-1].x[0]=i;
    g[i-1].x[1]=1;
    }
    for(i=N-K;i<N-1;i++){
    h[i-N+K].x[0]=i;
    h[i-N+K].x[1]=1;
    }
    for(i=0;i<K;i++)
    gg=vmul(gg,g[i],N);
    for(i=0;i<N-K-1;i++)
    hh=vmul(hh,h[i],N);
    mod.x[0]=N-1;
    mod.x[N-1]=1;
    printpoln(vmul(gg,hh,N));
    printpoln(vmod(vmul(gg,hh,N),mod));
    kubi.g=hh;
    kubi.h=gg;
    //exit(1);

return kubi;
}

vec msm(vec err){
int i,l;
vec syn[K]={0};
vec sin={0};

printpoln(err);
for(int j=1;j<M;j++)
{
l=0;
    for(i=1;i<K+1;i++)
    {
        if(err.x[j]>0)
        {
        l=mltn(i,j);
            sin.x[i-1]+=l;
            sin.x[i-1]%=N;
        //sin.x[i-1]=trace(err,l);
        }
    }
    printf("l=%d %d %d\n",sin.x[j-1],j,l%N);
    }
        printf("Uh!\n");
        printpoln(sin);
        ymo y=bm_itr(sin.x);
        chen(y.f);

        return sin;
}


vec keygen(){
    int i;
    vec g0={0},gg[K+1]={0};

    for(i=0;i<N-K;i++){
    gg[i].x[0]=N-mltn(i,Prim);
    gg[i].x[1]=1;
    printf("%d\n",trace(gg[i],i));
    }
    for(i=0;i<4;i++)
    printf("%d,",trace(gg[0],i));
    printf("\n");
    //exit(1);

    g0.x[0]=1;
    for(i=0;i<N-K;i++)
    g0=vmul(g0,gg[i],N);
    printpoln(g0);
    for(i=1;i<K+1;i++)
    printf("%d,",trace(g0,i));
    printf("\n");
    //exit(1);

    return g0;
}


/**
 * struct state - represents the 320-bit state of ascon
 *
 * @x: array containing the five 64-bit registers of the state
 */
struct state {
	uint64_t x[5];
};

static inline unsigned rotr(uint32_t x, uint8_t n)
{
	return x >> n | x << (32 - n);
}



int vor(vec v){
    for(int i=0;i<N;i++){
    if(v.x[256]>0){
        printf("baka\n");
        exit(1);
    }
    }
}

uni coda(uni on,vec a1,vec a2){
    uni hola={0};
    int i;

    /*
    for(i=0;i<N-1;i++)
    printf("%d=%d\n",i,ink(mltn(i,3)));
    printf("%d\n",v2i(i2v(rotr(1234,17))));
    //exit(1);
    */
    //for(i=0;i<K;i++)
    //    on.c[i]=i+1;
    for(i=0;i<K;i++)
    on.d[i]=rotr(on.d[i],17);
    //vec a[8]={0},inv_a[8]={0};
    //vec aa[30]={0},o[30]={0};

    //exit(1);
    /*
    vec ass={0},them={0};
    for(i=0;i<32;i++){
        ass.x[i]=(1<<a.x[i]);
        them.x[i]=(1<<inv_a.x[i]);
    }
    */
    uni und={0};
    for(i=0;i<K*4;i++)
        und.c[i]=on.c[a1.x[i]];
    for(i=0;i<K;i++)
    printf("z%d,",und.d[i]);
    printf("\n");
    for(i=0;i<K;i++)
    printf("%b\n,",und.d[i]);
    printf("\n");
    uni dd={0};
    uni tar={0};
    for(i=0;i<K;i++){
            vec v=i2v(und.d[i]);
            vec vv={0};
        for(int j=0;j<32;j++){
            vv.x[j]=v.x[a2.x[j]];
            //tar.d[i]^=und.d[i]&them.x[j];
        }
        tar.d[i]=v2i(vv);
        //tar.d[i]=rotr(tar.d[i],17);
    }
    //exit(1);

return tar;
}


uni koda(uni und,vec inv_a0,vec inv_a1){
    uni hola={0};
    int i;

    /*
    for(i=0;i<N-1;i++)
    printf("%d=%d\n",i,ink(mltn(i,3)));
    printf("%d\n",v2i(i2v(rotr(1234,17))));
    //exit(1);
    */
    for(i=0;i<K;i++)
    printf("%b\n",und.d[i]);
    //exit(1);

    uni dd={0};
    uni tar={0};

    //for(i=0;K/4;i++)
    //tar.d[i]=rotr(tar.d[i],15);

    for(i=0;i<K;i++){
            vec v=i2v(und.d[i]);
            vec vv={0};
        for(int j=0;j<32;j++){
            vv.x[j]=v.x[inv_a1.x[j]];
            //tar.d[i]^=und.d[i]&them.x[j];
        }
        tar.d[i]=v2i(vv);
    }
    for(i=0;i<K*4;i++)
    dd.c[i]=tar.c[inv_a0.x[i]];
    for(i=0;i<K;i++){
    dd.d[i]=rotr(dd.d[i],15);
    }
    printf("\n");

    return dd;
}


vec zind(vec e){
int i;
vec sin={0};

printf("blahut\n");
for(i=1;i<K+1;i++)
    sin.x[i-1]=trace(e,mltn(i,Prim));

    return sin;
}

vec zynd(vec e) {
    vec sin = {0};
    int t = (N - K) / 2;   // 訂正能力
    int Slen = 2 * t;

    for (int i = 1; i <= Slen; i++) {
        sin.x[i-1] = trace(e, mltn(i, Prim));
    }
    return sin;
}


vec L2(){
    vec v={0};
    int i,count=0,l=0;
    vec L={0};
    vec f={0},ff[K]={0};
    vec g0={0};

    while(1){
        l=rand()%N;
        if(L.x[l]==0 && l>0){
            f.x[count++]=l;
            L.x[l]=1;
        }
        if(count==K)
        break;
    }

    return f;
}


vec L3(vec f){
    vec v={0};
    int i,count=0,l=0;
    vec L={0};
    vec ff[K]={0};
    vec g0={0};


for (count = 0; count < K; ++count) {
    int alpha_pow = mltn(count+1, Prim); // alpha^(count+1)
    // (x - alpha_pow) -> constant term = (q - alpha_pow) % q
    ff[count].x[0] = (N - alpha_pow) % N;
    ff[count].x[1] = 1;
}
g0.x[0] = 1;
for (i = 0; i < K; ++i) g0 = vmul(g0, ff[i], N);

/*
    while(1){
        ff[count].x[0]=N-mltn(f.x[count],Prim);
        ff[count++].x[1]=1;

        if(count==N-K)
        break;
    }
    g0.x[0]=1;
    for(i=0;i<N-K;i++)
    g0=vmul(g0,ff[i],N);
*/
    return g0;
}





void cipher(){
    int i,j;
    uni on={0};
    vec mv={0},us={0};
    vec g0={0},a[2]={0},inv_a[2]={0};
    unsigned gol=0b101011100011;

    vec L={0};
    //L=L2(); //keygen();
    for(i=0;i<K/2;i++)
    L.x[i]=i+1;
    g0=L3(L);


    for(i=0;i<K*4;i++)
    a[0].x[i]=i;
    for(i=0;i<32;i++)
    a[1].x[i]=i;

    random_shuffle(a[0].x,K*4);
    random_shuffle(a[1].x,32);

    for(i=0;i<K*4;i++){
        printf("%d,",a[0].x[i]);
    inv_a[0].x[a[0].x[i]]=i;
    }
    for(i=0;i<32;i++)
    inv_a[1].x[a[1].x[i]]=i;


    printf("ooky\n");
    for(i=0;i<K/2;i++)
    mv.x[i]=on.d[i]=i+1;
    mv=vmul(mv,g0,N);
    printpoln(mv);
    //exit(1);
    //for(i=0;i<K;i++)
    //on.d[i]=mm.x[i]; //rotr(on.d[i],17);
    for(i=0;i<K;i++)
    on.d[i]=mv.x[i]; //=v2i(vdiv(i2v(mm.x[i]),i2v(gol)));
    on=coda(on,a[0],a[1]);

    vec cc={0};
    //for(i=0;i<K*4;i++)
    //cc.x[i]=m(on.c[i],gol);
    printf("y&t=");
    printpoln(cc);
    //exit(1);
    for(i=0;i<K*4;i++)
    on.c[i]=v2i(bdiv(i2v(cc.x[i]),i2v(gol)));
    printf("\n");
    //printpoln(us);
    //int y=m(0b11111111,gol);
    //int p=v2i(bdiv(i2v(y),i2v(gol)));
    printf("%b \n",p);
    //exit(1);

    on=koda(on,inv_a[0],inv_a[1]);
    printf("\n");
    //for(i=0;i<K;i++)
    //us.x[i]=v2i(vdiv(i2v(on.d[i]),i2v(gol)));

    for(i=0;i<K*2;i++)
    us.x[i]=on.d[i];
    us=vdiv(us,g0);
    printpoln(us);
    for(i=0;i<K/2;i++)
    printf("u%d,",us.x[i]);
    printf("\n");
    //on.d[i]=rotr(on.d[i],15);
    //fugo();
    //exit(1);
}


vec itiji(vec c[]){
int i,j,x[N+1]={0},count=0;
vec a={0};

i=0;
x[0]=1;
while(1){
  j=0;
  count=1;
while(1){
    int d=deg(c[N-i-1]);
if(c[N-i-1].x[j]>0)
{
  printf("Ah=%d %d %d\n",c[N-i-1].x[j],j,d);
x[count]=(x[count]+N-(c[N-i-1].x[j]*x[j]))%N;
//x[count]=c[N-i-1].x[d]*x[count];
count++;
}
j++;
//printf("j=%d\n",j);
if(j==d)
break;
}

i++;
if(i>N-1)
break;
}
for(i=0;i<5;i++){
printf("%d=%d\n",i,x[i]);
}

}


vec printrev(vec a,int count){
int i,j=deg(a);
vec x={0};

for(i=0;i<N+1;i++){
    if(a.x[i]!=0){
//printf("%dQ^%d+",a.x[i],j-i);
x.x[j-i]=a.x[i];
    }
}
//printf("\n");

return x;
}

#include "dmd.c"
vec dmd(MTX bb){
    int i,j;

    oterm xa[N+1][N+1]={0};
    xa[0][0].a=1;
    xa[0][0].n=1;
    for(i=0;i<N;i++){
      for(int j=0;j<N+1;j++){
      xa[N-1-i][j].a=bb.z[i].x[N-j];
      //printf("%d=x[%d]\n",b.z[i].x[j],N-j);
      }
    }
    for(i=0;i<N;i++){
      for(int j=0;j<N+1;j++){
      printf("%dx[%d]=%d\n",xa[i][j].a,j,i);
      }
    }
    int yx[N+1]={0};
    yx[0]=1;
    xa[0][0].n=1;
    //x[0][0].a=1;
    for(i=1;i<N+1;i++){
      for(j=0;j<i;j++){
        xa[0][i].n+=xa[i-1][j].a*xa[0][j].n;
        printf("y[%d]=%d*%d\n",i-1,xa[i-1][j].a,xa[0][j].n);
      }
        xa[0][i].n=yx[i]=N-xa[0][i].n%N;
    }
    for(i=0;i<N+1;i++)
    printf("%d,%d[%d]=%d\n",yx[i]%N,xa[0][i].a,i,xa[0][i].n);
    printf("\n");

    vec ucc={0};
    for(i=0;i<T+1;i++){
    ucc.x[T-i]=yx[i]%N;
    printf("%d,",ucc.x[T-i]);
    }
    printf("\n");
    vec uec={0};
    for(i=T+1;i<N+1;i++){
    uec.x[N-i]=yx[i]%N;
    printf("%d,",uec.x[N-i]);
    }
    printf("\n");
    printpoln(uec);
    printpoln(ucc);
    printf("\n");
    //exit(1);

    printf("%d\n",deg(vdiv(uec,ucc)));
    vec b=vdiv(uec,ucc);
    printf("timole1= ");
    printpoln(b);
    vec a=vmod(uec,ucc);
    printf("timole2= ");
    printpoln(a);
    if(deg(a)>0){
        printf("baka\n");
        exit(1);
    }
    vec ans={0};
    //if(deg(a)==0 && vLT(a).a>0)
    //b=kof2(vLT(a).a,b);
    
    for(i=0;i<N;i++){
    ans.x[i]=((N-trace(b,(i+1)%N))%N);
    printf("sy[%d],",ans.x[i]);
    }
    printf("\n");
    
    //printpoln(ans);
    //exit(1);

    return ans;
}

vec vecky(vec mm){
int i,j,k;
vec c={0};

for(j=0;j<N;j++){
    for(k=0;k<K;k++)
        c.x[j]+=mm.x[k]*vb[k][j];
        c.x[j]%=N;
        printf("%d,",c.x[j]);
}
    printf("\n veck= ");;
    printpoln(c);

return c;
}


vec konv(vec mm){
    vec f={0};
    int i,j;

    for(i=0;i<N;i++){
        printf("tr%d,",(trace(mm,((i+1)%N)))%N);
        f.x[i]=(trace(mm,((i+1)%N)))%N;
    }
    printf("\n konv= ");
    printpoln(f);
    //exit(1);

    return f;
}


vec vmul_2(vec a, vec b)
{
    int i, j, k, l;
    vec c = {0};

    k = deg(a);
    l = deg(b);

    if(l+k>N*2){
        printf("blake %d a=%d b=%d\n",l+k,deg(a),deg(b));
        //exit(1);
    }
    i = 0;
    while (i < k + 1)
    {
        for (j = 0; j < l + 1; j++)
        {
            if (a.x[i] > 0)
                c.x[i + j] = (c.x[i + j] + a.x[i] * b.x[j]);
        }
        i++;
    }

    return c;
}


unsigned int to_m(unsigned int y, unsigned int z)
{
  unsigned int c;
  
  c=0;
  while(y!=0){
    if(y&1) c ^=z;
    z<<=1; y>>=1;
  }
  return c;
  
}

void naoya(void){

    int i,count=0,cnt=0;
    vec mm={0},f={0};
    int zk[8192]={0};

    while(1){
        cnt=rand()%N;
        if(zk[cnt]==0){
        zk[cnt]=1;
        count++;
        }
        if(count>86)
        break;
    }
    // s=129
    for(i=0;i<42;i++)
    mm.x[i]=rand()%N;
    printf("mm= ");
    printpoln(mm);
    //s+1=130
    vec r=mkpol2(42);
    vec ma=vadd(r,mm);
    //s+d=193
    vec rr=mkpol2(85);
    vec cd[256]={0};
    //s+d+1=194
    for(i=0;i<86;i++){
    cd[i].x[1]=1;
    //cd[i].x[0]=rand()%N;
    }
    count=0;
    for(i=0;i<8192;i++){
    if(zk[i]>0){
        cd[count].x[0]=i;
        count++;
    }
    if(count==86)
    break;
    }
    
    vec ed={0};
    ed.x[0]=1;
    //d=63
    vec ra=mkpol2(43);
    for(i=0;i<86;i++)
    ed=vmul(ed,cd[i],N);
    //s+d+1=194
    vec pk=vadd(ed,ra);
    vec ca=vadd(vmul(pk,ma,N),rr);
    for(i=0;i<86;i++)
    f.x[i]=trace(ca,-(cd[i].x[0]));
    vec num={0};
    vec out={0};
    for(i=0;i<86;i++)
    num.x[i]=cd[i].x[0];
    lagrange_interpolate(&num.x, &f.x,40,&out.x);
    vec md=vsub(ca,out);
    vec mc=vdiv(md,ed);
    vec a=vsub(mc,r);
    printf("plain=\n");
    printpoln(a);
    printf("original=\n");
    printpoln(mm);
    printf("\n\n pubkey=\n");
    printpoln(pk);
    printf("cipher=\n");
    printpoln(ca);
    printf("secret point of eval==\n");
    printpoln(ed);
    exit(1);
    //f.x[i]=trace(mm,i+1);

}

int main()
{
    int i,j, u = 0;
    int s[K + 1] = {0}, z1[N] = {0};
    fair ff20={0};
    unsigned gol=0b101011100011;

    vec x={0};
    vec us={0},mm={0},cc={0};
    vec gg[K+1]={0},g0={0};
    int nn=3;
    int jj=1;
    uni on={0};

    srand(clock());

    for(i=0;i<K;i++)
    on.x[i]=i+1;


    vec dd[100]={0};

    x=mkpol2(10);
    vec rr=mkpol2(10);
    for(i=0;i<10;i++){
    dd[i].x[1]=1;
    dd[i].x[0]= (i+1); //rand()%N;
    }
    vec l={0};
    l.x[0]=1;
    for(i=0;i<10;i++)
    l=vmul(l,dd[i],N);
    vec pk=vadd(vmul(x,rr,N),l);
    printpoln(pk);
    printpoln(l);
    for(i=0;i<10;i++)
    printf("%d,",trace(l,-(i+1)));
    printf("\n");
    vec m={0};
    for(int i=0;i<10;i++)
    m.x[i]=i+1;
    
    vec d2[100]={0};
    for(int i=0;i<10;i++){
    d2[i].x[1]=1;
    d2[i].x[0]=i+10;
    }
    vec lcl={0};
    lcl.x[0]=1;
    for(i=0;i<10;i++)
    lcl=vmul(lcl,d2[i],N);
    vec lr=mkpol2(10);
    vec cb=vadd(vmul(pk,vadd(lcl,m),N),lr);
    printpoln(cb);


    //exit(1);


/*
MTX A={0},inv_A={0};
    for(i=0;i<N;i++){
        for(int j=0;j<N+1;j++)
        A.x[i][j]=rand()%N;
        }
MTA B={0};
    for(i=0;i<K;i++){
        for(int j=0;j<K+1;j++)
        B.x[i][j]=rand()%N;
    }

    matinv(A,&inv_A,K);
    //sankaku(A,N);
    //renritu(monde(N),N-1);
    //cipher();
    exit(1);
    van(K);
    //vv(K);
    vec v={0};
    //for(i=0;i<N;i++)
    v.x[3]=1;
    for(i=0;i<N;i++)
    printf("%d,",trace(v,mltn(i,Prim)));
    printf("\n");
    //list(v);
    //exit(1);
    //for(i=0;i<K-1;i++)
    //mm.x[i]=rand()%2;
    */
    int y=to_m(0b11111111,gol);
    int p=v2i(bdiv(i2v(y),i2v(gol)));
    printf("%b \n",p);
    //exit(1);
    vec vc={0}; //vmul(mm,v,N);
    //vc=vadd(vc,vmod(vc,g0));

    //mm.x[0]=rand()%N;
    //mm.x[1]=rand()%N;
    vec L={0};
    for(i=0;i<K;i++)
    L.x[i]=i+1;
    g0=L3(L); //
    //keygen();
    //exit(1);
    
    printf("g0=\n");
    printpoln(g0);
    printf("code=\n");
    printpoln(vmul(mm,g0,N));
    printf("m=\n");
    printpoln(mm);
    //(vecky(mm));
    //(konv(mm));
    //exit(1);

    naoya();
    exit(1);

    /*
    for(i=0;i<N;i++){
        vc.x[i]=0;
        for(j=0;j<K;j++){
            vc.x[i]+=mm.x[j]*mat[j][i]%N;
        printf("%d,",mat[i][j]);
        }
        vc.x[i]%=N;
    printf("\n");
    }


    printpoln(vc);
    //printpoln(mm);
    //exit(1);
    */


    for(i=0;i<129;i++)
    mm.x[i]=i+1;
    vec cv=vmul(mm,g0,N); //vecky(mm);
    printf("符号語\n");
    printpoln(cv);
    printf("平文\n");
    printpoln(mm);
    for(i=0;i<N;i++){
        printf("trf %d,",trace(mm,(i+1)%N));
    }
    printf("\n");
    //exit(1);

    vec t={0};
    mkerr(t.x,T);

    //for(i=0;i<N;i++)
    //t.x[i]=1;
    vec f={0}; //konv(mm);
    //exit(1);
    vec vx=vadd(f,t);
    //list(vx);
    /*
    //exit(1);
    for(i=0;i<N;i++){
        for(int j=0;j<N+1;j++){
        A.x[i][j]=mat[i][j];
        printf("%d,",mat[i][j]);
    }
    printf("\n");
    }
    //exit(1);
    
    MTX bb=sankaku(A,N);
    //exit(1);

    vec ans=dmd2(bb);
    printf("ans= ");
    printpoln(ans);
    
    printf("f+e= ");
    printpoln(vx);
    printf("e= ");
    printpoln(t);
    printf("err= ");
    printpoln(vsub(vx,ans));
    */


    //cv.x[0]+=1;
    
    
    vc=vadd(vc,vmod(vc,g0));
    //vmul(mm,g0,N);
    printpoln(mm);
    printf("encode=\n");
    printpoln(vc);
    vec xv=vmul(mm,g0,N);
    printpoln(xv);
    //exit(1);
    
    /*
    vec vv={0};
    for(i=0;i<N;i++)
    vv.x[N-i-1]=vc.x[i];
    //vc=vadd(xv,t);
    //vec zx=monde(mm);
    printf("monde=\n");
    printpoln(vc);
    //printpoln(zx);
    printpoln(xv);
    printpoln(t);
    //list(vc);
    //exit(1);

    vec vb=xv; //vadd(t,vc);
    list(vc);
    for(i=0;i<N;i++){
        for(int j=0;j<N+1;j++){
        A.x[i][j]=mat[i][j];
        printf("%d,",mat[i][j]);
    }
    printf("\n");
    }
    //exit(1);
     bb=sankaku(A,N);

    printpoln(mm);
    dmd(bb);
    printf("code=\n");
    printpoln(vc);
    printpoln(xv);
    printpoln(t);
    //monde(N);
    //exit(1);
    */
    
    //printpoln(vc);
    //exit(1);
    //mkerr(cc.x,T);
    for(i=0;i<T;i++)
    cc.x[i]=1;
    printf("AI\n");
    vec r=vadd(cc,xv); //zind(vadd(cc,vc));
    printpoln(r);
    printf("Ah!\n");
    MTX Z={0}; //sankaku(monde(r,N),N);
    printpoln(vc);
    printf("aa!\n");
    /*
    for(i=0;i<N;i++){
   Z.z[i]=printrev(Z.z[i],i);
   printpoln(Z.z[i]);
    }
    //sankaku(Z,N);
    //itiji(Z.z);
    //exit(1);
    */
    vec ss=zynd(r);
    vec vtmp={0};
    for(i=0;i<N-K;i++)
    vtmp.x[N-K-1-i]=r.x[i];
    ymo o0={0}; //=bm_itr(zind(r).x);
    bm_vec(vtmp.x,N-K,&o0.f);
    chen(o0.f);
    for(i=0;i<N;i++)
    printf("%d,%d\n",i,cc.x[i]);
    printf("\n");
    printpoln(r);
    printf("odamari\n");
    exit(1);

    for(i=0;i<N;i++)
    vc.x[i]=to_m(vc.x[i],gol);
    for(i=0;i<N;i++)
    vc.x[i]=v2i(bdiv(i2v(vc.x[i]),i2v(gol)));
    vc=vdiv(vc,g0);
    printpoln(vc);
    //exit(1);

    vec c=vmul(mm,g0,N);
    printpoln(c);

    vec b={0}; //vmod(c,ff20.h);
    //exit(1);
    for(i=0;i<deg(c);i++)
    b.x[i]=to_m(c.x[i],gol);

    unsigned int P[N]={0},inv_P[N]={0};
    for(i=0;i<N;i++)
    P[i]=i;
    random_shuffle(P,N);
    for(i=0;i<N;i++)
    printf("%d,",P[i]);
    printf("\n");
    //exit(1);



    unsigned plain=0b10000001;
    //printf("%b %b\n",v2i(bdiv(i2v(m(gol,plain)),i2v(gol))),plain);

    printpoln(b);
    printf(" ==bbc\n");
    //exit(1);
    for(i=0;i<N;i++){
    b.x[i]=v2i(bdiv(i2v(b.x[i]),i2v(gol)));
    }
    printpoln(b);
    printf(" =vvc\n");
    b=vdiv(b,g0);
    printpoln(b);
    printpoln(c);
    plain=to_m(plain,gol);
    int geb=v2i(bdiv(i2v(plain),i2v(gol)));
    printf("%b\n",geb);
    //exit(1);

    vec err={0},sin={0},d={0};
    mkerr(err.x,T);
    vec e=vadd(err,c);
    printf("dioscroites=");
    printpoln(e);

    printf("e=");
    printpoln(e);
    printf("b=");
    printpoln(b);
    vec vvc=i2v(gol);
    printpoln(e);

    // trace を使って受信後からシンドロームを計算する
    sin=zind(e);

    printpoln(sin);

    printf("\n");
    printpoln(err);

    vec ea={0};

    ymo yy=bm_itr(sin.x);
    x=chen(yy.f);

    for(i=0;i<N;i++){
        if(x.x[i]>0)
        ea.x[x.x[i]]=1;
    }

    vec e2=vsub(e,(ea));
    for(i=0;i<N;i++){
        if(err.x[i]!=ea.x[i]){
            printf("Ah! %d %d %d\n",i,err.x[i],ea.x[i]);
            exit(1);
        }
    }
    //exit(1);
    for(i=0;i<T;i++)
    printf("%d,",x.x[i]);
    printf("\n\n");
    //exit(1);

    //printpoln(d);
    d=vdiv(e2,g0);

    printpol(d);
    printf(" ==Uh!\n");
    //exit(1);


    //生成行列を使いたいときはこれを使う。主にGoppa符号の時
    van(K);
    //mkd(g0, K);

    return 0;
}
