//学校で習った知識を使ったら問題解決できました。大学行っといてよかったｗ

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define F K *E
#define X (K) * E
#define AY (K / 2 + 1) * E

static MTX inv_S = {0};
static MTX S = {0};
static MTX SS = {0};

extern void makeS();


int is_reg(MTX cc, MTX *R)
{
    int i, j, k, l;
    int retry;
    unsigned char inv_a[K*E][K*E] = {0}; // 逆行列
    unsigned char n = K*E;

    // 元の行列をコピー
    unsigned char A[K*E][K*E];
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            A[i][j] = cc.x[i][j];

    // 単位行列作成
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            inv_a[i][j] = (i==j)?1:0;

    retry = 0;
    for(i=0;i<n;i++)
    {
        // ピボット探索
        if(A[i][i]==0)
        {
            j = i+1;
            while(j<n && A[j][i]==0) j++;

            if(j==n)
            {
                // 0ベクトル or rank不足の行を再生成
                do {
                    //srand(clock() + i + retry);
                    //A[i][i]=1;
                    for(k=0;k<n;k++){
                    //  if(i!=k)
                        A[i][k] = rand()%2;
                    }
                } while (A[i][i]==0); // ピボットが立つまで

                // 逆行列行もリセット
                for(k=0;k<n;k++)
                    inv_a[i][k] = (i==k)?1:0;

                retry++;
                if(retry>10) return -1; // 無限ループ防止
                i--; // 同じ行を再処理
                continue;
            }

            // 行入れ替え
            for(k=0;k<n;k++)
            {
                unsigned char tmp = A[i][k];
                A[i][k] = A[j][k];
                A[j][k] = tmp;

                tmp = inv_a[i][k];
                inv_a[i][k] = inv_a[j][k];
                inv_a[j][k] = tmp;
            }
        }

        // 下三角の消去
        for(l=i+1;l<n;l++)
        {
            if(A[l][i])
            {
                for(k=0;k<n;k++)
                {
                    A[l][k] ^= A[i][k];
                    inv_a[l][k] ^= inv_a[i][k];
                }
            }
        }
    }

    // 上三角の消去（逆行列完成）
    for(i=n-1;i>=0;i--)
    {
        for(k=i-1;k>=0;k--)
        {
            if(A[k][i])
            {
                for(j=0;j<n;j++)
                {
                    A[k][j] ^= A[i][j];
                    inv_a[k][j] ^= inv_a[i][j];
                }
            }
        }
    }

    // 検算
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
        {
            unsigned char sum = 0;
            for(k=0;k<n;k++)
                sum ^= (cc.x[i][k] & inv_a[k][j]);
            if(i==j && sum!=1) return -1;
            if(i!=j && sum!=0) return -1;
        }

    // 成功なら逆行列をコピー
    for(i=0;i<n;i++)
        for(j=0;j<n;j++)
            R->x[i][j] = inv_a[i][j];

    return 0;
}

int is_reg2(MTX cc, MTX *R)
{
    int i, j, k;
    unsigned char inv[F][F] = {0};

    /* 単位行列 */
    for (i = 0; i < F; i++)
        inv[i][i] = 1;

    /* ガウス消去 */
    for (i = 0; i < F; i++)
    {
        if (cc.x[i][i] == 0)
        {
            for (j = i + 1; j < F; j++)
            {
                if (cc.x[j][i])
                    break;
            }
            if (j == F)
                return -1;

            /* 行交換 */
            for (k = 0; k < F; k++)
            {
                unsigned char t;
                t = cc.x[i][k]; cc.x[i][k] = cc.x[j][k]; cc.x[j][k] = t;
                t = inv[i][k];  inv[i][k]  = inv[j][k];  inv[j][k]  = t;
            }
        }

        /* 下消去 */
        for (j = i + 1; j < F; j++)
        {
            if (cc.x[j][i])
            {
                for (k = 0; k < F; k++)
                {
                    cc.x[j][k] ^= cc.x[i][k];
                    inv[j][k]  ^= inv[i][k];
                }
            }
        }
    }

    /* 上消去 */
    for (i = F - 1; i >= 0; i--)
    {
        for (j = 0; j < i; j++)
        {
            if (cc.x[j][i])
            {
                for (k = 0; k < F; k++)
                {
                    cc.x[j][k] ^= cc.x[i][k];
                    inv[j][k]  ^= inv[i][k];
                }
            }
        }
    }

    memcpy(R->x, inv, sizeof(inv));
    return 0;
}

int is_leg(MTX cc, MTX *R)
{

  int i, j, k, l;
  unsigned char b[K * E][K * E] = {0};
  unsigned char dd[K * E] = {0};
  unsigned int flg = 0, count = 0;
  unsigned char cl[K * E][K * E];
  time_t t;
  FILE *fq;
  unsigned char inv_a[K * E][K * E] = {0}; //ここに逆行列が入る
  unsigned char buf;                       //一時的なデータを蓄える
  int n = K * E;                           //配列の次数
  static MTX O = {0};

  while (flg != F)
  {
  labo:
    // memset(cc,0,sizeof(cc));
    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    // g2();

    //#pragma omp parallel for private(j)
    for (i = 0; i < F; i++)
    {

      for (j = 0; j < F; j++)
      {
        // printf("%d,",cc.x[i][j]);
        cl[i][j] = cc.x[i][j];
        dd[j] = cc.x[i][j];
      }
      // printf("\n");
    }

    // memset(inv_a,0,sizeof(inv_a));

    //単位行列を作る
    //#pragma omp parallel for
    for (i = 0; i < F; i++)
    {
      for (j = 0; j < F; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    //掃き出し法

    for (i = 0; i < F; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;

        while (cc.x[j][i] == 0 && j < F)
        {
          j++;
        }

        //#pragma omp parallel for
        if (j > F)
        {
          printf("S is not reg in is_reg %d\n", j);
          // exit(1);
          // cc.reg=-1;
          return -1;
        }
        for (k = 0; k < F; k++)
        {
          cc.x[i][k] ^= cc.x[j][k];
          inv_a[i][k] ^= inv_a[j][k];
        }

        cc.x[i][i] = 1;
      }
      //  exit(1);

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < F; l++)
        {
          if (cc.x[l][i] == 1)
          {
            //#pragma omp parallel for
            for (k = 0; k < F; k++)
            {
              cc.x[l][k] ^= cc.x[i][k];
              inv_a[l][k] ^= inv_a[i][k];
            }
          }
        }

        // printf("@%d\n",i);
      }
      // printf("@i=%d\n",i);
    }

    //  exit(1);

    for (i = 1; i < F; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          for (j = 0; j < F; j++)
          {
            // if(a[k][i]==1){
            cc.x[k][j] ^= cc.x[i][j];
            inv_a[k][j] ^= inv_a[i][j];
            //}
          }
        }
      }
    }

    //検算
    for (i = 0; i < F; i++)
    {

      {
        for (j = 0; j < F; j++)
        {
          l = 0;
          for (k = 0; k < F; k++)
          {
            b[i][j] ^= (cl[i][k] & inv_a[k][j]);
            // l^=(cl[i][k]&inv_a[k][j]);
          }
          // b[i][j]=l;
        }
      }
    }

    for (i = 0; i < F; i++)
    {
      //   printf("%d",b[i][i]);
      // printf("==\n");
      if (b[i][i] == 1)
      {
        // printf("baka");
        //    exit(1);
        flg++;
      }
    }
    count = 0;

    for (i = 0; i < F; i++)
    {
      for (j = 0; j < F; j++)
      {
        if (b[i][j] == 0 && i != j)
          count++;
      }
    }
    printf("%d,%d %d,%d", flg, F, count, F * F - F);
    // wait();
    // if(cl[0][0]>0)
    //   goto labo;
    //
    printf("S[K][K]=\n{\n");
    if (flg == F && count == (F * F - F))
    // if(flg==F)
    {

      printf("inv_S[K][K]=\n{\n");
      for (i = 0; i < F; i++)
      {
        printf("{");
        for (j = 0; j < F; j++)
        {
          R->x[i][j] = inv_a[i][j];
          // O.reg=0;
          printf("%d,", inv_S.x[i][j]);
        }
         printf("},\n");
      }
       printf("};\n");
      return 0;
    }
    // O.reg= -1;
    return -1;
  }

  return -1;
}

int mkS(MTX cc, MTX *R)
{
  int i, j, k, l;
  unsigned char b[X][X] = {0};
  unsigned char dd[X] = {0};
  unsigned int flg = 0, count = 0;
  // unsigned char cc.x[X][X] = {0};
  unsigned char cl[X][X];
  time_t t;
  FILE *fq;
  unsigned char inv_a[X][X] = {0}; //ここに逆行列が入る
  unsigned char buf;               //一時的なデータを蓄える
  int n = K * E;                   //配列の次数

  // while(flg!=F || count!=F*F-F)
  // while(count!=F*F-F)
  while (flg != X)
  {
  labo:
    // memset(cc,0,sizeof(cc));
    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    // g2();
    /*
    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
        cc.x[i][j] = xor128() % 2;
    }
    */
    printf("end of g2\n");
    // exit(1);

    //#pragma omp parallel for private(j)
    for (i = 0; i < X; i++)
    {

      for (j = 0; j < X; j++)
      {
        // printf("%d,",cc.x[i][j]);
        cl[i][j] = cc.x[i][j];
        dd[j] = cc.x[i][j];
      }
      // printf("\n");
    }

    // memset(inv_a,0,sizeof(inv_a));

    //単位行列を作る
    //#pragma omp parallel for private(j)
    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    //掃き出し法

    for (i = 0; i < X; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;
        /*
  cc.x[i][i]=1;
  for(k=i+1;k<F;k++)
    cc.x[i][k]^=rand()%2;
  //printf("i=%d\n",i);
  */

        while (cc.x[j][i] == 0 && j < X)
        {
          j++;
          // buf=cc.x[j++][i];
        }

        //  cc.x[i][i]=1;
        //  printf("j=%d\n",j);

        //  exit(1);
        //#pragma omp parallel for
        if (j >= X)
        {
          printf("baka in mkS %d\n", j);
          // exit(1);
          return -1;
          // goto labo;
        }
        for (k = 0; k < X; k++)
        {
          cc.x[i][k] ^= cc.x[j][k] % 2;
          inv_a[i][k] ^= inv_a[j][k];
        }

        cc.x[i][i] = 1;
      }
      //  exit(1);

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < X; l++)
        {
          if (cc.x[l][i] == 1)
          {
            //#pragma omp parallel for
            for (k = 0; k < X; k++)
            {
              cc.x[l][k] ^= cc.x[i][k] % 2;
              inv_a[l][k] ^= inv_a[i][k];
            }
          }
        }

        // printf("@%d\n",i);
      }
      // printf("@i=%d\n",i);
    }

    //  exit(1);
    //#pragma omp parallel for private(j,k)
    for (i = 1; i < X; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          for (j = 0; j < X; j++)
          {
            // if(a[k][i]==1){
            cc.x[k][j] ^= cc.x[i][j] % 2;
            inv_a[k][j] ^= inv_a[i][j];
            //}
          }
        }
      }
    }

    /*
        //逆行列を出力
        for (i = 0; i < F; i++)
        {
          for (j = 0; j < F; j++)
          {
            printf("a %d,", inv_a[i][j]);
          }
          printf("\n");
        }
    */
    // exit(1);

    //検算
    //#pragma omp parallel for private(j, k) num_threads(16)
    for (i = 0; i < X; i++)
    {
      //#pragma omp parallel num_threads(8) //private(j,k)
      {
        for (j = 0; j < X; j++)
        {
          l = 0;
          //#pragma omp parallel for reduction(^:l)
          for (k = 0; k < X; k++)
          {
            b[i][j] ^= (cl[i][k] & inv_a[k][j]);
            // l^=(cl[i][k]&inv_a[k][j]);
          }
          // b[i][j]=l;
        }
      }
    }

    for (i = 0; i < X; i++)
    {
      //   printf("%d",b[i][i]);
      // printf("==\n");
      if (b[i][i] == 1)
      {
        // printf("baka");
        //    exit(1);
        flg++;
      }
    }
    count = 0;

    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
      {
        if (b[i][j] == 0 && i != j)
          count++;
      }
    }

    // if(cl[0][0]>0)
    //   goto labo;
    //
    printf("S[K][K]=\n{\n");
    if (flg == X && count == (X * X - X))
    // if(flg==F)
    {
      for (i = 0; i < X; i++)
      {
        // printf("{");
        for (j = 0; j < X; j++)
        {
          //
          dd[j] = cl[i][j];
          S.x[i][j] = cl[i][j];
          printf("%d,", S.x[i][j]);
        }

        printf("},\n");
      }
      printf("};\n");

      printf("inv_S[K][K]=\n{\n");
      for (i = 0; i < X; i++)
      {
        // printf("{");
        for (j = 0; j < X; j++)
        {
          dd[j] = inv_a[i][j];
          R->x[i][j] = inv_a[i][j];
          // printf("%d,", inv_S.w[i][j]);
        }
        // printf("},\n");
      }
      // printf("};\n");

      /*
            for (i = 0; i < F; i++)
            {
              for (j = 0; j < F; j++)
                printf("%d, ", b[i][j]);
              printf("\n");
            }
            //  exit(1);
            */
      return 0;
    }
    return -1;
  }
  return -1;
}

int mkRS(MTX cc, MTX *R)
{
  int i, j, k, l;
  unsigned char b[AY][AY] = {0};
  unsigned char dd[AY] = {0};
  unsigned int flg = 0, count = 0;
  // unsigned char cc.x[X][X] = {0};
  unsigned char cl[AY][AY];
  time_t t;
  FILE *fq;
  unsigned char inv_a[AY][AY] = {0}; //ここに逆行列が入る
  unsigned char buf;                 //一時的なデータを蓄える
  int n = K * E;                     //配列の次数

  // while(flg!=F || count!=F*F-F)
  // while(count!=F*F-F)
  while (flg != AY)
  {
  labo:
    // memset(cc,0,sizeof(cc));
    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    // g2();
    /*
    for (i = 0; i < X; i++)
    {
      for (j = 0; j < X; j++)
        cc.x[i][j] = xor128() % 2;
    }
    */
    printf("end of g2\n");
    // exit(1);

    //#pragma omp parallel for private(j)
    for (i = 0; i < AY; i++)
    {

      for (j = 0; j < AY; j++)
      {
        // printf("%d,",cc.x[i][j]);
        cl[i][j] = cc.x[i][j];
        dd[j] = cc.x[i][j];
      }
      // printf("\n");
    }

    // memset(inv_a,0,sizeof(inv_a));

    //単位行列を作る
    //#pragma omp parallel for private(j)
    for (i = 0; i < AY; i++)
    {
      for (j = 0; j < AY; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    //掃き出し法

    for (i = 0; i < AY; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;
        /*
  cc.x[i][i]=1;
  for(k=i+1;k<F;k++)
    cc.x[i][k]^=rand()%2;
  //printf("i=%d\n",i);
  */

        while (cc.x[j][i] == 0 && j < AY)
        {
          j++;
          // buf=cc.x[j++][i];
        }

        //  cc.x[i][i]=1;
        //  printf("j=%d\n",j);

        //  exit(1);
        //#pragma omp parallel for
        if (j >= AY)
        {
          printf("baka in mkS %d\n", j);
          // exit(1);
          return -1;
          // goto labo;
        }
        for (k = 0; k < AY; k++)
        {
          cc.x[i][k] ^= cc.x[j][k] % 2;
          inv_a[i][k] ^= inv_a[j][k];
        }

        cc.x[i][i] = 1;
      }
      //  exit(1);

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < AY; l++)
        {
          if (cc.x[l][i] == 1)
          {
            //#pragma omp parallel for
            for (k = 0; k < AY; k++)
            {
              cc.x[l][k] ^= cc.x[i][k] % 2;
              inv_a[l][k] ^= inv_a[i][k];
            }
          }
        }

        // printf("@%d\n",i);
      }
      // printf("@i=%d\n",i);
    }

    //  exit(1);
    //#pragma omp parallel for private(j,k)
    for (i = 1; i < AY; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          for (j = 0; j < AY; j++)
          {
            // if(a[k][i]==1){
            cc.x[k][j] ^= cc.x[i][j] % 2;
            inv_a[k][j] ^= inv_a[i][j];
            //}
          }
        }
      }
    }

    /*
        //逆行列を出力
        for (i = 0; i < F; i++)
        {
          for (j = 0; j < F; j++)
          {
            printf("a %d,", inv_a[i][j]);
          }
          printf("\n");
        }
    */
    // exit(1);

    //検算
    //#pragma omp parallel for private(j, k) num_threads(16)
    for (i = 0; i < AY; i++)
    {
      //#pragma omp parallel num_threads(8) //private(j,k)
      {
        for (j = 0; j < AY; j++)
        {
          l = 0;
          //#pragma omp parallel for reduction(^:l)
          for (k = 0; k < AY; k++)
          {
            b[i][j] ^= (cl[i][k] & inv_a[k][j]);
            // l^=(cl[i][k]&inv_a[k][j]);
          }
          // b[i][j]=l;
        }
      }
    }

    for (i = 0; i < AY; i++)
    {
      //   printf("%d",b[i][i]);
      // printf("==\n");
      if (b[i][i] == 1)
      {
        // printf("baka");
        //    exit(1);
        flg++;
      }
    }
    count = 0;

    for (i = 0; i < AY; i++)
    {
      for (j = 0; j < AY; j++)
      {
        if (b[i][j] == 0 && i != j)
          count++;
      }
    }

    // if(cl[0][0]>0)
    //   goto labo;
    //
    printf("S[K][K]=\n{\n");
    if (flg == AY && count == (AY * AY - AY))
    // if(flg==F)
    {
      for (i = 0; i < AY; i++)
      {
        // printf("{");
        for (j = 0; j < AY; j++)
        {
          //
          dd[j] = cl[i][j];
          S.x[i][j] = cl[i][j];
          printf("%d,", S.x[i][j]);
        }

        printf("},\n");
      }
      printf("};\n");

      printf("inv_S[K][K]=\n{\n");
      for (i = 0; i < AY; i++)
      {
        // printf("{");
        for (j = 0; j < AY; j++)
        {
          dd[j] = inv_a[i][j];
          R->x[i][j] = inv_a[i][j];
          // printf("%d,", inv_S.w[i][j]);
        }
        // printf("},\n");
      }
      // printf("};\n");

      /*
            for (i = 0; i < F; i++)
            {
              for (j = 0; j < F; j++)
                printf("%d, ", b[i][j]);
              printf("\n");
            }
            //  exit(1);
            */
      return 0;
    }
    return -1;
  }
  return -1;
}

int binv(MTX cc, MTX *L, int Y)
{
  int i, j, k, l;
  unsigned char b[N][N] = {0};
  unsigned char dd[N] = {0};
  unsigned int flg = 0, count = 0;
  // unsigned char cc.x[X][X] = {0};
  unsigned char cl[N][N];
  time_t t;
  FILE *fq;
  unsigned char inv_a[N][N] = {0}; //ここに逆行列が入る
  unsigned char buf;               //一時的なデータを蓄える
  int n = K * E;                   //配列の次数

  // while(flg!=F || count!=F*F-F)
  // while(count!=F*F-F)
  while (flg != Y)
  {
  labo:
    // memset(cc,0,sizeof(cc));
    flg = 0;
    count = 0;
    srand(clock() + time(&t));

    // g2();
    /*
    for (i = 0; i < Y; i++)
    {
      for (j = 0; j < Y; j++)
        cc.x[i][j] = xor128() % 2;
    }
    */
    printf("end of g2\n");
    // exit(1);

#pragma omp parallel for private(j)
    for (i = 0; i < Y; i++)
    {

      for (j = 0; j < Y; j++)
      {
        // printf("%d,",cc.x[i][j]);
        cl[i][j] = cc.x[i][j];
        dd[j] = cc.x[i][j];
      }
      // printf("\n");
    }

// memset(inv_a,0,sizeof(inv_a));

//単位行列を作る
#pragma omp parallel for private(j)
    for (i = 0; i < Y; i++)
    {
      for (j = 0; j < Y; j++)
      {
        inv_a[i][j] = (i == j) ? 1.0 : 0.0;
      }
    }

    //掃き出し法

    for (i = 0; i < Y; i++)
    {
      if (cc.x[i][i] == 0)
      {
        j = i;
        /*
  cc.x[i][i]=1;
  for(k=i+1;k<F;k++)
    cc.x[i][k]^=rand()%2;
  //printf("i=%d\n",i);
  */

        while (cc.x[j][i] == 0 && j < Y)
        {
          j++;
          // buf=cc.x[j++][i];
        }

        //  cc.x[i][i]=1;
        //  printf("j=%d\n",j);

        //  exit(1);
        //#pragma omp parallel for
        if (j >= Y)
        {
          printf("baka in binv %d\n", j);
          // exit(1);
          return -1;
          // goto labo;
        }
        for (k = 0; k < Y; k++)
        {
          cc.x[i][k] ^= cc.x[j][k] % 2;
          inv_a[i][k] ^= inv_a[j][k];
        }

        cc.x[i][i] = 1;
      }
      //  exit(1);

      if (cc.x[i][i] == 1)
      {
        for (l = i + 1; l < Y; l++)
        {
          if (cc.x[l][i] == 1)
          {
            //#pragma omp parallel for
            for (k = 0; k < Y; k++)
            {
              cc.x[l][k] ^= cc.x[i][k] % 2;
              inv_a[l][k] ^= inv_a[i][k];
            }
          }
        }

        // printf("@%d\n",i);
      }
      // printf("@i=%d\n",i);
    }

    //  exit(1);
    //#pragma omp parallel for private(j,k)
    for (i = 1; i < Y; i++)
    {
      for (k = 0; k < i; k++)
      {
        if (cc.x[k][i] == 1)
        {
          for (j = 0; j < Y; j++)
          {
            // if(a[k][i]==1){
            cc.x[k][j] ^= cc.x[i][j] % 2;
            inv_a[k][j] ^= inv_a[i][j];
            //}
          }
        }
      }
    }

    /*
        //逆行列を出力
        for (i = 0; i < F; i++)
        {
          for (j = 0; j < F; j++)
          {
            printf("a %d,", inv_a[i][j]);
          }
          printf("\n");
        }
    */
    // exit(1);

    //検算
#pragma omp parallel for private(j, k) num_threads(16)
    for (i = 0; i < Y; i++)
    {
      //#pragma omp parallel num_threads(8) //private(j,k)
      {
        for (j = 0; j < Y; j++)
        {
          l = 0;
          //#pragma omp parallel for reduction(^:l)
          for (k = 0; k < Y; k++)
          {
            b[i][j] ^= (cl[i][k] & inv_a[k][j]);
            // l^=(cl[i][k]&inv_a[k][j]);
          }
          // b[i][j]=l;
        }
      }
    }

    for (i = 0; i < Y; i++)
    {
      //   printf("%d",b[i][i]);
      // printf("==\n");
      if (b[i][i] == 1)
      {
        // printf("baka");
        //    exit(1);
        flg++;
      }
    }
    count = 0;

    for (i = 0; i < Y; i++)
    {
      for (j = 0; j < Y; j++)
      {
        if (b[i][j] == 0 && i != j)
          count++;
      }
    }

    // if(cl[0][0]>0)
    //   goto labo;
    //
    printf("S[K][K]=\n{\n");
    if (flg == Y && count == (Y * Y - Y))
    // if(flg==F)
    {
      for (i = 0; i < Y; i++)
      {
        // printf("{");
        for (j = 0; j < Y; j++)
        {
          //
          dd[j] = cl[i][j];
          S.x[i][j] = cl[i][j];
          printf("%d,", S.x[i][j]);
        }

        printf("},\n");
      }
      printf("};\n");

      printf("inv_S[K][K]=\n{\n");
      for (i = 0; i < Y; i++)
      {
        // printf("{");
        for (j = 0; j < Y; j++)
        {
          dd[j] = inv_a[i][j];
          L->x[i][j] = inv_a[i][j];
          // printf("%d,", inv_S.w[i][j]);
        }
        // printf("},\n");
      }
      // printf("};\n");

      /*
            for (i = 0; i < F; i++)
            {
              for (j = 0; j < F; j++)
                printf("%d, ", b[i][j]);
              printf("\n");
            }
            //  exit(1);
            */
      return 0;
    }
    return -1;
  }

  /*
  fq=fopen("S.key","wb");
  for(i=0;i<F;i++){
    for(j=0;j<F;j++)
      dd[j]=cl[i][j];
    fwrite(dd,1,n,fq);

  }
  fclose(fq);
  fq=fopen("inv_S.key","wb");
  for(i=0;i<F;i++){
    for(j=0;j<F;j++)
      dd[j]=inv_a[i][j];
    fwrite(dd,1,n,fq);
  }
  fclose(fq);
*/

  // free(b);
  return -1;
}
