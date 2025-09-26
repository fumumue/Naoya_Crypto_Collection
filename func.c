
// 20200816:正規化したいところだがうまく行かない
// 多項式の足し算
vec vsub(vec a, vec b)
{
    vec c = {0};
    // int i, j, k, l = 0;
    vec h = {0}, f2 = {0}, g2 = {0};

    for (int i = 0; i < DEG; i++)
    {
        if (a.x[i] >= b.x[i])
            c.x[i] = (a.x[i] - b.x[i]) % N;
        if (a.x[i] < b.x[i])
            c.x[i] = (N + a.x[i] - b.x[i]) % N;
    }

    return c;
}


int mul = 0, mul2 = 0;
vec vmul(vec a, vec b,int R)
{
    int i, j, k, l;
    vec c = {0};

    k = deg(a);
    l = deg(b);

    if(l+k>N){
        printf("blake %d\n",l+k);        
        exit(1);
    }
    i = 0;
    while (i < k + 1)
    {
        for (j = 0; j < l + 1; j++)
        {
            if (a.x[i] > 0)
                c.x[i + j] = (c.x[i + j] + a.x[i] * b.x[j]) % R;
        }
        i++;
    }

    return c;
}

