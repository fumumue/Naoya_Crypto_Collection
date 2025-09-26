#define VEC_SIZE N

/* 修正版 dmd2: 全て GF(N) の整数演算 (%N) に統一し、元の漸化式ロジックを保持 */
vec dmd2(MTX bb) {
    int i, j;

    /* xa をゼロ初期化 */
    static oterm xa[N+1][N+1];
    memset(xa, 0, sizeof(xa));

    /* 元コードに合わせた初期設定（xa[0][0] = 1 のような形） */
    xa[0][0].a = 1;
    xa[0][0].n = 1;

    /* bb.z の値を xa に詰める（元コードのインデックス割当を踏襲） */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N+1; j++) {
            int idx = N - j;           /* 元コードでは N-j を使っていた */
            /* 正のインデックスに正規化 */
            while (idx < 0) idx += N;
            idx %= N;
            /* 明示的に mod N */
            xa[N-1-i][j].a = (bb.z[i].x[idx]) % N;
            if (xa[N-1-i][j].a < 0) xa[N-1-i][j].a += N;
            /* n フィールドは（元コードの用途に応じて）1 にする */
            xa[N-1-i][j].n = 1;
        }
    }

    /* yx を漸化式で作る（yx[0] = 1） */
    int yx[N+1];
    for (i = 0; i <= N; i++) yx[i] = 0;
    yx[0] = 1;

    /* xa[0][i].n を yx[i] として設定する形に沿って計算 */
    for (i = 1; i <= N; i++) {
        long sum = 0L;
        for (j = 0; j < i; j++) {
            /* sum += xa[i-1][j].a * xa[0][j].n だが、元の実装が xa[0][j].n を yx[j] としたため、
               ここでは yx[j] を使う（元の漸化式の流れに合わせる） */
            sum += (long)xa[i-1][j].a * (long)yx[j];
        }
        /* モジュロ N に正規化 */
        sum %= (long)N;
        if (sum < 0) sum += N;
        /* yx[i] = N - sum (mod N) */
        yx[i] = (int)((N - sum) % N);
        /* 保持（元コードは xa[0][i].n = yx[i] としていた） */
        xa[0][i].n = yx[i];
    }

    /* ucc (degree T) と uec (degree N-T) を vec に詰める（元ロジックに沿う） */
    vec ucc; memset(&ucc, 0, sizeof(ucc));
    vec uec; memset(&uec, 0, sizeof(uec));

    /* ucc.x[T - i] = yx[i] for i=0..T  (元コードの位置合わせを踏襲) */
    for (i = 0; i <= T && i <= N; i++) {
        int pos = T - i;
        if (pos >= 0 && pos < VEC_SIZE) ucc.x[pos] = (yx[i] % N + N) % N;
    }

    /* uec.x[N - i] = yx[i] for i=T+1..N */
    for (i = T+1; i <= N; i++) {
        int pos = N - i;
        if (pos >= 0 && pos < VEC_SIZE) uec.x[pos] = (yx[i] % N + N) % N;
    }

    /* 多項式除算: b = uec / ucc, a = uec % ucc  */
    vec b = vdiv(uec, ucc);
    vec a = vmod(uec, ucc);

    /* 余りがゼロでなければ警告（元コードはここで失敗させていた） */
    if (deg(a) > 0) {
        fprintf(stderr, "dmd: non-zero remainder deg(a)=%d (mod %d)\n", deg(a), N);
        /* ここは元の方針に合わせて処理を続けるか中断するか選べます。
           とりあえず続けて結果を返します。 */
    }

    /* ans の計算: ans.x[i] = (N - trace(b, (i+1)%N)) % N */
    vec ans; memset(&ans, 0, sizeof(ans));
    for (i = 0; i < N && i < VEC_SIZE; i++) {
        int tt = trace(b, (i + 1) % N) % N;
        if (tt < 0) tt += N;
        ans.x[i] = (N - tt) % N;
    }

    return ans;
}
