# USAGE
ulimit -s unlimited

gcc -O3 naoya.c

./a.out

# 小型軽量公開鍵暗号

# 1. 設計思想

1. **秘密点の利用**
   多項式の解 $$(P = \\{ a_1, ..., a_n \\}, a_i \in F_q),l=(x-a_0)(x-a_1)...(x-a_n)$$ を秘密鍵とする。
   ここで、$q$は256ビット素数であるとする。
2. **平文多項式の消去**
   暗号化では、平文多項式を$m,deg(m)=s$であり、$lm$を秘密点$a_i$で評価した際に消えるようにしたい。

---

# 2. 暗号化手順（1変数版）
以下、qを256ビット素数とし$Z_q$上の多項式を考える。$a_i \in Z_q$とする。
1. 秘密多項式:　$$l=(x-a_0)(x-a_1)...(x-a_{d+s}),n=d+s+1$$とする。

  ここで$a_i$は、重複なくランダムに選んだ値$a \in F_q$とし、これが秘密鍵である。
# 概要

本稿では、多項式の０点を用いた新しい公開鍵暗号方式「消える平文暗号（Vanishing Plaintext Cipher）」を提案する。この方式は、平文を暗号文の中で一時的に「消す」ことで、秘密情報を保護しつつ効率的な復号を可能にする。従来のNTRUやMcEliece系の符号暗号とは異なり、公開鍵に秘密多項式を組み込むことで、解析耐性を高める構造となっている。

---

# 1. 設計思想

1. **秘密点の利用**
   多項式の解 $$(P = \\{ a_1, ..., a_n \\}, a_i \in F_q),l=(x-a_0)(x-a_1)...(x-a_n)$$ を秘密鍵とする。
   ここで、$q$は256ビット素数であるとする。
2. **平文多項式の消去**
   暗号化では、平文多項式を$m,deg(m)=s$であり、$lm$を秘密点$a_i$で評価した際に消えるようにしたい。

---

# 2. 暗号化手順（1変数版）
以下、qを256ビット素数とし$Z_q$上の多項式を考える。$a_i \in Z_q$とする。
1. 秘密多項式:$$l=(x-a_0)(x-a_1)...(x-a_{d+s}),n=d+s+1$$とする。

  ここで$a_i$は、重複なくランダムに選んだ値$a \in F_q$とし、これが秘密鍵である。

2. ランダム多項式 $R,deg(R)=d$ を生成する。

3. 公開鍵生成:$C_a$を秘密多項式$l$にかける秘密の値とする。公開鍵は、
   $$
   \text{pub}(x) = R(x) + l'(x), \quad l'=C_al
   $$
  であるとする。（これによって秘密多項式ｌがスケーリングされる）
  
4. 暗号化:
   平文多項式 $m$ とランダム多項式 $r,R'$ を用いて
   $$
   c = pub(m)+R' = (R + l')(m+r) + R'
   $$
   とし、$c_0=(c,r)$を暗号文として送信する。
  ここで、$deg(R')=d+s,deg(r)=s+1,deg(m)=s$であるとする。

---

# 3. 復号手順

### 1. 秘密点 $(a_i \in F_q)$ を用いて、暗号文から平文成分を消去。
  $$c=(R+l')(m+r)+R'=Rm+Rr+l'm+l'r+R'$$

  より、$deg(l')=d+s+1$であり、$d+s+1$個の点で消える。
  つまり$l'(a)=0$から$(l'm+l'r)(a)=0$であり、$c_s=Rm+Rr+R'$とすると、 
  $deg(c_s)=d+s$なので、$a$の$c_s$への代入値よりラグランジュ補間法によって$c_s$が計算で求まる。
  $$c(a)=(Rm+Rr+R')(a)$$
  ここで、復号は秘密点、$$a_i \in F_q,(0 \leq i \leq d+s)$$を知っている受信者にしかできない。

 
### 2. 消去後、残りの多項式を l'(x) で割ることでmを再構成。
$$c-c'=l'm+l'r$$
さらに、受信者は$l',r$を知っていることから、$$m'=(c-c')/l'=m+r$$

最後に$m'-r=m$として平文が復号できた。

---

# 4. 安全性の特徴

* 公開鍵だけでは R と l' を分離できない。（加算により潰れる）
* 秘密点を使わないと平文を消した後の多項式はランダム多項式に見えるため、攻撃者には情報が得られない。
* もし公開鍵$pk=R+l'$から$l'$の上位$s+1$次の項がわかっても、全ての評価点$a_i$は未提携数行が多いために決定できない。（自由度が大きすぎる）

* 攻撃者は公開鍵からRを消去してｌを因数分解して評価点を求めようとすることができる。
したがって、Rの総探索やバースデーアタックによって秘密鍵を探すことができる。
おすすめのパラメーターは、$d=42,s=41$となるであろう。

---

# 5. 実装の工夫

* 1変数版でも十分安全性が確保可能。
* 素体$F_q$もしくは$Z/oZ$を大きく取る（256ビットなど）ことで解析コストを増大。
* ラグランジュ補間や評価多項式を用いることで、復号操作を効率化。
* 小規模実装が可能なため、作りやすい。
* 有限体上の多項式の四則演算ができれば実現可能。（加算、減算、商、積）

---

# 6. 考察

* NTRUのような格子ベース暗号に依存せず、平文消去による独自の安全性を実現。
* 公開鍵に秘密多項式を組み込むという手法は、従来の公開鍵方式とは異なる新しいアプローチ。
* 今後、2変数版やMQ版などへの拡張も検討可能。

---

# まとめ

本方式「消える平文暗号」は、秘密点と多項式構造を活用することで、平文を暗号文の中で消去しながら復号を可能にする新しい公開鍵暗号である。先頭係数の秘密化により公開鍵から秘密情報を抽出できず、解析耐性が高い構造となっていると思われる。

---
# 感想
これは以前もやっていたんですが、できないことがわかってすぐ消しました。
なぜ以前はできなかったのかというと、余計なランダム式を消そうとして、ランダム多項式側に$l$をかけていたからです。
今回は公開鍵に秘密鍵を内蔵する形にして、復号の時、先に暗号化した平文を秘密点で消すようにトラップドアを工夫したところうまく行きました。

今までもいろいろ暗号のアイデアを出してきたんですが、今回ほど出来のいいのは見つかりませんでした。
自分は今まで、強くて小型で安全性が証明できる暗号の開発を続けていたのですが、ついにその候補が出てきた感じです。

そして、HQCの仕様を見て驚きました。
符号に過剰なエラーを入れた時どうなるのか、それをうまく消すことができるのか？
実はこの取り組みは解決法が見つからずに断念したのですが、HQCがまさにその方法だったのです。

最近はNTRUに見られるように実装の簡単な公開鍵暗号を作ろうとしていたのですが、符号ともMQ暗号とも違う平文を消すという方法を思いついて、それをうまく実現できる方法がひらめいたのです。
それが今回のこの暗号です。
まさにHQCの発送の双対と言った感じです。
**今の所解読できてません。興味のある方はぜひ解読してみてください。**

評価点を多項式表現したものがｌなわけですが、実際ランダムな係数を持つRによって隠れるのは低次数のところだけということで、半分しかわかっていない分離多項式の根はどうやって求めるのかが、この暗号の肝になそうです。

一般的に全ての係数がわからないと評価のしようがないので、式の下半分がわからないときに評価点を求めるのは難しいらしいです。
なので最初高次の部分だけが見えている状態が気持ち悪かったのですが、公開鍵を$$pk=lr+R$$というように、新たに低次のランダム多項式$r$をかけることで係数の露出はなくせるようです。


# 付録：ChatGPTによるラグランジュ補間プログラム

**ラグランジュ補間法のサンプルプロフラム**

```c
#include <stdio.h>
#include <stdint.h>

#define N 1129

#define DEG 1024+1
typedef struct {
    int x[DEG];
} vec;

// ---- GF(17) 演算 ----
#define Q 17  // 素体の法

int gf_add(int a, int b) {
    int r = (a + b) % Q;
    if (r < 0) r += Q;
    return r;
}

int gf_mul(int a, int b) {
    long long r = (long long)a * b % Q;
    if (r < 0) r += Q;
    return (int)r;
}

int gf_pow(int a, int e) {
    long long res = 1, base = a % Q;
    while (e > 0) {
        if (e & 1) res = (res * base) % Q;
        base = (base * base) % Q;
        e >>= 1;
    }
    return (int)res;
}

// 逆元 (a^(Q-2) mod Q) — フェルマーの小定理
int gf_inv(int a) {
    if (a == 0) {
        fprintf(stderr, "division by zero in GF(17)!\n");
        return 0;
    }
    return gf_pow(a, Q - 2);
}

int gf_div(int a, int b) {
    return gf_mul(a, gf_inv(b));
}

// ---- ラグランジュ補間 ----
void lagrange_interpolate(const int *xs, const int *ys, int k, int *result) {
    for (int i = 0; i < k; i++) result[i] = 0;

    for (int i = 0; i < k; i++) {
        // 基底多項式 l_i(x) の分子
        int numer[64] = {0}; // k が小さい前提
        numer[0] = 1;
        int numer_deg = 0;
        int denom = 1;

        for (int j = 0; j < k; j++) {
            if (j == i) continue;

            int new_numer[64] = {0};
            for (int d = 0; d <= numer_deg; d++) {
                int c = numer[d];
                if (c == 0) continue;

                // (x - xj) = x + (-xj) だが mod Q
                new_numer[d]   = gf_add(new_numer[d], gf_mul(c, (Q - xs[j]) % Q));
                new_numer[d+1] = gf_add(new_numer[d+1], c);
            }
            for (int d = 0; d <= numer_deg+1; d++) numer[d] = new_numer[d];
            numer_deg++;

            denom = gf_mul(denom, gf_add(xs[i], (Q - xs[j]) % Q)); // xi - xj
        }

        int scale = gf_div(ys[i], denom);

        for (int d = 0; d <= numer_deg; d++) {
            if (numer[d] == 0) continue;
            result[d] = gf_add(result[d], gf_mul(numer[d], scale));
        }
    }
}


int gf_sub(int a, int b) { return (a - b + Q) % Q; }
// 多項式の表示
void poly_print(const uint8_t *p, int deg) {
    for (int i = deg; i >= 0; i--) {
        if (p[i]) printf("1");
        else printf("0");
    }
    printf("\n");
}



// 逆元（フェルマーの小定理 a^(Q-2) mod Q）
int gf_inv_q(int a) {
    int res = 1, base = a % Q, e = Q - 2;
    while (e > 0) {
        if (e & 1) res = gf_mul(res, base);
        base = gf_mul(base, base);
        e >>= 1;
    }
    return res;
}

// 多項式の掛け算 mod Q
void poly_mul_q(const int *a, int deg_a, const int *b, int deg_b, int *res) {
    for (int i = 0; i <= deg_a + deg_b; i++) res[i] = 0;
    for (int i = 0; i <= deg_a; i++)
        for (int j = 0; j <= deg_b; j++)
            res[i + j] = gf_add(res[i + j], gf_mul(a[i], b[j]));
}

// 多項式の割り算（商）mod Q
void poly_div_q(const int *num, int deg_num, const int *den, int deg_den, int *quot) {
    int rem[DEG*2+1];
    for (int i = 0; i <= deg_num; i++) rem[i] = num[i];
    for (int i = 0; i <= deg_num - deg_den; i++) quot[i] = 0;

    for (int i = deg_num; i >= deg_den; i--) {
        if (rem[i] != 0) {
            int factor = gf_div(rem[i], den[deg_den]);
            quot[i - deg_den] = factor;
            for (int j = 0; j <= deg_den; j++)
                rem[i - deg_den + j] = gf_sub(rem[i - deg_den + j], gf_mul(factor, den[j]));
        }
    }
}



int mltn(int n,int x){
    int t=1;
    for(int i=0;i<n;i++){
    t*=x;
    t%=N;
    }
    return t;
}


// 多項式の代入値
long long int
trace(vec f,int deg, long long int x)
{
    long long int u = 0;
    vec v = (f);

    for (int i = 0; i < deg+1; i++)
    {
        if (v.x[i] > 0){
            u = (u + (v.x[i] * mltn(i, x))) % N;
            //printf("u%d %d\n",u,i);
        }
    }

    return u;
}




// ---- main ----
int main(void) {
    // 元の多項式 f(x) = 3 + 5x + 2x^2
    vec f = {3, 5, 2};
    int deg = 2;

    // 評価点
    int xs[3] = {1, 2, 4};
    int ys[3];

    for (int i = 0; i < 3; i++) {
        ys[i] = trace(f,deg, xs[i]);
        printf("f(%d) = %d\n", xs[i], ys[i]);
    }

    // 補間
    int rec[3];
    lagrange_interpolate(xs, ys, 3, rec);

    printf("\nrecovered polynomial:\n");
    for (int i = 0; i <= deg; i++) {
        printf("x^%d: %d\n", i, rec[i]);
    }

    return 0;
}
```

2. ランダム多項式 $R,deg(R)=d$ を生成する。

3. 公開鍵生成:$C_a$を秘密多項式$l$にかける秘密の値とする。公開鍵は、
   $$
   \text{pub}(x) = R(x) + l'(x), \quad l'=C_al
   $$
  であるとする。（これによって秘密多項式ｌがスケーリングされる）
  
4. 暗号化:
   平文多項式 $m$ とランダム多項式 $r,R'$ を用いて
   $$
   c = pub(m)+R' = (R + l')(m+r) + R'
   $$
   とし、$c_0=(c,r)$を暗号文として送信する。
  ここで、$deg(R')=d+s,deg(r)=s+1,deg(m)=s$であるとする。

---

# 3. 復号手順

### 1. 秘密点 $(a_i \in F_q)$ を用いて、暗号文から平文成分を消去。
  $$c=(R+l')(m+r)+R'=Rm+Rr+l'm+l'r+R'$$

  より、$deg(l')=d+s+1$であり、$d+s+1$個の点で消える。
  つまり$l'(a)=0$から$(l'm+l'r)(a)=0$であり、$c_s=Rm+Rr+R'$とすると、 
  $deg(c_s)=d+s$なので、$a$の$c_s$への代入値よりラグランジュ補間法によって$c_s$が計算で求まる。
  $$c(a)=(Rm+Rr+R')(a)$$
  ここで、復号は秘密点、$$a_i \in F_q,(0 \leq i \leq d+s)$$を知っている受信者にしかできない。

 
### 2. 消去後、残りの多項式を l'(x) で割ることでmを再構成。
$$c-c'=l'm+l'r$$
さらに、受信者は$l',r$を知っていることから、$$m'=(c-c')/l'=m+r$$

最後に$m'-r=m$として平文が復号できた。

