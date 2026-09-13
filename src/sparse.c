/*
sparse.c

疎行列の LU 分解 (前処理の Schur 補元用)。

precond.c の Schur 補元 S = S11 - sum_b B_b D_b^-1 C_b は「回路未知数どうしの
行列」で、葉ごとのクリークを足し合わせた**局所的な**パターンになる。実測では
1 行あたりの非ゼロ数が規模によらず ~36 に飽和する (12x12 の板で 27.6、
20x20 で 35.4、28x28 で 36.3、40x40 で 35.8) ので、nnz は O(N) しかない。
密に持つと n1^2 = O(N^2) で、40x40 の板では前処理 45 MB のうち 43 MB が
これになる (全体 128 MB の 1/3)。**圧縮経路に残った最後の O(N^2)** なので、
ここを疎にする。

方針は「前処理そのものは一切変えず、格納と分解だけ疎にする」。分解は
近似ではなく厳密な LU なので、前処理の性質 (= GMRES の反復数) は変わらない。

構成 :
  1. 三つ組 (COO) から CSC を組む (重複は加算)
  2. 最小次数による充填低減順序付け (mindeg_order)
  3. Gilbert-Peierls の左向き LU + 閾値部分ピボット (slu_build)
  4. 疎な前進・後退代入 (slu_solve)

順序付けの選択は実測による (S のパターンを取り出して充填を数えた) :

  | n1   | nnz   | 密        | 自然順      | RCM        | 最小次数    |
  |------|-------|-----------|-------------|------------|-------------|
  |  441 | 15591 |   194481  | 52655 (3.4) | 43037(2.8) | 20175 (1.3) |
  |  841 | 30531 |   707281  |139007 (4.6) | 88817(2.9) | 43175 (1.4) |
  | 1681 | 60193 |  2825761  |384035 (6.4) |293493(4.9) | 97039 (1.6) |

括弧内は nnz 比。最小次数は充填が nnz の 1.3〜1.6 倍にとどまり、密比では
9.6x -> 16x -> 29x と**規模が大きいほど有利**になる。RCM は 2.8〜4.9 倍で
伸びも速いので採用しない。

最小次数は消去グラフを陽に持つ素朴な実装。消去グラフの大きさは充填そのもの
なので、上の実測どおり nnz の数倍で収まる。それでも病的な入力はありうるので
MD_MAXFILL で打ち切り、超えたら自然順に落とす (順序が悪くなるだけで結果は
変わらない)。

スレッド数不変性 : 分解も求解も直列。GMRES から前処理として呼ばれるだけで、
呼び出し順は固定なのでスレッド数を変えても結果はビット単位で一致する。
*/

#include "peec.h"

#define SLU_TOL     0.1               // 閾値ピボット (対角を |x| >= tol*max で優先)
#define MD_MAXFILL  64                // 消去グラフの上限 (nnz の倍数)

// 疎な LU : L*U = A(p, q)。L は単位下三角
struct sparse_lu_t {
	int n;
	int *lp, *li;  d_complex_t *lx;    // L (CSC、対角は含まない)
	int *up, *ui;  d_complex_t *ux;    // U (CSC、対角を最後に持つ)
	int *pinv;                         // 行 i -> ピボット位置
	int *q;                            // 列順序 (q[k] = k 番目に消去する列)
	d_complex_t *w;                    // [n] 作業領域
	int *xi;                           // [2n] DFS スタック / パターン
};

static double mag(d_complex_t z)
{
	return fabs(z.r) + fabs(z.i);
}

/* ───────────────────────── 最小次数順序付け ─────────────────────────
消去グラフを陽に持つ素朴な最小次数。A + A^T のパターンを使う (LU の充填は
非対称だが、順序付けの目安としては対称化したもので十分)。

次数バケット (deg -> 頂点の双方向リスト) で最小次数の頂点を O(1) で取り出し、
消去のたびに隣接集合をマージする。隣接は頂点ごとの可変長配列。

戻り値 0 = 正常 (order[] に消去順)、1 = 打ち切り (呼び出し側は自然順に落とす)
*/
typedef struct {
	int *a;                            // 隣接 (容量 cap、長さ len)
	int len, cap;
} mdadj_t;

static int md_push(mdadj_t *v, int w, size_t *total, size_t limit)
{
	if (v->len >= v->cap) {
		const int nc = (v->cap > 0) ? (2 * v->cap) : 8;
		int *na = (int *)realloc(v->a, (size_t)nc * sizeof(int));
		if (na == NULL) return 1;
		// realloc 済みのポインタは先に引き取る (打ち切り時に解放できるように)
		*total += (size_t)(nc - v->cap);
		v->a = na;
		v->cap = nc;
		if (*total > limit) return 1;
	}
	v->a[v->len++] = w;

	return 0;
}

static int mindeg_order(int n, const int *ap, const int *ai, int *order)
{
	// A + A^T のパターン (自己ループは除く)
	mdadj_t *adj = (mdadj_t *)calloc((size_t)(n > 0 ? n : 1), sizeof(mdadj_t));
	int *mark = (int *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
	int *deg  = (int *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
	// 次数バケット (双方向リスト)
	int *head = (int *)malloc((size_t)(n + 1) * sizeof(int));
	int *next = (int *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
	int *prev = (int *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
	int *nbr  = (int *)malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
	int fail = ((adj == NULL) || (mark == NULL) || (deg == NULL) || (head == NULL)
		|| (next == NULL) || (prev == NULL) || (nbr == NULL));

	const size_t limit = ((size_t)MD_MAXFILL * (size_t)ap[n]) + 1024;
	size_t total = 0;

	if (!fail) {
		for (int i = 0; i < n; i++) mark[i] = -1;
		for (int j = 0; (j < n) && !fail; j++) {
			for (int e = ap[j]; e < ap[j + 1]; e++) {
				const int i = ai[e];
				if (i == j) continue;
				// 重複を避けるため、片側の走査で両方向に入れる
				if (mark[i] != j) {
					mark[i] = j;
					if (md_push(&adj[j], i, &total, limit)) { fail = 1; break; }
					if (md_push(&adj[i], j, &total, limit)) { fail = 1; break; }
				}
			}
		}
	}

	// **重複の除去は必須** : 上の走査は列ごとにしか重複を見ていないので、
	// (i,j) と (j,i) が両方格納されていると同じ辺を 2 回積んでしまう。
	// すると次数 adj[i].len が真の次数の 2 倍まで膨らみ、次数バケットの
	// 添字 head[deg] (head の長さは n+1) を突き抜けてヒープを壊す。
	// 実際に踏んだ : n=18 の Schur で adj[0].len=34 -> head[34] を書き、
	// 隣の確保領域 (MNA の三つ組) が壊れて 3 周波数目で結果が狂った。
	// ASan でも再現しなかった (確保のレイアウト次第) ので、ここで必ず潰す。
	// 除去後は「自分以外の相異なる頂点」なので len <= n-1 が保証される。
	int tag = 0;                       // mark[] の世代 (頂点番号は使い回せない)
	if (!fail) {
		for (int i = 0; i < n; i++) mark[i] = -1;
		for (int v = 0; v < n; v++) {
			const int t = ++tag;
			int len = 0;
			for (int e = 0; e < adj[v].len; e++) {
				const int x = adj[v].a[e];
				if (mark[x] != t) {
					mark[x] = t;
					adj[v].a[len++] = x;
				}
			}
			adj[v].len = len;
		}
	}

	if (!fail) {
		for (int i = 0; i <= n; i++) head[i] = -1;
		for (int i = 0; i < n; i++) mark[i] = -1;
		for (int i = 0; i < n; i++) {
			deg[i] = adj[i].len;
			next[i] = head[deg[i]];
			prev[i] = -1;
			if (head[deg[i]] >= 0) prev[head[deg[i]]] = i;
			head[deg[i]] = i;
		}
	}

	int mind = 0;
	for (int k = 0; (k < n) && !fail; k++) {
		while ((mind <= n) && (head[mind] < 0)) mind++;
		if (mind > n) { fail = 1; break; }
		const int v = head[mind];

		// v をバケットから外す
		head[mind] = next[v];
		if (next[v] >= 0) prev[next[v]] = -1;
		deg[v] = -1;                   // 消去済みの印
		order[k] = v;

		// 生きている隣接を集める (消去済みを外し、重複も落とす。
		// これで m <= 生存頂点数 <= n が保証され、nbr[] があふれない)
		int m = 0;
		const int tv = ++tag;
		for (int e = 0; e < adj[v].len; e++) {
			const int w = adj[v].a[e];
			if ((deg[w] >= 0) && (mark[w] != tv)) {
				mark[w] = tv;
				nbr[m++] = w;
			}
		}
		free(adj[v].a);
		adj[v].a = NULL;
		adj[v].len = 0;
		adj[v].cap = 0;

		// クリーク化 : 各 w に、まだ隣接していない他の隣接を足す
		for (int b = 0; (b < m) && !fail; b++) {
			const int w = nbr[b];
			const int t = ++tag;
			// w の隣接を圧縮しつつ印を付ける (消去済みと重複を落とす)
			int len = 0;
			for (int e = 0; e < adj[w].len; e++) {
				const int x = adj[w].a[e];
				if ((deg[x] >= 0) && (mark[x] != t)) {
					mark[x] = t;
					adj[w].a[len++] = x;
				}
			}
			adj[w].len = len;
			mark[w] = t;               // 自分自身は入れない
			for (int c = 0; c < m; c++) {
				const int x = nbr[c];
				if (mark[x] != t) {
					mark[x] = t;
					if (md_push(&adj[w], x, &total, limit)) { fail = 1; break; }
				}
			}
			if (fail) break;

			// 次数を更新してバケットを張り替える
			const int nd = adj[w].len;
			if (nd != deg[w]) {
				if (prev[w] >= 0) next[prev[w]] = next[w];
				else head[deg[w]] = next[w];
				if (next[w] >= 0) prev[next[w]] = prev[w];
				deg[w] = nd;
				next[w] = head[nd];
				prev[w] = -1;
				if (head[nd] >= 0) prev[head[nd]] = w;
				head[nd] = w;
				if (nd < mind) mind = nd;
			}
		}
	}

	for (int i = 0; i < n; i++) free(adj[i].a);
	free(adj); free(mark); free(deg);
	free(head); free(next); free(prev); free(nbr);

	return fail;
}

/* ───────────────────────── Gilbert-Peierls LU ───────────────────────── */

// L の到達可能集合を深さ優先で求める (再帰なしの明示スタック)。
// xi[top..n-1] に、消去順で降順のパターンが入る。戻り値 = top
static int reach(int n, const int *lp, const int *li, const int *pinv,
	const int *ap, const int *ai, int j, int *xi, int *mark, int k)
{
	int top = n;

	// xi[0..head] を DFS スタック、xi[top..n-1] を出力に使う。両者は
	// 「まだ辿っている頂点」と「もう辿り終えた頂点」で互いに素なので、
	// 合わせて n を超えず、同じ配列を両端から使っても衝突しない。
	// 再開位置 (pstack) は xi[n..2n-1] に持つ。
	for (int e = ap[j]; e < ap[j + 1]; e++) {
		const int i = ai[e];
		if (mark[i] == k) continue;

		int head = 0;
		xi[0] = i;
		mark[i] = k;
		xi[n + 0] = (pinv[i] >= 0) ? lp[pinv[i]] : 0;
		while (head >= 0) {
			const int p = xi[head];
			const int pj = pinv[p];
			int done = 1;
			if (pj >= 0) {
				for (int e2 = xi[n + head]; e2 < lp[pj + 1]; e2++) {
					const int i2 = li[e2];
					if (mark[i2] == k) continue;
					mark[i2] = k;
					xi[n + head] = e2 + 1;   // 親の再開位置
					head++;
					xi[head] = i2;
					xi[n + head] = (pinv[i2] >= 0) ? lp[pinv[i2]] : 0;
					done = 0;
					break;
				}
			}
			if (done) {
				// p は辿り終えた : 逆行順 = 位相順で出力側に積む
				head--;
				top--;
				xi[top] = p;
			}
		}
	}

	return top;
}

// x = L \ A(:,j) を疎に解く。パターンは xi[top..n-1]
static int spsolve(int n, const int *lp, const int *li, const d_complex_t *lx,
	const int *pinv, const int *ap, const int *ai, const d_complex_t *ax,
	int j, int *xi, int *mark, int k, d_complex_t *x)
{
	const int top = reach(n, lp, li, pinv, ap, ai, j, xi, mark, k);

	for (int p = top; p < n; p++) x[xi[p]] = d_complex(0, 0);
	for (int e = ap[j]; e < ap[j + 1]; e++) x[ai[e]] = d_add(x[ai[e]], ax[e]);

	for (int p = top; p < n; p++) {
		const int i = xi[p];
		const int pj = pinv[i];
		if (pj < 0) continue;          // まだピボットになっていない行
		for (int e = lp[pj]; e < lp[pj + 1]; e++) {
			x[li[e]] = d_sub(x[li[e]], d_mul(lx[e], x[i]));
		}
	}

	return top;
}

static int grow(int **idx, d_complex_t **val, int need, int *cap)
{
	if (need <= *cap) return 0;
	int nc = (*cap > 0) ? (2 * (*cap)) : 1024;
	while (nc < need) nc *= 2;
	int *ni = (int *)realloc(*idx, (size_t)nc * sizeof(int));
	if (ni == NULL) return 1;
	*idx = ni;
	d_complex_t *nv = (d_complex_t *)realloc(*val, (size_t)nc * sizeof(d_complex_t));
	if (nv == NULL) return 1;
	*val = nv;
	*cap = nc;

	return 0;
}

void slu_free(struct sparse_lu_t *s)
{
	if (s == NULL) return;
	free(s->lp); free(s->li); free(s->lx);
	free(s->up); free(s->ui); free(s->ux);
	free(s->pinv); free(s->q); free(s->w); free(s->xi);
	free(s);
}

// 三つ組 (ri, ci, val) から A を組み、LU 分解する。重複は加算する。
// 戻り値 NULL = 失敗 (特異またはメモリ不足)
struct sparse_lu_t *slu_build(int n, int nnz, const int *ri, const int *ci,
	const d_complex_t *val)
{
	struct sparse_lu_t *s = (struct sparse_lu_t *)calloc(1, sizeof(struct sparse_lu_t));
	if (s == NULL) return NULL;
	s->n = n;
	if (n <= 0) return s;

	// ── COO -> CSC (列ごとに集約、重複は加算) ──────────────────────
	int *ap = (int *)calloc((size_t)(n + 1), sizeof(int));
	int *ai = (int *)malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof(int));
	d_complex_t *ax = (d_complex_t *)malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof(d_complex_t));
	int *mark = (int *)malloc((size_t)n * sizeof(int));
	int *hpos = (int *)malloc((size_t)n * sizeof(int));
	if ((ap == NULL) || (ai == NULL) || (ax == NULL) || (mark == NULL) || (hpos == NULL)) {
		free(ap); free(ai); free(ax); free(mark); free(hpos);
		slu_free(s);
		return NULL;
	}
	{
		// 列ごとの個数 (重複込み) を数えて詰め、そのあと列内で重複を潰す
		int *cnt = (int *)calloc((size_t)(n + 1), sizeof(int));
		int *tri = (int *)malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof(int));
		d_complex_t *tvl = (d_complex_t *)malloc((size_t)(nnz > 0 ? nnz : 1) * sizeof(d_complex_t));
		if ((cnt == NULL) || (tri == NULL) || (tvl == NULL)) {
			free(cnt); free(tri); free(tvl);
			free(ap); free(ai); free(ax); free(mark); free(hpos);
			slu_free(s);
			return NULL;
		}
		for (int e = 0; e < nnz; e++) cnt[ci[e] + 1]++;
		for (int j = 0; j < n; j++) cnt[j + 1] += cnt[j];
		int *pos = (int *)malloc((size_t)(n + 1) * sizeof(int));
		if (pos == NULL) {
			free(cnt); free(tri); free(tvl);
			free(ap); free(ai); free(ax); free(mark); free(hpos);
			slu_free(s);
			return NULL;
		}
		memcpy(pos, cnt, (size_t)(n + 1) * sizeof(int));
		for (int e = 0; e < nnz; e++) {
			const int j = ci[e];
			tri[pos[j]] = ri[e];
			tvl[pos[j]] = val[e];
			pos[j]++;
		}
		for (int i = 0; i < n; i++) mark[i] = -1;
		int m = 0;
		for (int j = 0; j < n; j++) {
			ap[j] = m;
			for (int e = cnt[j]; e < cnt[j + 1]; e++) {
				const int i = tri[e];
				if (mark[i] == j) {
					ax[hpos[i]] = d_add(ax[hpos[i]], tvl[e]);
				}
				else {
					mark[i] = j;
					hpos[i] = m;
					ai[m] = i;
					ax[m] = tvl[e];
					m++;
				}
			}
		}
		ap[n] = m;
		free(cnt); free(tri); free(tvl); free(pos);
	}

	// ── 充填低減順序 ────────────────────────────────────────────
	s->q = (int *)malloc((size_t)n * sizeof(int));
	s->pinv = (int *)malloc((size_t)n * sizeof(int));
	s->w = (d_complex_t *)malloc((size_t)n * sizeof(d_complex_t));
	s->xi = (int *)malloc((size_t)(2 * n) * sizeof(int));
	s->lp = (int *)malloc((size_t)(n + 1) * sizeof(int));
	s->up = (int *)malloc((size_t)(n + 1) * sizeof(int));
	if ((s->q == NULL) || (s->pinv == NULL) || (s->w == NULL) || (s->xi == NULL)
	 || (s->lp == NULL) || (s->up == NULL)) {
		free(ap); free(ai); free(ax); free(mark); free(hpos);
		slu_free(s);
		return NULL;
	}
	if (mindeg_order(n, ap, ai, s->q)) {
		for (int k = 0; k < n; k++) s->q[k] = k;   // 打ち切り : 自然順
	}

	// ── Gilbert-Peierls 左向き LU + 閾値部分ピボット ─────────────
	int lcap = 0, ucap = 0, lnz = 0, unz = 0;
	int bad = 0;
	for (int i = 0; i < n; i++) {
		s->pinv[i] = -1;
		mark[i] = -1;
	}

	for (int k = 0; k < n; k++) {
		s->lp[k] = lnz;
		s->up[k] = unz;
		if (grow(&s->li, &s->lx, lnz + n, &lcap)
		 || grow(&s->ui, &s->ux, unz + n, &ucap)) { bad = 1; break; }

		const int col = s->q[k];
		const int top = spsolve(n, s->lp, s->li, s->lx, s->pinv,
			ap, ai, ax, col, s->xi, mark, k, s->w);

		// ピボット選択 : 未確定行の最大値。ただし対角 (行 col) が
		// SLU_TOL 倍以上あれば対角を採る (充填を抑える)
		int ipiv = -1;
		double amax = -1;
		for (int p = top; p < n; p++) {
			const int i = s->xi[p];
			if (s->pinv[i] >= 0) continue;
			const double a = mag(s->w[i]);
			if (a > amax) {
				amax = a;
				ipiv = i;
			}
		}
		if (ipiv < 0) { bad = 1; break; }
		// 対角 (行 col) がパターンにあり、最大値の SLU_TOL 倍以上あれば
		// そちらを採る。mark[col] == k がパターンに入っている印 (これを
		// 見ないと w[col] が前の列の残骸になる)
		if ((mark[col] == k) && (s->pinv[col] < 0)
		 && (mag(s->w[col]) >= SLU_TOL * amax) && (mag(s->w[col]) > 0)) {
			ipiv = col;
		}
		const d_complex_t piv = s->w[ipiv];
		if (mag(piv) < 1e-300) { bad = 1; break; }

		s->pinv[ipiv] = k;
		const d_complex_t pinvz = d_inv(piv);
		for (int p = top; p < n; p++) {
			const int i = s->xi[p];
			if (s->pinv[i] < 0) {
				// L : 行番号は未確定なので元の行番号で入れ、最後に直す
				s->li[lnz] = i;
				s->lx[lnz] = d_mul(s->w[i], pinvz);
				lnz++;
			}
			else if (s->pinv[i] < k) {
				// pinv[i] == k はピボット自身なので除く (対角は下で入れる)
				s->ui[unz] = s->pinv[i];
				s->ux[unz] = s->w[i];
				unz++;
			}
		}
		// U の対角は列の末尾に置く (後退代入がここを見る)
		s->ui[unz] = k;
		s->ux[unz] = piv;
		unz++;
	}
	s->lp[n] = lnz;
	s->up[n] = unz;

	free(ap); free(ai); free(ax); free(mark); free(hpos);
	if (bad) {
		slu_free(s);
		return NULL;
	}

	// L の行番号をピボット位置に直す (これで L は単位下三角になる)
	for (int e = 0; e < lnz; e++) s->li[e] = s->pinv[s->li[e]];

	return s;
}

// b <- A^-1 b
void slu_solve(const struct sparse_lu_t *s, d_complex_t *b)
{
	const int n = s->n;
	if (n <= 0) return;
	d_complex_t *x = s->w;

	// 行の置換 : x[pinv[i]] = b[i]
	for (int i = 0; i < n; i++) x[s->pinv[i]] = b[i];

	// 前進代入 (L は単位下三角、対角は格納していない)
	for (int k = 0; k < n; k++) {
		for (int e = s->lp[k]; e < s->lp[k + 1]; e++) {
			x[s->li[e]] = d_sub(x[s->li[e]], d_mul(s->lx[e], x[k]));
		}
	}

	// 後退代入 (U の対角は各列の末尾)
	for (int k = n - 1; k >= 0; k--) {
		const int d = s->up[k + 1] - 1;
		x[k] = d_div(x[k], s->ux[d]);
		for (int e = s->up[k]; e < d; e++) {
			x[s->ui[e]] = d_sub(x[s->ui[e]], d_mul(s->ux[e], x[k]));
		}
	}

	// 列の置換 : b[q[k]] = x[k]
	for (int k = 0; k < n; k++) b[s->q[k]] = x[k];
}

// LU の非ゼロ数 (メモリ表示用)
size_t slu_nnz(const struct sparse_lu_t *s)
{
	if ((s == NULL) || (s->n <= 0)) return 0;

	return (size_t)s->lp[s->n] + (size_t)s->up[s->n];
}
