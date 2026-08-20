/*
hmatrix.c

部分インダクタンス行列の階層的低ランク圧縮 (compression = 1 のときのみ)

密行列は O(N^2) のメモリと充填を要する。遠く離れたセル群どうしの相互作用は
カーネル 1/R が滑らかなので数値的に低ランクになる。これを使って

  1. クラスタツリー : セルをバウンディングボックスの最長軸で二分し、
     葉のセル数が HM_LEAF 以下になるまで再帰する。セルは perm[] で
     並べ替えられ、各クラスタは perm 上の連続区間になる。
  2. ブロッククラスタツリー : 行クラスタ s と列クラスタ t の対を再帰的に
     見て、許容条件 (admissibility)
        min(diam(s), diam(t)) <= HM_ETA * dist(s, t)
     を満たせば「遠方」= 低ランク近似、満たさず片方が葉なら「近傍」= 密。
  3. ACA (Adaptive Cross Approximation) : 遠方ブロックを部分ピボットで
     A ~ U V^T (rank k) に近似する。要素は lp_entry() を呼んで必要な行・列
     だけ評価するので、ブロック全体を作らずに済む (充填も O(k(m+n)))。

密行列との差は ACA の許容誤差 ctol で制御する。近傍ブロックは密のままなので
自己項・隣接項 (物理の不変条件が効く部分) は一切近似されない。

スレッド数不変性 : 充填はブロックごとに独立 (共有配列への += なし)。
行列ベクトル積は行クラスタごとにグループ化して並列化するので、各スレッドは
出力の互いに素な区間だけを書き、1 つの区間内のブロック加算順は固定される。
したがって結果はスレッド数によらずビット単位で一致する。
*/

#include "peec.h"

// 要素コールバックの Lp 版 (ctx = peec_t *)
static d_complex_t lp_entry_cb(void *ctx, int i, int j, double kw)
{
	return lp_entry((const peec_t *)ctx, i, j, kw);
}

#define HM_LEAF   32      // 葉クラスタの最大セル数
#define HM_ETA    2.0     // 許容条件のパラメータ (大きいほど圧縮が積極的)
#define HM_MAXIT  200     // ACA の最大ランク (超えたら密に落とす)
#define HM_PROBE  4       // 収束後に「近似に現れていない行」を試す回数の上限

// ── クラスタツリー ─────────────────────────────────────────────
typedef struct {
	int    *lo, *hi;          // [nc] perm 上の区間 [lo, hi)
	int    *left, *right;     // [nc] 子 (-1 = 葉)
	double *cen;              // [3*nc] バウンディングボックス中心
	double *rad;              // [nc] 対角の半分 (diam = 2*rad)
	int     nc, cap;
} ctree_t;

struct hmat_t {
	int      n;               // 行列サイズ (= nseg)
	int     *perm;            // [n] クラスタ順 -> 元のセル番号
	ctree_t  ct;
	// 遠方ブロック (低ランク U V^T)
	int      nfar, fcap;
	int     *fr, *fc, *frank;
	d_complex_t **fu, **fv;   // fu[b] : [m x k] 行優先、fv[b] : [nn x k]
	// 近傍ブロック (密、行優先 [m x nn])
	int      nnear, ncap;
	int     *nr, *nc2;
	d_complex_t **nd;
	// 葉行クラスタごとのブロック並び (matvec のスレッド不変な並列化用)
	int      nleaf;           // 葉クラスタの個数 ([0, n) を重複なく分割する)
	int     *leafid;          // [nleaf] 葉クラスタ id (perm 上の位置順)
	int     *rfoff, *rfidx;   // 葉ごとの遠方ブロック index (CSR 風)
	int     *rnoff, *rnidx;   // 葉ごとの近傍ブロック index
	int     *ftoff;           // [nfar+1] 遠方ブロックのランク前置和 (V^T x の置き場)
	int      ftlen;
};

static void ct_reserve(ctree_t *t, int need)
{
	if (need < t->cap) return;
	t->cap = t->cap ? (2 * t->cap) : 64;
	if (t->cap < need + 1) t->cap = need + 1;
	t->lo    = (int *)realloc(t->lo, (size_t)t->cap * sizeof(int));
	t->hi    = (int *)realloc(t->hi, (size_t)t->cap * sizeof(int));
	t->left  = (int *)realloc(t->left, (size_t)t->cap * sizeof(int));
	t->right = (int *)realloc(t->right, (size_t)t->cap * sizeof(int));
	t->cen   = (double *)realloc(t->cen, (size_t)t->cap * 3 * sizeof(double));
	t->rad   = (double *)realloc(t->rad, (size_t)t->cap * sizeof(double));
}


/*
クラスタを作る (perm[lo..hi) を最長軸で二分)。戻り値はクラスタ id。
分割点は中点分割 (バウンディングボックス中央) で、偏った場合は中央値に
落とす。どちらも決定的なのでスレッド数・実行回数に依らず同じ木になる。
*/
static int ct_build(ctree_t *t, const double *pts, int *perm, int lo, int hi)
{
	ct_reserve(t, t->nc);
	const int id = t->nc++;
	t->lo[id] = lo;
	t->hi[id] = hi;
	t->left[id] = t->right[id] = -1;

	// バウンディングボックス
	double bmin[3], bmax[3];
	for (int c = 0; c < 3; c++) {
		bmin[c] = bmax[c] = pts[(3 * perm[lo]) + c];
	}
	for (int i = lo + 1; i < hi; i++) {
		const double *x = &pts[3 * perm[i]];
		for (int c = 0; c < 3; c++) {
			if (x[c] < bmin[c]) bmin[c] = x[c];
			if (x[c] > bmax[c]) bmax[c] = x[c];
		}
	}
	double d2 = 0;
	for (int c = 0; c < 3; c++) {
		t->cen[(3 * id) + c] = 0.5 * (bmin[c] + bmax[c]);
		const double d = bmax[c] - bmin[c];
		d2 += d * d;
	}
	t->rad[id] = 0.5 * sqrt(d2);

	if ((hi - lo) <= HM_LEAF) return id;

	// 最長軸
	int ax = 0;
	for (int c = 1; c < 3; c++) {
		if ((bmax[c] - bmin[c]) > (bmax[ax] - bmin[ax])) ax = c;
	}
	if ((bmax[ax] - bmin[ax]) <= 0) return id;    // 全点が同一 : これ以上割れない

	// 中点で分割 (in-place パーティション)
	const double cut = 0.5 * (bmin[ax] + bmax[ax]);
	int i = lo, j = hi - 1;
	while (i <= j) {
		if (pts[(3 * perm[i]) + ax] < cut) {
			i++;
		}
		else {
			const int tmp = perm[i];
			perm[i] = perm[j];
			perm[j] = tmp;
			j--;
		}
	}
	int mid = i;
	// 片側が空になったら中央で割る (決定的なフォールバック)
	if ((mid == lo) || (mid == hi)) mid = (lo + hi) / 2;

	const int l = ct_build(t, pts, perm, lo, mid);
	const int r = ct_build(t, pts, perm, mid, hi);
	t->left[id] = l;
	t->right[id] = r;

	return id;
}

// クラスタ間の距離 (バウンディング球の隙間。重なりは 0)
static double ct_dist(const ctree_t *t, int a, int b)
{
	double d2 = 0;
	for (int c = 0; c < 3; c++) {
		const double d = t->cen[(3 * a) + c] - t->cen[(3 * b) + c];
		d2 += d * d;
	}
	const double d = sqrt(d2) - t->rad[a] - t->rad[b];

	return (d > 0) ? d : 0;
}

// 許容条件 : min(diam) <= eta * dist
static int ct_admissible(const ctree_t *t, int a, int b)
{
	const double dmin = 2 * ((t->rad[a] < t->rad[b]) ? t->rad[a] : t->rad[b]);
	const double dist = ct_dist(t, a, b);

	return (dist > 0) && (dmin <= HM_ETA * dist);
}

// ── ブロックリストの構築 ───────────────────────────────────────
static void hm_add_far(struct hmat_t *h, int a, int b)
{
	if (h->nfar >= h->fcap) {
		h->fcap = h->fcap ? (2 * h->fcap) : 256;
		h->fr    = (int *)realloc(h->fr, (size_t)h->fcap * sizeof(int));
		h->fc    = (int *)realloc(h->fc, (size_t)h->fcap * sizeof(int));
		h->frank = (int *)realloc(h->frank, (size_t)h->fcap * sizeof(int));
		h->fu    = (d_complex_t **)realloc(h->fu, (size_t)h->fcap * sizeof(d_complex_t *));
		h->fv    = (d_complex_t **)realloc(h->fv, (size_t)h->fcap * sizeof(d_complex_t *));
	}
	h->fr[h->nfar] = a;
	h->fc[h->nfar] = b;
	h->frank[h->nfar] = 0;
	h->fu[h->nfar] = NULL;
	h->fv[h->nfar] = NULL;
	h->nfar++;
}

static void hm_add_near(struct hmat_t *h, int a, int b)
{
	if (h->nnear >= h->ncap) {
		h->ncap = h->ncap ? (2 * h->ncap) : 256;
		h->nr  = (int *)realloc(h->nr, (size_t)h->ncap * sizeof(int));
		h->nc2 = (int *)realloc(h->nc2, (size_t)h->ncap * sizeof(int));
		h->nd  = (d_complex_t **)realloc(h->nd, (size_t)h->ncap * sizeof(d_complex_t *));
	}
	h->nr[h->nnear] = a;
	h->nc2[h->nnear] = b;
	h->nd[h->nnear] = NULL;
	h->nnear++;
}

// ブロッククラスタツリーを再帰的に降りて遠方 / 近傍に振り分ける
static void hm_split(struct hmat_t *h, int a, int b)
{
	const ctree_t *t = &h->ct;

	if (ct_admissible(t, a, b)) {
		hm_add_far(h, a, b);
		return;
	}
	// どちらかが葉なら密ブロックにする
	if ((t->left[a] < 0) || (t->left[b] < 0)) {
		hm_add_near(h, a, b);
		return;
	}
	hm_split(h, t->left[a],  t->left[b]);
	hm_split(h, t->left[a],  t->right[b]);
	hm_split(h, t->right[a], t->left[b]);
	hm_split(h, t->right[a], t->right[b]);
}

// ── ACA (部分ピボット) ────────────────────────────────────────
/*
ブロック A(m x nn) を A ~ sum_k u_k v_k^T に近似する。
行 i* を選んで残差行を取り、その最大成分の列 j* から残差列を取る、を繰り返す
古典的な部分ピボット ACA。停止判定は
  |u_k| |v_k| <= tol * ||A_k||_F        (A_k = ここまでの近似)
戻り値 : ランク k (0 < k)。HM_MAXIT を超えたら -1 (呼び出し側が密に落とす)。

**可約なブロックへの対策** : Lp の要素は方向余弦 (t_i・t_j) を含むので、
直交するセルどうしの要素は**厳密に 0** になる。格子状の面導体では x 向きと
y 向きのセルが混在するため、素朴な部分ピボット ACA は x 群だけを追いかけて
残差列が y 群の行で恒等的に 0 になり、**y-y 部分ブロックを一度も見ないまま
「収束した」と判断する** (実際に踏んだ : ランクを 8 → 13 → 23 と増やしても
ブロック相対誤差が 0.76 のまま動かない)。

そこで各行が u_l に現れた最大の大きさ urep[i] を記録し、収束したように
見えたときに**まだ一度も現れていない行**が残っていればそこから再開する。
可約でないブロックでは urep はすべて有限なので余分な行評価は起きない
(面導体でも群の切り替わりの 1 行だけ)。
*/
// 近似にまだ現れていない未使用行を選ぶ (無ければ -1)
static int aca_unrepresented(const int *rused, const double *urep, int m, double thr)
{
	int nxt = -1;
	double best = 0;
	for (int i = 0; i < m; i++) {
		if (rused[i]) continue;
		if ((nxt < 0) || (urep[i] < best)) {
			best = urep[i];
			nxt = i;
		}
	}

	return ((nxt >= 0) && (best <= thr)) ? nxt : -1;
}

static int aca_block(hmat_entry_fn fn, void *ctx, double kw, const int *perm,
	int r0, int m, int c0, int nn, double tol,
	d_complex_t **uout, d_complex_t **vout)
{
	const int kmax = (m < nn) ? m : nn;
	const int lim = (kmax < HM_MAXIT) ? kmax : HM_MAXIT;
	d_complex_t *u = (d_complex_t *)malloc((size_t)m * lim * sizeof(d_complex_t));
	d_complex_t *v = (d_complex_t *)malloc((size_t)nn * lim * sizeof(d_complex_t));
	int *rused = (int *)calloc((size_t)m, sizeof(int));
	double *urep = (double *)calloc((size_t)m, sizeof(double));
	if ((u == NULL) || (v == NULL) || (rused == NULL) || (urep == NULL)) {
		free(u); free(v); free(rused); free(urep);
		return -1;
	}

	double fro2 = 0;               // ||A_k||_F^2
	double urmax = 0;              // max_l max_i |u_l[i]|
	int k = 0;
	int irow = 0;
	int probe = 0;                 // 連続して空振りした再開試行の回数

	while (k < lim) {
		// 残差行 v = A(irow, :) - sum_l u_l[irow] v_l
		// (同時に行ノルムと最大成分の列を取る)
		int jmax = -1;
		double amax = 0, rn2 = 0;
		for (int j = 0; j < nn; j++) {
			d_complex_t a = fn(ctx, perm[r0 + irow], perm[c0 + j], kw);
			for (int l = 0; l < k; l++) {
				a = d_sub(a, d_mul(u[(l * m) + irow], v[(l * nn) + j]));
			}
			v[(k * nn) + j] = a;
			const double t = d_abs(a);
			rn2 += t * t;
			if (t > amax) {
				amax = t;
				jmax = j;
			}
		}
		rused[irow] = 1;

		if ((jmax < 0) || (amax <= 0) || (sqrt(rn2) <= tol * sqrt(fro2))) {
			// この行には有効な残差が無い : 近似に現れていない行から再開する
			if (++probe > HM_PROBE) break;
			const int nxt = aca_unrepresented(rused, urep, m, tol * urmax);
			if (nxt < 0) break;
			irow = nxt;
			continue;
		}
		probe = 0;

		// 行を正規化
		const d_complex_t piv = v[(k * nn) + jmax];
		for (int j = 0; j < nn; j++) {
			v[(k * nn) + j] = d_div(v[(k * nn) + j], piv);
		}

		// 残差列 u = A(:, jmax) - sum_l u_l v_l[jmax]
		for (int i = 0; i < m; i++) {
			d_complex_t a = fn(ctx, perm[r0 + i], perm[c0 + jmax], kw);
			for (int l = 0; l < k; l++) {
				a = d_sub(a, d_mul(u[(l * m) + i], v[(l * nn) + jmax]));
			}
			u[(k * m) + i] = a;
		}

		// ノルムと停止判定 (||u_k|| ||v_k|| が近似全体に対して十分小さいか)
		double un = 0, vn = 0;
		for (int i = 0; i < m; i++) un += d_norm(u[(k * m) + i]);
		for (int j = 0; j < nn; j++) vn += d_norm(v[(k * nn) + j]);
		un = sqrt(un);
		vn = sqrt(vn);

		// ||A_k||_F^2 の更新 (交差項込み)
		double cross = 0;
		for (int l = 0; l < k; l++) {
			d_complex_t du = d_complex(0, 0), dv = d_complex(0, 0);
			for (int i = 0; i < m; i++) {
				du = d_add(du, d_mul(d_conj(u[(l * m) + i]), u[(k * m) + i]));
			}
			for (int j = 0; j < nn; j++) {
				dv = d_add(dv, d_mul(d_conj(v[(l * nn) + j]), v[(k * nn) + j]));
			}
			cross += 2 * d_abs(d_mul(du, dv));
		}
		fro2 += (un * un * vn * vn) + cross;
		k++;

		// 各行が近似にどれだけ現れたか (可約なブロックの検出用)
		for (int i = 0; i < m; i++) {
			const double t = d_abs(u[((k - 1) * m) + i]);
			if (t > urep[i]) urep[i] = t;
			if (t > urmax) urmax = t;
		}

		// 次のピボット行 : u の最大成分 (未使用のもの)
		int inext = -1;
		if ((un * vn) > tol * sqrt((fro2 > 0) ? fro2 : 1)) {
			double umax = 0;
			for (int i = 0; i < m; i++) {
				if (rused[i]) continue;
				const double t = d_abs(u[((k - 1) * m) + i]);
				if (t > umax) {
					umax = t;
					inext = i;
				}
			}
		}
		if (inext < 0) {
			// 収束した、または残差列が未使用行で恒等的に 0 (直交群の混在)。
			// 近似に現れていない行が残っていればそこから再開する。
			if (++probe > HM_PROBE) break;
			inext = aca_unrepresented(rused, urep, m, tol * urmax);
			if (inext < 0) break;
		}
		irow = inext;
	}

	free(rused);
	free(urep);
	if ((k <= 0) || (k >= lim)) {
		// 収束しなかった (低ランクでない) : 呼び出し側で密にする
		free(u);
		free(v);
		return -1;
	}
	*uout = u;
	*vout = v;

	return k;
}

// ── 構築 ──────────────────────────────────────────────────────
/*
汎用の構築 : カーネルは要素コールバック fn(ctx, i, j, kw)、幾何は代表点の
配列 pts[3n] だけで決まる。これで部分インダクタンス Lp (区間) と電位係数 P
(電荷セル) の両方を同じ機構で圧縮できる。
*/
struct hmat_t *hmat_build_gen(int n, const double *pts,
	hmat_entry_fn fn, void *ctx, double kw, double tol,
	const char *label, FILE *fp_log)
{
	if (n <= 0) return NULL;

	struct hmat_t *h = (struct hmat_t *)calloc(1, sizeof(struct hmat_t));
	if (h == NULL) return NULL;
	h->n = n;
	h->perm = (int *)malloc((size_t)n * sizeof(int));
	if (h->perm == NULL) {
		free(h);
		return NULL;
	}
	for (int i = 0; i < n; i++) {
		h->perm[i] = i;
	}
	ct_build(&h->ct, pts, h->perm, 0, n);
	hm_split(h, 0, 0);

	// 遠方ブロックを ACA で圧縮 (ブロックごとに独立 = 並列化してよい)
	// MSVC の OpenMP 2.0 のためループ変数は事前宣言
	int b;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (b = 0; b < h->nfar; b++) {
		const int a = h->fr[b], c = h->fc[b];
		const int r0 = h->ct.lo[a], m = h->ct.hi[a] - r0;
		const int c0 = h->ct.lo[c], nn = h->ct.hi[c] - c0;
		d_complex_t *u = NULL, *v = NULL;
		const int k = aca_block(fn, ctx, kw, h->perm, r0, m, c0, nn, tol, &u, &v);
		if (k > 0) {
			h->frank[b] = k;
			h->fu[b] = u;
			h->fv[b] = v;
		}
		else {
			// 低ランクにならなかった : 密で持つ (frank = 0 の印)
			h->frank[b] = 0;
			d_complex_t *d = (d_complex_t *)malloc((size_t)m * nn * sizeof(d_complex_t));
			if (d != NULL) {
				for (int i = 0; i < m; i++) {
					for (int j = 0; j < nn; j++) {
						d[(i * nn) + j] = fn(ctx, h->perm[r0 + i], h->perm[c0 + j], kw);
					}
				}
			}
			h->fu[b] = d;
			h->fv[b] = NULL;
		}
	}

	// 近傍ブロックは密 (自己項・隣接項は一切近似しない)
	int q;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (q = 0; q < h->nnear; q++) {
		const int a = h->nr[q], c = h->nc2[q];
		const int r0 = h->ct.lo[a], m = h->ct.hi[a] - r0;
		const int c0 = h->ct.lo[c], nn = h->ct.hi[c] - c0;
		d_complex_t *d = (d_complex_t *)malloc((size_t)m * nn * sizeof(d_complex_t));
		if (d != NULL) {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < nn; j++) {
					d[(i * nn) + j] = fn(ctx, h->perm[r0 + i], h->perm[c0 + j], kw);
				}
			}
		}
		h->nd[q] = d;
	}

	/*
	matvec 用にブロックを「葉行クラスタ」ごとにまとめる。

	ブロックの行クラスタは木の途中の節点でもよいので、**異なる階層のブロック
	どうしで行区間が重なる**。行クラスタ単位で並列化すると複数スレッドが y の
	同じ要素を read-modify-write して静かに寄与を落とす (実際に踏んだ : 840
	セルで相対誤差 ~0.5、実行ごとに値が変わる)。葉クラスタは perm 上の [0, n)
	を重複なく分割するので、葉ごとに並列化すれば各スレッドは y の互いに素な
	区間だけを書く。各葉の中のブロック順は登録順で固定なので加算順も決定的。
	*/
	int *leafat = (int *)malloc((size_t)n * sizeof(int));
	int *lidx   = (int *)malloc((size_t)n * sizeof(int));
	if ((leafat == NULL) || (lidx == NULL)) {
		free(leafat); free(lidx);
		hmat_free(h);
		return NULL;
	}
	for (int i = 0; i < n; i++) {
		leafat[i] = -1;
	}
	for (int i = 0; i < h->ct.nc; i++) {
		if (h->ct.left[i] < 0) leafat[h->ct.lo[i]] = i;
	}
	h->nleaf = 0;
	for (int i = 0; i < n; i++) {
		if (leafat[i] >= 0) h->nleaf++;
	}
	h->leafid = (int *)malloc((size_t)h->nleaf * sizeof(int));
	h->rfoff  = (int *)calloc((size_t)(h->nleaf + 1), sizeof(int));
	h->rnoff  = (int *)calloc((size_t)(h->nleaf + 1), sizeof(int));
	if ((h->leafid == NULL) || (h->rfoff == NULL) || (h->rnoff == NULL)) {
		free(leafat); free(lidx);
		hmat_free(h);
		return NULL;
	}
	{
		int e = 0;
		for (int i = 0; i < n; i++) {
			if (leafat[i] >= 0) h->leafid[e++] = leafat[i];
		}
		e = -1;
		for (int i = 0; i < n; i++) {
			if (leafat[i] >= 0) e++;
			lidx[i] = e;
		}
	}
	free(leafat);

	// 葉が [0, n) を隙間なく重複なく敷き詰めていることの自己検査。
	// (これが崩れると matvec が競合するか寄与を落とす)
	{
		int bad = (h->nleaf <= 0) || (h->ct.lo[h->leafid[0]] != 0)
			|| (h->ct.hi[h->leafid[h->nleaf - 1]] != n);
		for (int e = 0; !bad && (e + 1 < h->nleaf); e++) {
			if (h->ct.hi[h->leafid[e]] != h->ct.lo[h->leafid[e + 1]]) bad = 1;
		}
		if (bad) {
			if (fp_log != NULL) {
				fprintf(fp_log, "*** compression : leaf clusters do not tile [0, %d)\n", n);
			}
			free(lidx);
			hmat_free(h);
			return NULL;
		}
	}

	// 葉ごとのブロック本数 (1 ブロックはその行区間に含まれる全葉に載る)
	for (b = 0; b < h->nfar; b++) {
		const int a = h->fr[b];
		for (int e = lidx[h->ct.lo[a]]; e <= lidx[h->ct.hi[a] - 1]; e++) {
			h->rfoff[e + 1]++;
		}
	}
	for (q = 0; q < h->nnear; q++) {
		const int a = h->nr[q];
		for (int e = lidx[h->ct.lo[a]]; e <= lidx[h->ct.hi[a] - 1]; e++) {
			h->rnoff[e + 1]++;
		}
	}
	for (int e = 0; e < h->nleaf; e++) {
		h->rfoff[e + 1] += h->rfoff[e];
		h->rnoff[e + 1] += h->rnoff[e];
	}
	h->rfidx = (int *)malloc((size_t)((h->rfoff[h->nleaf] > 0) ? h->rfoff[h->nleaf] : 1) * sizeof(int));
	h->rnidx = (int *)malloc((size_t)((h->rnoff[h->nleaf] > 0) ? h->rnoff[h->nleaf] : 1) * sizeof(int));
	if ((h->rfidx == NULL) || (h->rnidx == NULL)) {
		free(lidx);
		hmat_free(h);
		return NULL;
	}
	{
		int *fp = (int *)malloc((size_t)(h->nleaf + 1) * sizeof(int));
		int *np = (int *)malloc((size_t)(h->nleaf + 1) * sizeof(int));
		if ((fp == NULL) || (np == NULL)) {
			free(fp); free(np); free(lidx);
			hmat_free(h);
			return NULL;
		}
		memcpy(fp, h->rfoff, (size_t)(h->nleaf + 1) * sizeof(int));
		memcpy(np, h->rnoff, (size_t)(h->nleaf + 1) * sizeof(int));
		for (b = 0; b < h->nfar; b++) {
			const int a = h->fr[b];
			for (int e = lidx[h->ct.lo[a]]; e <= lidx[h->ct.hi[a] - 1]; e++) {
				h->rfidx[fp[e]++] = b;
			}
		}
		for (q = 0; q < h->nnear; q++) {
			const int a = h->nr[q];
			for (int e = lidx[h->ct.lo[a]]; e <= lidx[h->ct.hi[a] - 1]; e++) {
				h->rnidx[np[e]++] = q;
			}
		}
		free(fp);
		free(np);
	}
	free(lidx);

	// 遠方ブロックの V^T x の置き場 (ランクの前置和)
	h->ftoff = (int *)malloc((size_t)(h->nfar + 1) * sizeof(int));
	if (h->ftoff == NULL) {
		hmat_free(h);
		return NULL;
	}
	h->ftoff[0] = 0;
	for (b = 0; b < h->nfar; b++) {
		h->ftoff[b + 1] = h->ftoff[b] + ((h->frank[b] > 0) ? h->frank[b] : 0);
	}
	h->ftlen = h->ftoff[h->nfar];

	if (fp_log != NULL) {
		const double mb = hmat_memory_mb(h);
		const double dm = hmat_dense_mb(n);
		fprintf(fp_log, "compression %s : %d cells, %d far blocks + %d near blocks, tol = %.1e, %.2f MB (dense %.2f MB, %.1fx)\n",
			label, n, h->nfar, h->nnear, tol, mb, dm, (mb > 0) ? (dm / mb) : 0);
		fflush(fp_log);
	}

	return h;
}

// 部分インダクタンス Lp (区間セル、代表点は区間の中点)
struct hmat_t *hmat_build(peec_t *p, double f, FILE *fp_log)
{
	const int n = p->nseg;
	if (n <= 0) return NULL;

	double *pts = (double *)malloc((size_t)n * 3 * sizeof(double));
	if (pts == NULL) return NULL;
	for (int i = 0; i < n; i++) {
		for (int c = 0; c < 3; c++) {
			pts[(3 * i) + c] = 0.5 * (p->seg[i].x1[c] + p->seg[i].x2[c]);
		}
	}
	const double kw = p->retardation ? (2 * PI * f / C0) : 0;
	struct hmat_t *h = hmat_build_gen(n, pts, lp_entry_cb, p, kw, p->ctol, "Lp", fp_log);
	free(pts);

	return h;
}

void hmat_free(struct hmat_t *h)
{
	if (h == NULL) return;
	for (int b = 0; b < h->nfar; b++) {
		free(h->fu[b]);
		free(h->fv[b]);
	}
	for (int q = 0; q < h->nnear; q++) {
		free(h->nd[q]);
	}
	free(h->fr); free(h->fc); free(h->frank); free(h->fu); free(h->fv);
	free(h->nr); free(h->nc2); free(h->nd);
	free(h->leafid); free(h->rfoff); free(h->rfidx); free(h->rnoff); free(h->rnidx);
	free(h->ftoff);
	free(h->ct.lo); free(h->ct.hi); free(h->ct.left); free(h->ct.right);
	free(h->ct.cen); free(h->ct.rad);
	free(h->perm);
	free(h);
}

double hmat_memory_mb(const struct hmat_t *h)
{
	if (h == NULL) return 0;
	size_t nz = 0;
	for (int b = 0; b < h->nfar; b++) {
		const int m = h->ct.hi[h->fr[b]] - h->ct.lo[h->fr[b]];
		const int nn = h->ct.hi[h->fc[b]] - h->ct.lo[h->fc[b]];
		nz += (h->frank[b] > 0) ? ((size_t)h->frank[b] * (m + nn)) : ((size_t)m * nn);
	}
	for (int q = 0; q < h->nnear; q++) {
		const int m = h->ct.hi[h->nr[q]] - h->ct.lo[h->nr[q]];
		const int nn = h->ct.hi[h->nc2[q]] - h->ct.lo[h->nc2[q]];
		nz += (size_t)m * nn;
	}

	return (double)nz * sizeof(d_complex_t) / (1024.0 * 1024.0);
}

double hmat_dense_mb(int n)
{
	return (double)n * n * sizeof(d_complex_t) / (1024.0 * 1024.0);
}

// ── 行列ベクトル積 ─────────────────────────────────────────────
void hmat_matvec(const struct hmat_t *h, const d_complex_t *x, d_complex_t *y)
{
	const int n = h->n;
	memset(y, 0, (size_t)n * sizeof(d_complex_t));
	if (h->nleaf <= 0) return;

	/*
	第 1 段 : 遠方ブロックごとに t = V^T x を作る。ブロックごとに独立で、
	書き先も互いに素 (ftoff の区間) なので競合しない。第 2 段でブロックの
	行区間を葉に切って処理するため、ここで一度だけ作って共有する。
	*/
	d_complex_t *ft = NULL;
	if (h->ftlen > 0) {
		ft = (d_complex_t *)malloc((size_t)h->ftlen * sizeof(d_complex_t));
		if (ft == NULL) return;
		int b;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
		for (b = 0; b < h->nfar; b++) {
			const int k = h->frank[b];
			if (k <= 0) continue;
			const int c0 = h->ct.lo[h->fc[b]];
			const int nn = h->ct.hi[h->fc[b]] - c0;
			for (int l = 0; l < k; l++) {
				const d_complex_t *vl = &h->fv[b][(size_t)l * nn];
				d_complex_t t = d_complex(0, 0);
				for (int j = 0; j < nn; j++) {
					t = d_add(t, d_mul(vl[j], x[h->perm[c0 + j]]));
				}
				ft[h->ftoff[b] + l] = t;
			}
		}
	}

	// 第 2 段 : 葉行クラスタごとに並列化。葉は perm 上の [0, n) を重複なく
	// 分割するので各スレッドは y の互いに素な区間だけを書き、1 区間内の
	// 加算順は rfidx / rnidx の並び順で固定される。したがってスレッド数を
	// 変えても結果はビット単位で一致する。
	int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (r = 0; r < h->nleaf; r++) {
		const int a = h->leafid[r];
		const int i0 = h->ct.lo[a], i1 = h->ct.hi[a];

		// 遠方 : y += U (V^T x)
		for (int e = h->rfoff[r]; e < h->rfoff[r + 1]; e++) {
			const int b = h->rfidx[e];
			const int r0 = h->ct.lo[h->fr[b]];
			const int m = h->ct.hi[h->fr[b]] - r0;
			const int c0 = h->ct.lo[h->fc[b]];
			const int nn = h->ct.hi[h->fc[b]] - c0;
			const int k = h->frank[b];
			if (k <= 0) {
				// 密で持っている遠方ブロック
				const d_complex_t *d = h->fu[b];
				if (d == NULL) continue;
				for (int i = i0; i < i1; i++) {
					const d_complex_t *row = &d[(size_t)(i - r0) * nn];
					d_complex_t s = d_complex(0, 0);
					for (int j = 0; j < nn; j++) {
						s = d_add(s, d_mul(row[j], x[h->perm[c0 + j]]));
					}
					y[h->perm[i]] = d_add(y[h->perm[i]], s);
				}
				continue;
			}
			for (int l = 0; l < k; l++) {
				const d_complex_t t = ft[h->ftoff[b] + l];
				const d_complex_t *ul = &h->fu[b][(size_t)l * m];
				for (int i = i0; i < i1; i++) {
					y[h->perm[i]] = d_add(y[h->perm[i]], d_mul(ul[i - r0], t));
				}
			}
		}

		// 近傍 : y += A x (密)
		for (int e = h->rnoff[r]; e < h->rnoff[r + 1]; e++) {
			const int q = h->rnidx[e];
			const int r0 = h->ct.lo[h->nr[q]];
			const int c0 = h->ct.lo[h->nc2[q]];
			const int nn = h->ct.hi[h->nc2[q]] - c0;
			const d_complex_t *d = h->nd[q];
			if (d == NULL) continue;
			for (int i = i0; i < i1; i++) {
				const d_complex_t *row = &d[(size_t)(i - r0) * nn];
				d_complex_t s = d_complex(0, 0);
				for (int j = 0; j < nn; j++) {
					s = d_add(s, d_mul(row[j], x[h->perm[c0 + j]]));
				}
				y[h->perm[i]] = d_add(y[h->perm[i]], s);
			}
		}
	}

	free(ft);
}

// 近傍ブロックの要素を列挙する (ILU(0) のパターン構築用、直列)
void hmat_near_each(const struct hmat_t *h,
	void (*cb)(int i, int j, d_complex_t v, void *arg), void *arg)
{
	for (int q = 0; q < h->nnear; q++) {
		const int r0 = h->ct.lo[h->nr[q]];
		const int m = h->ct.hi[h->nr[q]] - r0;
		const int c0 = h->ct.lo[h->nc2[q]];
		const int nn = h->ct.hi[h->nc2[q]] - c0;
		const d_complex_t *d = h->nd[q];
		if (d == NULL) continue;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < nn; j++) {
				cb(h->perm[r0 + i], h->perm[c0 + j], d[(i * nn) + j], arg);
			}
		}
	}
}

/*
葉クラスタの列挙 (前処理の構築用)

葉は perm 上の [0, n) を重複なく分割する。葉 r について
  *idx = そのセル番号の配列 (元の番号)、戻り値 = 対角ブロック (m x m 行優先)
を返す。対角ブロック (a, a) は dist = 0 で許容条件を満たさないため必ず近傍
ブロックとして密に持たれている (= 近似されていない)。
*/
int hmat_nleaf(const struct hmat_t *h)
{
	return (h != NULL) ? h->nleaf : 0;
}

const d_complex_t *hmat_leaf_block(const struct hmat_t *h, int r,
	const int **idx, int *m)
{
	const int a = h->leafid[r];
	*idx = &h->perm[h->ct.lo[a]];
	*m = h->ct.hi[a] - h->ct.lo[a];
	for (int e = h->rnoff[r]; e < h->rnoff[r + 1]; e++) {
		const int q = h->rnidx[e];
		if ((h->nr[q] == a) && (h->nc2[q] == a)) return h->nd[q];
	}

	return NULL;
}
