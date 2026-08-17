/*
iterative.c

周波数掃引の反復解法 (acceleration = 1) と行列フリー解法 (compression = 1)

acceleration = 1 : 掃引では隣接周波数の MNA 行列 A(f) = S + j omega T が
近いことを利用して、直近に LU 分解した行列 M = A(f_LU) を右前処理に使う
GMRES で解く :

  A M^-1 u = b,   x = M^-1 u

収束しない (= 前処理が古い) ときは呼び出し側 (solve.c) が現在の周波数で
LU を取り直す。最悪でも従来の毎周波数 LU と同等で、掃引点が密なほど
LU の回数が減り、O(n^3) が O(反復 x n^2) に置き換わる。

compression = 1 : 密行列を一切持たず、行列ベクトル積を疎部 (mna_apply) +
H 行列 (hmat_matvec) で、前処理を葉ブロック消去 + Schur 補元 (precond.c)
で行う。前処理が周波数ごとに作り直されるので通常は 1 サイクルで収束するが、
念のため GMRES(GMRES_MAXIT) のリスタートを GMRES_CYCLES 回まで許す
(基底のメモリを GMRES_MAXIT x n に抑えるため)。

スレッド数不変性 : 行列ベクトル積は行ごとに独立 (各行の内積は直列)、
Gram-Schmidt / ノルム / Givens 回転はすべて直列なので、結果はスレッド数に
依らずビット単位で一致する。順序依存の並列リダクションは入れないこと。
*/

#include "peec.h"

#define GMRES_MAXIT  60
#define GMRES_CYCLES 10

// conj(u)・v (複素内積、直列)
static d_complex_t zdotc(int n, const d_complex_t *u, const d_complex_t *v)
{
	double sr = 0, si = 0;
	for (int i = 0; i < n; i++) {
		sr += (u[i].r * v[i].r) + (u[i].i * v[i].i);
		si += (u[i].r * v[i].i) - (u[i].i * v[i].r);
	}
	return d_complex(sr, si);
}

static double znorm(int n, const d_complex_t *u)
{
	double s = 0;
	for (int i = 0; i < n; i++) {
		s += (u[i].r * u[i].r) + (u[i].i * u[i].i);
	}
	return sqrt(s);
}

// y = A x (密行列。行ごとに独立なので並列化してもビット一致)
static void zmatvec(int n, const d_complex_t *a, const d_complex_t *x, d_complex_t *y)
{
	// MSVC の OpenMP 2.0 のためループ変数は事前宣言
	int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (i = 0; i < n; i++) {
		const d_complex_t *row = &a[(size_t)i * n];
		double sr = 0, si = 0;
		for (int j = 0; j < n; j++) {
			sr += (row[j].r * x[j].r) - (row[j].i * x[j].i);
			si += (row[j].r * x[j].i) + (row[j].i * x[j].r);
		}
		y[i] = d_complex(sr, si);
	}
}

/*
右前処理 GMRES の 1 サイクル (最大 GMRES_MAXIT 反復)

  aop(actx, x, y) : y = A x
  mop(mctx, z)    : z <- M^-1 z (in-place)
  r    : 初期残差 (呼び出し側が b - A x0 を渡す)
  x    : 出力。このサイクルの補正 M^-1 V y を**常に**書く (未収束でも
         最小二乗解を返すので、リスタートの初期値に使える)
  tolabs : 残差の絶対収束判定 (相対判定は呼び出し側が ||b|| を掛けて渡す)

戻り値 : 反復回数 (>= 1、収束時) / -1 = このサイクルでは未収束
*/
static int gmres_cycle(int n,
	void (*aop)(void *ctx, const d_complex_t *x, d_complex_t *y), void *actx,
	void (*mop)(void *ctx, d_complex_t *z), void *mctx,
	const d_complex_t *r, d_complex_t *x, double tolabs)
{
	const int m = GMRES_MAXIT;
	// Krylov 基底 V[(m+1) x n]、Hessenberg h[(m+1) x m]、Givens 係数、右辺 g
	d_complex_t *v = (d_complex_t *)malloc((size_t)(m + 1) * n * sizeof(d_complex_t));
	d_complex_t *w = (d_complex_t *)malloc((size_t)n * sizeof(d_complex_t));
	d_complex_t *z = (d_complex_t *)malloc((size_t)n * sizeof(d_complex_t));
	d_complex_t *h = (d_complex_t *)calloc((size_t)(m + 1) * m, sizeof(d_complex_t));
	d_complex_t *cs = (d_complex_t *)malloc((size_t)m * sizeof(d_complex_t));
	d_complex_t *sn = (d_complex_t *)malloc((size_t)m * sizeof(d_complex_t));
	d_complex_t *g = (d_complex_t *)calloc((size_t)(m + 1), sizeof(d_complex_t));
	d_complex_t *y = (d_complex_t *)malloc((size_t)m * sizeof(d_complex_t));
	if ((v == NULL) || (w == NULL) || (z == NULL) || (h == NULL)
	 || (cs == NULL) || (sn == NULL) || (g == NULL) || (y == NULL)) {
		free(v); free(w); free(z); free(h); free(cs); free(sn); free(g); free(y);
		return -1;
	}

	const double beta = znorm(n, r);
	int iters = -1;
	if (beta <= 0) {
		// 残差 0 : 補正も 0
		memset(x, 0, (size_t)n * sizeof(d_complex_t));
		free(v); free(w); free(z); free(h); free(cs); free(sn); free(g); free(y);
		return 1;
	}
	for (int i = 0; i < n; i++) {
		v[i] = d_rmul(1 / beta, r[i]);
	}
	g[0] = d_complex(beta, 0);

	int j;
	for (j = 0; j < m; j++) {
		// w = A M^-1 v_j
		memcpy(z, &v[(size_t)j * n], (size_t)n * sizeof(d_complex_t));
		mop(mctx, z);
		aop(actx, z, w);

		// 修正 Gram-Schmidt (直列)
		for (int i = 0; i <= j; i++) {
			const d_complex_t hij = zdotc(n, &v[(size_t)i * n], w);
			h[i + ((size_t)j * (m + 1))] = hij;
			for (int k = 0; k < n; k++) {
				w[k] = d_sub(w[k], d_mul(hij, v[(size_t)i * n + k]));
			}
		}
		const double hn = znorm(n, w);
		h[(j + 1) + ((size_t)j * (m + 1))] = d_complex(hn, 0);
		if (hn > 0) {
			for (int k = 0; k < n; k++) {
				v[(size_t)(j + 1) * n + k] = d_rmul(1 / hn, w[k]);
			}
		}

		// これまでの Givens 回転を新しい列に適用
		for (int i = 0; i < j; i++) {
			const d_complex_t h1 = h[i + ((size_t)j * (m + 1))];
			const d_complex_t h2 = h[(i + 1) + ((size_t)j * (m + 1))];
			// [h1'; h2'] = [[conj(c), conj(s)], [-s, c]] [h1; h2]
			h[i + ((size_t)j * (m + 1))] = d_add(
				d_mul(d_complex(cs[i].r, -cs[i].i), h1),
				d_mul(d_complex(sn[i].r, -sn[i].i), h2));
			h[(i + 1) + ((size_t)j * (m + 1))] = d_sub(
				d_mul(cs[i], h2), d_mul(sn[i], h1));
		}

		// 新しい回転で h[j+1][j] を消す
		const d_complex_t h1 = h[j + ((size_t)j * (m + 1))];
		const d_complex_t h2 = h[(j + 1) + ((size_t)j * (m + 1))];
		const double rho = sqrt((h1.r * h1.r) + (h1.i * h1.i)
		                      + (h2.r * h2.r) + (h2.i * h2.i));
		if (rho <= 0) break;                            // 破綻 (既に厳密解)
		cs[j] = d_rmul(1 / rho, h1);
		sn[j] = d_rmul(1 / rho, h2);
		h[j + ((size_t)j * (m + 1))] = d_complex(rho, 0);
		h[(j + 1) + ((size_t)j * (m + 1))] = d_complex(0, 0);
		g[j + 1] = d_rmul(-1, d_mul(sn[j], g[j]));
		g[j] = d_mul(d_complex(cs[j].r, -cs[j].i), g[j]);

		// 残差の絶対値で判定
		if (d_abs(g[j + 1]) <= tolabs) {
			iters = j + 1;
			break;
		}
	}

	// 完了した反復数 k での最小二乗解 (未収束でもリスタート用に返す)
	const int k = (iters > 0) ? iters : j;
	memset(x, 0, (size_t)n * sizeof(d_complex_t));
	if (k > 0) {
		// 後退代入 y = H^-1 g
		for (int i = k - 1; i >= 0; i--) {
			d_complex_t s = g[i];
			for (int l = i + 1; l < k; l++) {
				s = d_sub(s, d_mul(h[i + ((size_t)l * (m + 1))], y[l]));
			}
			y[i] = d_div(s, h[i + ((size_t)i * (m + 1))]);
		}
		// u = V y、x = M^-1 u
		for (int l = 0; l < k; l++) {
			for (int i = 0; i < n; i++) {
				x[i] = d_add(x[i], d_mul(y[l], v[(size_t)l * n + i]));
			}
		}
		mop(mctx, x);
	}

	free(v); free(w); free(z); free(h); free(cs); free(sn); free(g); free(y);

	return iters;
}

// ── 密行列 + LU 前処理 (acceleration = 1) ─────────────────────────
typedef struct {
	int n;
	const d_complex_t *a;
} dense_a_t;

typedef struct {
	int n;
	const d_complex_t *alu;
	const int *piv;
} dense_m_t;

static void dense_aop(void *ctx, const d_complex_t *x, d_complex_t *y)
{
	const dense_a_t *c = (const dense_a_t *)ctx;
	zmatvec(c->n, c->a, x, y);
}

static void dense_mop(void *ctx, d_complex_t *z)
{
	const dense_m_t *c = (const dense_m_t *)ctx;
	lu_solve(c->n, c->alu, c->piv, z);
}

/*
右前処理 GMRES (1 サイクルのみ、従来の acceleration = 1 用)

  b : 入力 = 右辺 / 出力 = 解 (収束したときのみ上書き)

戻り値 : 反復回数 (>= 1) / -1 = 未収束 (b は変更しない。呼び出し側が
現在の周波数で LU を取り直す)
*/
int gmres_solve(int n, const d_complex_t *a, const d_complex_t *alu, const int *piv,
	d_complex_t *b, double tol)
{
	dense_a_t ac;
	dense_m_t mc;
	ac.n = n;
	ac.a = a;
	mc.n = n;
	mc.alu = alu;
	mc.piv = piv;

	d_complex_t *x = (d_complex_t *)malloc((size_t)n * sizeof(d_complex_t));
	if (x == NULL) return -1;
	const double beta = znorm(n, b);
	const int it = gmres_cycle(n, dense_aop, &ac, dense_mop, &mc,
		b, x, tol * beta);
	if (it > 0) {
		memcpy(b, x, (size_t)n * sizeof(d_complex_t));
	}
	free(x);

	return it;
}

// ── 行列フリー + 葉ブロック前処理 (compression = 1) ────────────────
typedef struct {
	const peec_t *p;
	const mna_sparse_t *sp;
	const hblock_t *blk;
	int nblk;
	d_complex_t *work;            // [max blk[].n]
} hmat_a_t;

typedef struct {
	const struct precond_t *pc;
} hmat_m_t;

static void hmat_aop(void *ctx, const d_complex_t *x, d_complex_t *y)
{
	const hmat_a_t *c = (const hmat_a_t *)ctx;
	mna_apply(c->p, c->sp, c->blk, c->nblk, x, y, c->work);
}

static void hmat_mop(void *ctx, d_complex_t *z)
{
	const hmat_m_t *c = (const hmat_m_t *)ctx;
	pc_apply(c->pc, z);
}

/*
行列フリー GMRES (compression = 1)。リスタートごとに真の残差 r = b - A x を
取り直すので、サイクル内の残差推定の丸めが蓄積しない。

戻り値 : 総反復数 (>= 1、収束時。b に解を上書き) / -1 = 未収束 (b は不定)
*/
int gmres_hmat(const peec_t *p, const mna_sparse_t *sp,
	const hblock_t *blk, int nblk,
	const struct precond_t *pc, d_complex_t *b, double tol)
{
	const int n = p->nunknown;
	hmat_a_t ac;
	hmat_m_t mc;
	int maxn = 1;
	for (int i = 0; i < nblk; i++) {
		if (blk[i].n > maxn) maxn = blk[i].n;
	}
	ac.p = p;
	ac.sp = sp;
	ac.blk = blk;
	ac.nblk = nblk;
	ac.work = (d_complex_t *)malloc((size_t)maxn * sizeof(d_complex_t));
	mc.pc = pc;

	d_complex_t *x = (d_complex_t *)calloc((size_t)n, sizeof(d_complex_t));
	d_complex_t *dx = (d_complex_t *)malloc((size_t)n * sizeof(d_complex_t));
	d_complex_t *r = (d_complex_t *)malloc((size_t)n * sizeof(d_complex_t));
	if ((ac.work == NULL) || (x == NULL) || (dx == NULL) || (r == NULL)) {
		free(ac.work); free(x); free(dx); free(r);
		return -1;
	}

	const double bnorm = znorm(n, b);
	if (bnorm <= 0) {
		// 右辺 0 : 解も 0
		memset(b, 0, (size_t)n * sizeof(d_complex_t));
		free(ac.work); free(x); free(dx); free(r);
		return 1;
	}
	const double tolabs = tol * bnorm;

	int total = 0;
	int done = 0;
	for (int cyc = 0; (cyc < GMRES_CYCLES) && !done; cyc++) {
		// 真の残差 r = b - A x (初回は x = 0 なので r = b)
		if (cyc == 0) {
			memcpy(r, b, (size_t)n * sizeof(d_complex_t));
		}
		else {
			hmat_aop(&ac, x, r);
			for (int i = 0; i < n; i++) {
				r[i] = d_sub(b[i], r[i]);
			}
			if (znorm(n, r) <= tolabs) {
				done = 1;
				break;
			}
		}
		const int it = gmres_cycle(n, hmat_aop, &ac, hmat_mop, &mc, r, dx, tolabs);
		for (int i = 0; i < n; i++) {
			x[i] = d_add(x[i], dx[i]);
		}
		total += (it > 0) ? it : GMRES_MAXIT;
		if (it > 0) done = 1;
	}

	if (done) {
		memcpy(b, x, (size_t)n * sizeof(d_complex_t));
	}
	free(ac.work);
	free(x);
	free(dx);
	free(r);

	return done ? total : -1;
}
