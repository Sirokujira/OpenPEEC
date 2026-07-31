/*
iterative.c

周波数掃引の反復解法 (acceleration = 1 のときのみ使用)

掃引では隣接周波数の MNA 行列 A(f) = S + j omega T が近いことを利用して、
直近に LU 分解した行列 M = A(f_LU) を右前処理に使う GMRES で解く :

  A M^-1 u = b,   x = M^-1 u

収束しない (= 前処理が古い) ときは呼び出し側 (solve.c) が現在の周波数で
LU を取り直す。最悪でも従来の毎周波数 LU と同等で、掃引点が密なほど
LU の回数が減り、O(n^3) が O(反復 x n^2) に置き換わる。

スレッド数不変性 : 行列ベクトル積は行ごとに独立 (各行の内積は直列)、
Gram-Schmidt / ノルム / Givens 回転はすべて直列なので、結果はスレッド数に
依らずビット単位で一致する。順序依存の並列リダクションは入れないこと。
*/

#include "peec.h"

#define GMRES_MAXIT 60

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
右前処理 GMRES (再スタートなし、最大 GMRES_MAXIT 反復)

  a    : 現在の周波数の MNA 行列 (n x n、未分解のまま)
  alu  : 前処理の LU 因子 (lu_decomp 済み) と piv
  b    : 入力 = 右辺 / 出力 = 解 (収束したときのみ上書き)
  tol  : 相対残差の収束判定

戻り値 : 反復回数 (>= 1) / -1 = 未収束 (b は変更しない)
*/
int gmres_solve(int n, const d_complex_t *a, const d_complex_t *alu, const int *piv,
	d_complex_t *b, double tol)
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

	// 初期残差 = b (x0 = 0)
	const double beta = znorm(n, b);
	int iters = -1;
	if (beta <= 0) {
		// 右辺 0 : 解も 0
		free(v); free(w); free(z); free(h); free(cs); free(sn); free(g); free(y);
		return 1;
	}
	for (int i = 0; i < n; i++) {
		v[i] = d_rmul(1 / beta, b[i]);
	}
	g[0] = d_complex(beta, 0);

	int j;
	for (j = 0; j < m; j++) {
		// w = A M^-1 v_j
		memcpy(z, &v[(size_t)j * n], (size_t)n * sizeof(d_complex_t));
		lu_solve(n, alu, piv, z);
		zmatvec(n, a, z, w);

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

		// 相対残差
		if (d_abs(g[j + 1]) <= tol * beta) {
			iters = j + 1;
			break;
		}
	}

	if (iters > 0) {
		// 後退代入 y = H^-1 g
		const int k = iters;
		for (int i = k - 1; i >= 0; i--) {
			d_complex_t s = g[i];
			for (int l = i + 1; l < k; l++) {
				s = d_sub(s, d_mul(h[i + ((size_t)l * (m + 1))], y[l]));
			}
			y[i] = d_div(s, h[i + ((size_t)i * (m + 1))]);
		}
		// u = V y、x = M^-1 u
		memset(z, 0, (size_t)n * sizeof(d_complex_t));
		for (int l = 0; l < k; l++) {
			for (int i = 0; i < n; i++) {
				z[i] = d_add(z[i], d_mul(y[l], v[(size_t)l * n + i]));
			}
		}
		lu_solve(n, alu, piv, z);
		memcpy(b, z, (size_t)n * sizeof(d_complex_t));
	}

	free(v); free(w); free(z); free(h); free(cs); free(sn); free(g); free(y);

	return iters;
}
