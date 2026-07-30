/*
lu.c

複素 LU 分解 (Doolittle, 部分ピボット)。フラット配列 a[i*n+j] を in-place で分解する。
外部ライブラリ (LAPACK 等) には依存しない。
*/

#include "peec.h"

static double mag(d_complex_t z)
{
	return fabs(z.r) + fabs(z.i);
}

// 戻り値 : -1 = 正常、>= 0 = 特異 (ゼロピボットの行番号)
int lu_decomp(int n, d_complex_t *a, int *piv)
{
	for (int k = 0; k < n; k++) {
		// pivot 探索 (列 k の最大絶対値)
		int    ip = k;
		double amax = mag(a[(size_t)k * n + k]);
		for (int i = k + 1; i < n; i++) {
			const double m = mag(a[(size_t)i * n + k]);
			if (m > amax) {
				amax = m;
				ip = i;
			}
		}
		piv[k] = ip;
		if (amax < 1e-30) return k;

		// 行交換
		if (ip != k) {
			for (int j = 0; j < n; j++) {
				const d_complex_t t = a[(size_t)k * n + j];
				a[(size_t)k * n + j] = a[(size_t)ip * n + j];
				a[(size_t)ip * n + j] = t;
			}
		}

		// 消去 (残余部分行列の更新)
		//
		// 行 i の更新は行 k しか読まないので行ごとに独立 = 並列化できる。
		// 各要素は「元の値 - lik * a[k][j]」を 1 回計算するだけで、順序に依存
		// する加算 (リダクション) が無いため、**スレッド数を変えても結果は
		// ビット単位で同一**になる (peec_check.sh が -n 1 と -n 4 の一致を判定)。
		//
		// 小さい行列ではスレッド生成のほうが高くつくので if 節で切り替える。
		// MSVC の OpenMP 2.0 は for 文内でのインデックス宣言を許さない (C3015)
		// ため、ループ変数は事前に宣言する。
		const d_complex_t pinv = d_inv(a[(size_t)k * n + k]);
		int i;
#ifdef _OPENMP
#pragma omp parallel for if ((n - k) > 64)
#endif
		for (i = k + 1; i < n; i++) {
			const d_complex_t lik = d_mul(a[(size_t)i * n + k], pinv);
			a[(size_t)i * n + k] = lik;
			for (int j = k + 1; j < n; j++) {
				a[(size_t)i * n + j] = d_sub(a[(size_t)i * n + j], d_mul(lik, a[(size_t)k * n + j]));
			}
		}
	}

	return -1;
}

void lu_solve(int n, const d_complex_t *a, const int *piv, d_complex_t *b)
{
	// pivot 適用
	for (int k = 0; k < n; k++) {
		if (piv[k] != k) {
			const d_complex_t t = b[k];
			b[k] = b[piv[k]];
			b[piv[k]] = t;
		}
	}

	// 前進代入 (L : 単位下三角)
	for (int i = 1; i < n; i++) {
		d_complex_t s = b[i];
		for (int j = 0; j < i; j++) {
			s = d_sub(s, d_mul(a[(size_t)i * n + j], b[j]));
		}
		b[i] = s;
	}

	// 後退代入
	for (int i = n - 1; i >= 0; i--) {
		d_complex_t s = b[i];
		for (int j = i + 1; j < n; j++) {
			s = d_sub(s, d_mul(a[(size_t)i * n + j], b[j]));
		}
		b[i] = d_mul(s, d_inv(a[(size_t)i * n + i]));
	}
}
