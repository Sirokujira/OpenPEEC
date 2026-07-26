/*
partial.c

部分インダクタンス (フィラメント近似・外部インダクタンス)
- 自己       : 平行フィラメント閉形式を距離 d = radius で評価
- 平行区間   : 解析式 (符号付き射影座標で向きも自動処理)
- 同一直線上 : d -> 0 の極限式 (特異回避)
- 一般       : Neumann 二重積分を 8 点 Gauss-Legendre x 2 で数値積分
- 遠方       : 中点近似
*/

#include "peec.h"

// 8 点 Gauss-Legendre ([-1,1])
static const double xg8[8] = {
	-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
	 0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
static const double wg8[8] = {
	 0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620};
// (対称なので重みは 4 個 : wg8[i] と wg8[7-i] が等しい)
static double weight8(int i)
{
	return (i < 4) ? wg8[i] : wg8[7 - i];
}

static double dot3(const double *a, const double *b)
{
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

// 平行フィラメントの原始関数 (偶関数)
static double fpar(double x, double d)
{
	return (x * asinh(x / d)) - sqrt((x * x) + (d * d));
}

// 同一直線上 (d -> 0) の原始関数 (偶関数、F0(0) = 0)
static double fcol(double x)
{
	const double ax = fabs(x);
	return (ax < 1e-300) ? 0 : ax * (log(ax) - 1);
}

// 自己部分インダクタンス (長さ l、半径 a)
double lp_self(double l, double a)
{
	const double u = l / a;
	return MU0 * l / (2 * PI) * (asinh(u) - sqrt(1 + (1 / (u * u))) + (1 / u));
}

// 2 区間の相互部分インダクタンス (符号付き : 各区間の n1 -> n2 向きを基準)
double lp_pair(const seg_t *s1, const seg_t *s2)
{
	double t1[3], t2[3];
	for (int c = 0; c < 3; c++) {
		t1[c] = (s1->x2[c] - s1->x1[c]) / s1->len;
		t2[c] = (s2->x2[c] - s2->x1[c]) / s2->len;
	}

	// 平行判定 : |t1 x t2|
	const double cx = (t1[1] * t2[2]) - (t1[2] * t2[1]);
	const double cy = (t1[2] * t2[0]) - (t1[0] * t2[2]);
	const double cz = (t1[0] * t2[1]) - (t1[1] * t2[0]);
	const double cross = sqrt((cx * cx) + (cy * cy) + (cz * cz));

	if (cross < 1e-9) {
		// s1 の軸に射影 (a1 = 0, a2 = len1)
		const double a1 = 0;
		const double a2 = s1->len;
		double r1[3], r2[3];
		for (int c = 0; c < 3; c++) {
			r1[c] = s2->x1[c] - s1->x1[c];
			r2[c] = s2->x2[c] - s1->x1[c];
		}
		const double b1 = dot3(r1, t1);
		const double b2 = dot3(r2, t1);
		double perp[3];
		for (int c = 0; c < 3; c++) {
			perp[c] = r1[c] - (b1 * t1[c]);
		}
		const double d = sqrt(dot3(perp, perp));
		const double lmax = (s1->len > s2->len) ? s1->len : s2->len;
		if (d < (1e-9 * lmax) + 1e-15) {
			// 同一直線上
			return MU0 / (4 * PI) *
				(fcol(a2 - b1) - fcol(a2 - b2) - fcol(a1 - b1) + fcol(a1 - b2));
		}
		return MU0 / (4 * PI) *
			(fpar(a2 - b1, d) - fpar(a2 - b2, d) - fpar(a1 - b1, d) + fpar(a1 - b2, d));
	}

	// 一般 (非平行)
	const double tdot = dot3(t1, t2);
	double c1[3], c2[3];
	for (int c = 0; c < 3; c++) {
		c1[c] = 0.5 * (s1->x1[c] + s1->x2[c]);
		c2[c] = 0.5 * (s2->x1[c] + s2->x2[c]);
	}
	const double dx = c2[0] - c1[0];
	const double dy = c2[1] - c1[1];
	const double dz = c2[2] - c1[2];
	const double rc = sqrt((dx * dx) + (dy * dy) + (dz * dz));
	const double lmax = (s1->len > s2->len) ? s1->len : s2->len;

	if (rc > 20 * lmax) {
		// 遠方 : 中点近似
		return MU0 / (4 * PI) * tdot * s1->len * s2->len / rc;
	}

	// Gauss-Legendre 8x8 (節点は内部にあるため端点共有でも特異にならない)
	double sum = 0;
	for (int i = 0; i < 8; i++) {
		const double u = 0.5 * (1 + xg8[i]);      // [0,1]
		double q1[3];
		for (int c = 0; c < 3; c++) {
			q1[c] = s1->x1[c] + ((s1->x2[c] - s1->x1[c]) * u);
		}
		for (int j = 0; j < 8; j++) {
			const double v = 0.5 * (1 + xg8[j]);
			double rx = s2->x1[0] + ((s2->x2[0] - s2->x1[0]) * v) - q1[0];
			double ry = s2->x1[1] + ((s2->x2[1] - s2->x1[1]) * v) - q1[1];
			double rz = s2->x1[2] + ((s2->x2[2] - s2->x1[2]) * v) - q1[2];
			const double r = sqrt((rx * rx) + (ry * ry) + (rz * rz));
			if (r > 1e-300) {
				sum += weight8(i) * weight8(j) / r;
			}
		}
	}
	// jacobian : (len1/2)(len2/2)
	return MU0 / (4 * PI) * tdot * s1->len * s2->len * 0.25 * sum;
}

// 部分インダクタンス行列 (対称密行列、周波数非依存なので 1 回だけ計算)
void lp_fill(peec_t *p, FILE *fp_log)
{
	const int n = p->nseg;
	if (n <= 0) return;

	p->lp = (double *)malloc((size_t)n * n * sizeof(double));

	// MSVC の OpenMP 2.0 は for 文内でのインデックス宣言を許さない (C3015) ため
	// ループ変数は事前に宣言する
	int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (i = 0; i < n; i++) {
		p->lp[(size_t)i * n + i] = lp_self(p->seg[i].len, p->seg[i].radius);
		for (int j = i + 1; j < n; j++) {
			p->lp[(size_t)i * n + j] = lp_pair(&p->seg[i], &p->seg[j]);
		}
	}
	// 下三角へミラー
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			p->lp[(size_t)i * n + j] = p->lp[(size_t)j * n + i];
		}
	}

	// 非物理的結合の検出 (|k| > 1 : 導体の重なり等)
	int nwarn = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			const double m = fabs(p->lp[(size_t)i * n + j]);
			const double s = sqrt(p->lp[(size_t)i * n + i] * p->lp[(size_t)j * n + j]);
			if ((m > s) && (nwarn < 10)) {
				fprintf(fp_log, "*** warning : |k| > 1 between segment %d and %d (overlapping conductors?)\n", i, j);
				nwarn++;
			}
		}
	}
}
