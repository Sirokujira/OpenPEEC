/*
partial.c

幾何二重積分 (Neumann 積分) と部分インダクタンス

  I(s1,s2) = ∫∫ dl1 dl2 / R          … 向きに依らない正の幾何量
  Lp(s1,s2) = (mu0/4pi) (t1・t2) I    … 部分インダクタンス
  P(c1,c2)  = I / (4 pi eps0 L1 L2)   … 電位係数 (potential.c で使用)

I の評価 :
- 自己       : 平行フィラメント閉形式を距離 d = radius で評価
- 平行区間   : 解析式
- 同一直線上 : d -> 0 の極限式 (特異回避)
- 一般       : 8 点 Gauss-Legendre の複合則 (近接時は 4x4 分割)
- 遠方       : 中点近似
*/

#include "peec.h"

// 8 点 Gauss-Legendre ([-1,1])
static const double xg8[8] = {
	-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
	 0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
static const double wg8[4] = {
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

// 自己二重積分 : ∫∫ dl dl' / sqrt((l-l')^2 + a^2)
double neumann_self(double l, double a)
{
	return 2 * ((l * asinh(l / a)) - sqrt((l * l) + (a * a)) + a);
}

// 一般配置の数値積分 (各区間を nsub 分割した複合 8 点 Gauss-Legendre)
static double gauss_pair(const seg_t *s1, const seg_t *s2, int nsub)
{
	double sum = 0;

	for (int i1 = 0; i1 < nsub; i1++) {
	for (int i2 = 0; i2 < nsub; i2++) {
		for (int i = 0; i < 8; i++) {
			const double u = (i1 + (0.5 * (1 + xg8[i]))) / nsub;
			double q1[3];
			for (int c = 0; c < 3; c++) {
				q1[c] = s1->x1[c] + ((s1->x2[c] - s1->x1[c]) * u);
			}
			for (int j = 0; j < 8; j++) {
				const double v = (i2 + (0.5 * (1 + xg8[j]))) / nsub;
				const double rx = s2->x1[0] + ((s2->x2[0] - s2->x1[0]) * v) - q1[0];
				const double ry = s2->x1[1] + ((s2->x2[1] - s2->x1[1]) * v) - q1[1];
				const double rz = s2->x1[2] + ((s2->x2[2] - s2->x1[2]) * v) - q1[2];
				const double r = sqrt((rx * rx) + (ry * ry) + (rz * rz));
				if (r > 1e-300) {
					sum += weight8(i) * weight8(j) / r;
				}
			}
		}
	}
	}

	// jacobian : (len1/2/nsub)(len2/2/nsub) を nsub^2 個の小領域について合計
	return s1->len * s2->len * 0.25 * sum / ((double)nsub * nsub);
}

// 2 区間の幾何二重積分 (常に正 : 向きの符号は呼び出し側で扱う)
double neumann_pair(const seg_t *s1, const seg_t *s2)
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
	const double lmax = (s1->len > s2->len) ? s1->len : s2->len;

	if (cross < 1e-9) {
		// s1 の軸に射影 (a1 = 0, a2 = len1)
		const double a1 = 0;
		const double a2 = s1->len;
		double r1[3], r2[3];
		for (int c = 0; c < 3; c++) {
			r1[c] = s2->x1[c] - s1->x1[c];
			r2[c] = s2->x2[c] - s1->x1[c];
		}
		double b1 = dot3(r1, t1);
		double b2 = dot3(r2, t1);
		if (b1 > b2) {
			// 逆向きの区間 : 積分区間を昇順にして正値を返す
			const double t = b1;
			b1 = b2;
			b2 = t;
		}
		double perp[3];
		for (int c = 0; c < 3; c++) {
			perp[c] = r1[c] - (dot3(r1, t1) * t1[c]);
		}
		// 細線カーネル : 電流は軸上、電位・磁束は導体表面で評価する。
		// 距離を半径で下限打ち切りすることで自己項と同じカーネルになり、
		// 同一直線上 (d -> 0) の特異性も消える。分割数によらず
		// sum_ij I_ij = I_self(全長) が厳密に成立する。
		const double aeff = 0.5 * (s1->radius + s2->radius);
		double d = sqrt(dot3(perp, perp));
		if (d < aeff) d = aeff;

		return fpar(a2 - b1, d) - fpar(a2 - b2, d) - fpar(a1 - b1, d) + fpar(a1 - b2, d);
	}

	// 一般 (非平行)
	double c1[3], c2[3];
	for (int c = 0; c < 3; c++) {
		c1[c] = 0.5 * (s1->x1[c] + s1->x2[c]);
		c2[c] = 0.5 * (s2->x1[c] + s2->x2[c]);
	}
	const double dx = c2[0] - c1[0];
	const double dy = c2[1] - c1[1];
	const double dz = c2[2] - c1[2];
	const double rc = sqrt((dx * dx) + (dy * dy) + (dz * dz));

	// 遠方 : 中点近似
	if (rc > 20 * lmax) {
		return s1->len * s2->len / rc;
	}

	// 近接 (端点共有を含む) は分割数を上げる。Gauss 節点は区間内部にあるため
	// 端点を共有していても被積分関数は有限にとどまる。
	return gauss_pair(s1, s2, (rc > 2 * lmax) ? 1 : 4);
}

// 自己部分インダクタンス (長さ l、半径 a)
double lp_self(double l, double a)
{
	return MU0 / (4 * PI) * neumann_self(l, a);
}

// 2 区間の相互部分インダクタンス (符号付き : 各区間の n1 -> n2 向きを基準)
double lp_pair(const seg_t *s1, const seg_t *s2)
{
	double t1[3], t2[3];
	for (int c = 0; c < 3; c++) {
		t1[c] = (s1->x2[c] - s1->x1[c]) / s1->len;
		t2[c] = (s2->x2[c] - s2->x1[c]) / s2->len;
	}

	return MU0 / (4 * PI) * dot3(t1, t2) * neumann_pair(s1, s2);
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
	for (i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			p->lp[(size_t)i * n + j] = p->lp[(size_t)j * n + i];
		}
	}

	// 非物理的結合の検出 (|k| > 1 : 導体の重なり等)
	int nwarn = 0;
	for (i = 0; i < n; i++) {
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
