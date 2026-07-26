/*
partial.c

幾何二重積分 (Neumann 積分) と部分インダクタンス

  I(s1,s2) = ∬ exp(-j k R) / R dl1 dl2,  R = sqrt(|dr|^2 + a^2)  (細線縮約カーネル)
  Lp(s1,s2) = (mu0/4pi) (t1・t2) I     … 部分インダクタンス
  P(c1,c2)  = I / (4 pi eps0 L1 L2)    … 電位係数 (potential.c で使用)

k = 0 (準静的) では I は実数で、次の解析式・数値積分で評価する :
- 自己       : 平行フィラメント閉形式 (距離 = 等価半径 a)
- 平行区間   : 解析式 (距離を sqrt(perp^2 + a^2) に正則化)
- 一般       : 8 点 Gauss-Legendre の複合則 (近接時は 4x4 分割)
- 遠方       : 中点近似

k > 0 (遅延あり) では特異性抽出を使う :
  I = I_static + ∬ (exp(-j k R) - 1) / R
第 2 項は R -> 0 で -jk に収束する正則な被積分関数なので、解析式を保ったまま
Gauss 求積で安全に評価できる。

全ての項で同一のカーネル (距離を等価半径 a で正則化) を使うことにより
sum_ij I_ij = I_self(全長) が厳密に成立する (CLAUDE.md の不変条件)。
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

// 自己二重積分 : ∬ dl dl' / sqrt((l-l')^2 + a^2)
double neumann_self(double l, double a)
{
	return 2 * ((l * asinh(l / a)) - sqrt((l * l) + (a * a)) + a);
}

// 数値積分の共通ループ。mode = 0 : 静的核 1/R、mode = 1 : 遅延補正 (exp(-jkR)-1)/R
static d_complex_t gauss_kernel(const seg_t *s1, const seg_t *s2,
	double aeff, double kw, int nsub, int mode)
{
	const double a2 = aeff * aeff;
	double sumr = 0, sumi = 0;

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
				const double dx = s2->x1[0] + ((s2->x2[0] - s2->x1[0]) * v) - q1[0];
				const double dy = s2->x1[1] + ((s2->x2[1] - s2->x1[1]) * v) - q1[1];
				const double dz = s2->x1[2] + ((s2->x2[2] - s2->x1[2]) * v) - q1[2];
				const double r = sqrt((dx * dx) + (dy * dy) + (dz * dz) + a2);
				const double w = weight8(i) * weight8(j);
				if (mode == 0) {
					sumr += w / r;
				}
				else {
					// (cos(kR) - 1 - j sin(kR)) / R : R -> 0 で -jk に収束する
					sumr += w * (cos(kw * r) - 1) / r;
					sumi -= w * sin(kw * r) / r;
				}
			}
		}
	}
	}

	// jacobian : (len1/2/nsub)(len2/2/nsub) を nsub^2 個の小領域について合計
	const double fac = s1->len * s2->len * 0.25 / ((double)nsub * nsub);

	return d_complex(fac * sumr, fac * sumi);
}

// 静的な幾何二重積分 (常に正)。a1, a2 は各区間の等価半径。
double neumann_pair(const seg_t *s1, const seg_t *s2, double a1, double a2)
{
	const double aeff = 0.5 * (a1 + a2);
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
		const double p1 = 0;
		const double p2 = s1->len;
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
		// 縮約カーネル : R = sqrt(軸方向^2 + perp^2 + a^2)
		const double d = sqrt(dot3(perp, perp) + (aeff * aeff));

		return fpar(p2 - b1, d) - fpar(p2 - b2, d) - fpar(p1 - b1, d) + fpar(p1 - b2, d);
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

	// 近接 (端点共有を含む) は分割数を上げる
	return gauss_kernel(s1, s2, aeff, 0, (rc > 2 * lmax) ? 1 : 4, 0).r;
}

// 遅延補正の分割数 : 位相変化 k*L と近接度から決める
static int corr_nsub(const seg_t *s1, const seg_t *s2, double kw)
{
	const double lmax = (s1->len > s2->len) ? s1->len : s2->len;
	int nsub = 2 + (int)(kw * lmax);
	if (nsub > 8) nsub = 8;

	return nsub;
}

// 平行 (= 解析式が厳密に使える) かどうか
static int is_parallel(const seg_t *s1, const seg_t *s2)
{
	double t1[3], t2[3];
	for (int c = 0; c < 3; c++) {
		t1[c] = (s1->x2[c] - s1->x1[c]) / s1->len;
		t2[c] = (s2->x2[c] - s2->x1[c]) / s2->len;
	}
	const double cx = (t1[1] * t2[2]) - (t1[2] * t2[1]);
	const double cy = (t1[2] * t2[0]) - (t1[0] * t2[2]);
	const double cz = (t1[0] * t2[1]) - (t1[1] * t2[0]);

	return (sqrt((cx * cx) + (cy * cy) + (cz * cz)) < 1e-9);
}

// 面導体セル (リボン) が絡むか
static int is_ribbon(const seg_t *s1, const seg_t *s2)
{
	return ((s1->wid > 0) || (s2->wid > 0));
}

// リボン用の分割数 : 近接ほど、また位相変化が大きいほど細かくする
static int ribbon_nsub(const seg_t *s1, const seg_t *s2, double kw)
{
	double c1[3], c2[3];
	for (int i = 0; i < 3; i++) {
		c1[i] = 0.5 * (s1->x1[i] + s1->x2[i]);
		c2[i] = 0.5 * (s2->x1[i] + s2->x2[i]);
	}
	const double dx = c2[0] - c1[0];
	const double dy = c2[1] - c1[1];
	const double dz = c2[2] - c1[2];
	const double rc = sqrt((dx * dx) + (dy * dy) + (dz * dz));
	double lmax = (s1->len > s2->len) ? s1->len : s2->len;
	if (s1->wid > lmax) lmax = s1->wid;
	if (s2->wid > lmax) lmax = s2->wid;

	int nsub = (rc > 2 * lmax) ? 1 : 2;
	nsub += (int)(kw * lmax);
	if (nsub > 4) nsub = 4;

	return nsub;
}

// 遅延を含む幾何二重積分 (kw = 0 なら静的値がそのまま返る)
//
// 放射抵抗は L 項と P 項の大きな値どうしの差として現れるため、静的部と
// 遅延補正部の評価法が食い違うとその差が汚染される。解析式が厳密な
// 平行区間以外では、静的部も補正部と同じ求積で評価する。
d_complex_t neumann_pair_k(const seg_t *s1, const seg_t *s2, double a1, double a2, double kw)
{
	if (is_ribbon(s1, s2)) {
		const int nsub = ribbon_nsub(s1, s2, kw);
		const d_complex_t is = d_complex(ribbon_static(s1, s2, nsub), 0);
		if (kw <= 0) return is;
		return d_add(is, ribbon_corr(s1, s2, kw, nsub));
	}

	if (kw <= 0) return d_complex(neumann_pair(s1, s2, a1, a2), 0);

	const double aeff = 0.5 * (a1 + a2);
	const int nsub = corr_nsub(s1, s2, kw);

	const double is = is_parallel(s1, s2)
		? neumann_pair(s1, s2, a1, a2)
		: gauss_kernel(s1, s2, aeff, 0, nsub, 0).r;

	return d_add(d_complex(is, 0), gauss_kernel(s1, s2, aeff, kw, nsub, 1));
}

// 遅延を含む自己二重積分
d_complex_t neumann_self_k(const seg_t *s, double a, double kw)
{
	if (s->wid > 0) return neumann_pair_k(s, s, a, a, kw);

	const d_complex_t is = d_complex(neumann_self(s->len, a), 0);
	if (kw <= 0) return is;

	return d_add(is, gauss_kernel(s, s, a, kw, corr_nsub(s, s, kw), 1));
}

// 自己部分インダクタンス (長さ l、等価半径 a)
double lp_self(double l, double a)
{
	return MU0 / (4 * PI) * neumann_self(l, a);
}

// 部分インダクタンス行列 (対称密行列)
// retardation = 0 なら周波数非依存なので 1 回だけ、1 なら周波数ごとに呼ぶ
void lp_fill(peec_t *p, double f, FILE *fp_log)
{
	const int n = p->nseg;
	if (n <= 0) return;

	if (p->lp == NULL) {
		p->lp = (d_complex_t *)malloc((size_t)n * n * sizeof(d_complex_t));
	}
	const double kw = p->retardation ? (2 * PI * f / C0) : 0;
	const double coef = MU0 / (4 * PI);

	// MSVC の OpenMP 2.0 は for 文内でのインデックス宣言を許さない (C3015) ため
	// ループ変数は事前に宣言する
	int i;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (i = 0; i < n; i++) {
		p->lp[(size_t)i * n + i] = d_rmul(coef, neumann_self_k(&p->seg[i], p->seg[i].aL, kw));
		for (int j = i + 1; j < n; j++) {
			double t1[3], t2[3];
			for (int c = 0; c < 3; c++) {
				t1[c] = (p->seg[i].x2[c] - p->seg[i].x1[c]) / p->seg[i].len;
				t2[c] = (p->seg[j].x2[c] - p->seg[j].x1[c]) / p->seg[j].len;
			}
			const double tdot = dot3(t1, t2);
			p->lp[(size_t)i * n + j] = d_rmul(coef * tdot,
				neumann_pair_k(&p->seg[i], &p->seg[j], p->seg[i].aL, p->seg[j].aL, kw));
		}
	}
	// 下三角へミラー
	for (i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			p->lp[(size_t)i * n + j] = p->lp[(size_t)j * n + i];
		}
	}

	// 非物理的結合の検出 (|k| > 1 : 導体の重なり等) は静的な初回のみ報告する
	if (p->lpwarn) return;
	p->lpwarn = 1;
	int nwarn = 0;
	for (i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			const double m = fabs(p->lp[(size_t)i * n + j].r);
			const double s = sqrt(fabs(p->lp[(size_t)i * n + i].r * p->lp[(size_t)j * n + j].r));
			if ((m > s) && (nwarn < 10)) {
				fprintf(fp_log, "*** warning : |k| > 1 between segment %d and %d (overlapping conductors?)\n", i, j);
				nwarn++;
			}
		}
	}
}
