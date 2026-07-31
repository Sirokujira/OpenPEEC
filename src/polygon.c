/*
polygon.c

平面多角形セル (三角形メッシュ由来) の幾何二重積分

一様面密度の平面多角形がつくるポテンシャルの閉形式 (辺ごとの和) :

  Wilton, Rao, Glisson, Schaubert, Al-Bundak, Butler,
  "Potential integrals for uniform and linear source distributions on
   polygonal and polyhedral domains", IEEE Trans. AP-32, pp.276-281, 1984.

  ∬_S dS'/R = Σ_edges [ t0 ln((R+ + l+)/(R- + l-))
                        - |z| (atan(t0 l+/(R0² + |z| R+))
                             - atan(t0 l-/(R0² + |z| R-))) ]

リボン (surface.c) と同じ構造で、内側 (解析) x 外側 (三角形 Gauss 求積)
の 2 段で対の積分を評価し、幅で規格化した Î を返す :

  Î = ∬∬ dS dS'/R / (w1 w2),   w = 面積/len (平均幅)

これにより Lp = (mu0/4pi)(t1・t2) Î、P = Î/(4 pi eps0 L1 L2) という細線・
リボンと同一の式がそのまま使える (physics-invariants (4) の一般化 :
len・wid = 面積 なので P の規格化は 1/(4 pi eps0 A1 A2) に一致する)。

セルは papex から星形 (papex + 頂点リング pv の隣接対で扇分割できる)。
papex がリング上の頂点と一致する場合は退化三角形 (面積 0) ができるが
重みが 0 なので無害。多角形の向き (巻き方) には依存しない (閉形式の和は
向きで全体符号が反転するだけなので絶対値を取る)。
*/

#include "peec.h"

static double dot3(const double *a, const double *b)
{
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

static void cross3(const double *a, const double *b, double *c)
{
	c[0] = (a[1] * b[2]) - (a[2] * b[1]);
	c[1] = (a[2] * b[0]) - (a[0] * b[2]);
	c[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

// Newell の公式による面法線ベクトル (大きさ = 2 x 面積)
static void poly_normal2(const seg_t *s, double *nv)
{
	nv[0] = nv[1] = nv[2] = 0;
	for (int k = 0; k < s->npv; k++) {
		const double *a = &s->pv[3 * k];
		const double *b = &s->pv[3 * ((k + 1) % s->npv)];
		double c[3];
		cross3(a, b, c);
		for (int i = 0; i < 3; i++) {
			nv[i] += c[i];
		}
	}
}

double poly_area(const seg_t *s)
{
	double nv[3];
	poly_normal2(s, nv);
	return 0.5 * sqrt(dot3(nv, nv));
}

// 点 pt から見た一様多角形セル s 上の ∬ dS'/R (解析式、常に正)
double poly_potential(const seg_t *s, const double *pt)
{
	double nv[3];
	poly_normal2(s, nv);
	const double n2 = sqrt(dot3(nv, nv));
	if (n2 <= 0) return 0;
	const double nh[3] = {nv[0] / n2, nv[1] / n2, nv[2] / n2};

	double d0[3];
	for (int i = 0; i < 3; i++) {
		d0[i] = pt[i] - s->pv[i];
	}
	const double z = dot3(d0, nh);
	const double az = fabs(z);
	// 面内投影点
	double rho[3];
	for (int i = 0; i < 3; i++) {
		rho[i] = pt[i] - (z * nh[i]);
	}

	double sum = 0;
	for (int k = 0; k < s->npv; k++) {
		const double *a = &s->pv[3 * k];
		const double *b = &s->pv[3 * ((k + 1) % s->npv)];
		double e[3];
		for (int i = 0; i < 3; i++) {
			e[i] = b[i] - a[i];
		}
		const double el = sqrt(dot3(e, e));
		if (el <= 0) continue;
		const double sh[3] = {e[0] / el, e[1] / el, e[2] / el};
		double mh[3];
		cross3(sh, nh, mh);              // 面内の辺法線 (向きは巻き方に従う)

		double ra[3], rb[3];
		for (int i = 0; i < 3; i++) {
			ra[i] = a[i] - rho[i];
			rb[i] = b[i] - rho[i];
		}
		const double lm = dot3(ra, sh);
		const double lp = dot3(rb, sh);
		const double t0 = dot3(ra, mh);  // 辺までの面内符号付き距離
		const double r02 = (t0 * t0) + (z * z);
		const double rm = sqrt((lm * lm) + r02);
		const double rp = sqrt((lp * lp) + r02);

		// 点が辺の線上 (r02 = 0) : t0 項も |z| 項も極限で 0
		if (r02 <= 0) continue;

		// ln((R+ + l+)/(R- + l-)) : l < 0 の側は R + l が桁落ちするので
		// 恒等式 R + l = R0²/(R - l) で計算する ((R+l)(R-l) = R0²)。
		// 隣接セルの求積点が辺の延長線近傍に乗っても 0/0 にならない。
		const double sp = (lp >= 0) ? (rp + lp) : (r02 / (rp - lp));
		const double sm = (lm >= 0) ? (rm + lm) : (r02 / (rm - lm));
		sum += t0 * log(sp / sm);
		if (az > 0) {
			sum -= az * (atan((t0 * lp) / (r02 + (az * rp)))
			           - atan((t0 * lm) / (r02 + (az * rm))));
		}
	}

	// 巻き方で全体符号が反転するだけなので絶対値 (積分値は常に正)
	return fabs(sum);
}

// 7 点三角形 Gauss 則 (5 次精度、重みは面積和 = 1 に規格化済み)
static const double tw7[7] = {
	0.225,
	0.1323941527885062, 0.1323941527885062, 0.1323941527885062,
	0.1259391805448271, 0.1259391805448271, 0.1259391805448271};
static const double tb7[7][3] = {
	{1.0 / 3, 1.0 / 3, 1.0 / 3},
	{0.0597158717897698, 0.4701420641051151, 0.4701420641051151},
	{0.4701420641051151, 0.0597158717897698, 0.4701420641051151},
	{0.4701420641051151, 0.4701420641051151, 0.0597158717897698},
	{0.7974269853530873, 0.1012865073234563, 0.1012865073234563},
	{0.1012865073234563, 0.7974269853530873, 0.1012865073234563},
	{0.1012865073234563, 0.1012865073234563, 0.7974269853530873}};

// 3 点三角形 Gauss 則 (2 次精度、遅延補正の滑らかな被積分関数用)
static const double tw3[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
static const double tb3[3][3] = {
	{2.0 / 3, 1.0 / 6, 1.0 / 6},
	{1.0 / 6, 2.0 / 3, 1.0 / 6},
	{1.0 / 6, 1.0 / 6, 2.0 / 3}};

static double tri_area(const double *a, const double *b, const double *c)
{
	double u[3], v[3], w[3];
	for (int i = 0; i < 3; i++) {
		u[i] = b[i] - a[i];
		v[i] = c[i] - a[i];
	}
	cross3(u, v, w);
	return 0.5 * sqrt(dot3(w, w));
}

// セルの求積点列を作る (最大数 : POLY_MAX 三角形 x 4 分割 x 7 点)
#define POLY_NQMAX (POLY_MAX * 4 * 7)

// 多角形セルの求積点 : papex からの扇 x (nsub = 2 なら各三角形を 4 分割) x n 点則
// 戻り値は点数。wt には点の重み x 微小面積が入る (合計 = セル面積)。
static int poly_qpts(const seg_t *s, int nsub, int n7,
	double *px, double *wt)
{
	const int nq = n7 ? 7 : 3;
	const double (*tb)[3] = n7 ? tb7 : tb3;
	const double *tw = n7 ? tw7 : tw3;
	int np = 0;

	for (int k = 0; k < s->npv; k++) {
		const double *b = &s->pv[3 * k];
		const double *c = &s->pv[3 * ((k + 1) % s->npv)];
		const double *a = s->papex;
		double v0[3], v1[3], v2[3];
		const double at = tri_area(a, b, c);
		if (at <= 0) continue;      // papex がリング上にある場合の退化扇
		// nsub = 1 : そのまま、nsub = 2 : 4 つの相似三角形に分割
		const int ns = (nsub >= 2) ? 4 : 1;
		for (int t = 0; t < ns; t++) {
			// 部分三角形の頂点 (重心座標で表現)
			static const double sub[4][3][3] = {
				{{1, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}},
				{{0.5, 0.5, 0}, {0, 1, 0}, {0, 0.5, 0.5}},
				{{0.5, 0, 0.5}, {0, 0.5, 0.5}, {0, 0, 1}},
				{{0.5, 0.5, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}}};
			const double (*sb)[3] = (ns == 4) ? sub[t] : NULL;
			for (int i = 0; i < 3; i++) {
				double bc0[3] = {1, 0, 0}, bc1[3] = {0, 1, 0}, bc2[3] = {0, 0, 1};
				if (sb != NULL) {
					for (int j = 0; j < 3; j++) {
						bc0[j] = sb[0][j];
						bc1[j] = sb[1][j];
						bc2[j] = sb[2][j];
					}
				}
				v0[i] = (bc0[0] * a[i]) + (bc0[1] * b[i]) + (bc0[2] * c[i]);
				v1[i] = (bc1[0] * a[i]) + (bc1[1] * b[i]) + (bc1[2] * c[i]);
				v2[i] = (bc2[0] * a[i]) + (bc2[1] * b[i]) + (bc2[2] * c[i]);
			}
			const double sat = at / ns;
			for (int q = 0; q < nq; q++) {
				for (int i = 0; i < 3; i++) {
					px[3 * np + i] = (tb[q][0] * v0[i]) + (tb[q][1] * v1[i])
					               + (tb[q][2] * v2[i]);
				}
				wt[np] = tw[q] * sat;
				np++;
			}
		}
	}

	return np;
}

// 細線 (wid = 0) / リボン / 多角形 いずれかのセルの求積点列
static int cell_qpts(const seg_t *s, int nsub, int n7, double *px, double *wt)
{
	if (s->npv > 0) return poly_qpts(s, nsub, n7, px, wt);

	// 細線 : 4 点 Gauss-Legendre x nsub 分割 (重み合計 = len)
	static const double xg4[4] = {
		-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
	static const double wg4[4] = {
		 0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
	int np = 0;
	if (s->wid <= 0) {
		for (int is = 0; is < nsub; is++) {
			for (int q = 0; q < 4; q++) {
				const double u = (is + (0.5 * (1 + xg4[q]))) / nsub;
				for (int i = 0; i < 3; i++) {
					px[3 * np + i] = s->x1[i] + ((s->x2[i] - s->x1[i]) * u);
				}
				wt[np] = 0.5 * wg4[q] * s->len / nsub;
				np++;
			}
		}
		return np;
	}

	// リボン : 軸 4 点 x nsub、横 4 点 (重み合計 = len x wid)
	double t[3], c[3];
	for (int i = 0; i < 3; i++) {
		t[i] = (s->x2[i] - s->x1[i]) / s->len;
		c[i] = 0.5 * (s->x1[i] + s->x2[i]);
	}
	for (int is = 0; is < nsub; is++) {
		for (int q = 0; q < 4; q++) {
			const double a = s->len * (((is + (0.5 * (1 + xg4[q]))) / nsub) - 0.5);
			for (int j = 0; j < 4; j++) {
				const double b = s->wid * (0.5 * xg4[j]);
				for (int i = 0; i < 3; i++) {
					px[3 * np + i] = c[i] + (a * t[i]) + (b * s->wv[i]);
				}
				wt[np] = 0.25 * wg4[q] * wg4[j] * s->len * s->wid / nsub;
				np++;
			}
		}
	}
	return np;
}

// 内側の解析ポテンシャル (多角形 / 矩形リボンを振り分け)
static double cell_potential(const seg_t *s, const double *pt)
{
	return (s->npv > 0) ? poly_potential(s, pt) : rect_potential(s, pt);
}

// 静的 : Î = (外側 Gauss)(内側 解析) / (w1 w2)
//   面セルは幅 (= 面積/len) で規格化、細線 (wid = 0) は規格化しない。
//   s2 が細線のときは解析式を持つ側 (多角形) を内側にするため入れ替える
//   (積分は対称なので値は変わらない)。
double poly_static(const seg_t *s1, const seg_t *s2, int nsub)
{
	if ((s2->npv <= 0) && (s2->wid <= 0)) {
		if ((s1->npv <= 0) && (s1->wid <= 0)) return 0;   // 細線どうしは来ない
		return poly_static(s2, s1, nsub);
	}

	double px[POLY_NQMAX * 3], wt[POLY_NQMAX];
	const int np = cell_qpts(s1, nsub, 1, px, wt);

	double sum = 0;
	for (int q = 0; q < np; q++) {
		sum += wt[q] * cell_potential(s2, &px[3 * q]);
	}

	const double w1 = (s1->wid > 0) ? s1->wid : 1;
	const double w2 = (s2->wid > 0) ? s2->wid : 1;

	return sum / (w1 * w2);
}

// 遅延補正 : (1/(w1 w2)) ∬∬ (exp(-jkR) - 1)/R  (両側 3 点三角形則 / 粗い格子)
d_complex_t poly_corr(const seg_t *s1, const seg_t *s2, double kw)
{
	double px1[POLY_NQMAX * 3], wt1[POLY_NQMAX];
	double px2[POLY_NQMAX * 3], wt2[POLY_NQMAX];
	const int n1 = cell_qpts(s1, 1, 0, px1, wt1);
	const int n2 = cell_qpts(s2, 1, 0, px2, wt2);

	double sr = 0, si = 0;
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++) {
			const double dx = px2[3 * j + 0] - px1[3 * i + 0];
			const double dy = px2[3 * j + 1] - px1[3 * i + 1];
			const double dz = px2[3 * j + 2] - px1[3 * i + 2];
			const double r = sqrt((dx * dx) + (dy * dy) + (dz * dz));
			const double w = wt1[i] * wt2[j];
			if (r > 1e-300) {
				sr += w * (cos(kw * r) - 1) / r;
				si -= w * sin(kw * r) / r;
			}
			else {
				si -= w * kw;
			}
		}
	}

	const double w1 = (s1->wid > 0) ? s1->wid : 1;
	const double w2 = (s2->wid > 0) ? s2->wid : 1;

	return d_complex(sr / (w1 * w2), si / (w1 * w2));
}
