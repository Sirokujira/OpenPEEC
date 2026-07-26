/*
surface.c

面導体セル (リボン = 幅 wid の帯) どうしの幾何二重積分

細線と同じ形の量になるよう、幅で規格化した積分を返す :

  Ihat(S1,S2) = 1/(w1 w2) ∬_{S1} ∬_{S2} exp(-j k R) / R dS dS'

こうすると w -> 0 で細線の ∬ dl dl'/R に一致し、
  Lp = (mu0/4pi)(t1・t2) Ihat,  P = Ihat/(4 pi eps0 A1 A2 /(w1 w2)) = Ihat/(4 pi eps0 L1 L2)
と、細線とまったく同じ式が使える。

静的部は内側の面積分を解析的に評価する :
  ∬_rect dS' / |p - r'| は矩形の閉形式 (Hess & Smith 型) で厳密に書ける。
  この値は矩形上でも有界で連続なので、外側だけを Gauss 求積すれば
  自己項・隣接項も安定に計算できる。

遅延補正 (exp(-jkR)-1)/R は R -> 0 で -jk に収束する正則関数なので、
4 次元 Gauss 求積をそのまま適用する。
*/

#include "peec.h"

// 8 点 / 4 点 Gauss-Legendre ([-1,1])
static const double xg8[8] = {
	-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
	 0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
static const double wg8[8] = {
	 0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620,
	 0.3626837833783620,  0.3137066458778873,  0.2223810344533745,  0.1012285362903763};
static const double xg4[4] = {
	-0.8611363115940526, -0.3399810435848563,  0.3399810435848563,  0.8611363115940526};
static const double wg4[4] = {
	 0.3478548451374538,  0.6521451548625461,  0.6521451548625461,  0.3478548451374538};

static double dot3(const double *a, const double *b)
{
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

// 矩形上の一様分布に対するポテンシャル核の原始関数
//   Phi(u,v,z) = u asinh(v/sqrt(u^2+z^2)) + v asinh(u/sqrt(v^2+z^2))
//                - |z| atan2(u v, |z| sqrt(u^2+v^2+z^2))
// asinh 形にすることで z = 0 の面上でも安全に評価できる。
static double phi_rect(double u, double v, double z)
{
	const double az = fabs(z);
	const double uz = (u * u) + (z * z);
	const double vz = (v * v) + (z * z);
	double t = 0;

	if (uz > 0) t += u * asinh(v / sqrt(uz));
	if (vz > 0) t += v * asinh(u / sqrt(vz));
	if (az > 0) t -= az * atan2(u * v, az * sqrt(uz + (v * v)));

	return t;
}

// セル s の局所座標系 (軸 t, 横 wv, 法線 n) と中心
static void local_frame(const seg_t *s, double t[3], double w[3], double nv[3], double c[3])
{
	for (int i = 0; i < 3; i++) {
		t[i] = (s->x2[i] - s->x1[i]) / s->len;
		w[i] = s->wv[i];
		c[i] = 0.5 * (s->x1[i] + s->x2[i]);
	}
	nv[0] = (t[1] * w[2]) - (t[2] * w[1]);
	nv[1] = (t[2] * w[0]) - (t[0] * w[2]);
	nv[2] = (t[0] * w[1]) - (t[1] * w[0]);
}

// 点 pt から見た矩形セル s 上の ∬ dS' / R (解析式)
double rect_potential(const seg_t *s, const double *pt)
{
	double t[3], w[3], nv[3], c[3], d[3];
	local_frame(s, t, w, nv, c);
	for (int i = 0; i < 3; i++) {
		d[i] = pt[i] - c[i];
	}
	const double u = dot3(d, t);
	const double v = dot3(d, w);
	const double z = dot3(d, nv);
	const double hl = 0.5 * s->len;
	const double hw = 0.5 * s->wid;

	return phi_rect(u + hl, v + hw, z) - phi_rect(u - hl, v + hw, z)
	     - phi_rect(u + hl, v - hw, z) + phi_rect(u - hl, v - hw, z);
}

// 静的 : Ihat = (1/(w1 w2)) ∬_{S1} [解析的な内側積分] dS
//   外側は S1 を nsub x nsub に分割した複合 8x8 点 Gauss 求積
double ribbon_static(const seg_t *s1, const seg_t *s2, int nsub)
{
	double t1[3], w1[3], n1[3], c1[3];
	local_frame(s1, t1, w1, n1, c1);

	double sum = 0;
	for (int ia = 0; ia < nsub; ia++) {
	for (int ib = 0; ib < nsub; ib++) {
		for (int i = 0; i < 8; i++) {
			// 軸方向 : 中心から測った [-len/2, len/2]
			const double a = s1->len * (((ia + (0.5 * (1 + xg8[i]))) / nsub) - 0.5);
			for (int j = 0; j < 8; j++) {
				const double b = s1->wid * (((ib + (0.5 * (1 + xg8[j]))) / nsub) - 0.5);
				double pt[3];
				for (int k = 0; k < 3; k++) {
					pt[k] = c1[k] + (a * t1[k]) + (b * w1[k]);
				}
				sum += wg8[i] * wg8[j] * rect_potential(s2, pt);
			}
		}
	}
	}

	// ∬_{S1} dS = (len1 wid1 / 4) sum(w_i w_j ...) / nsub^2、これを w1 w2 で割る
	return s1->len * sum / (4 * s2->wid * (double)nsub * nsub);
}

// 遅延補正 : (1/(w1 w2)) ∬∬ (exp(-jkR) - 1)/R dS dS'  (4 次元 Gauss、4 点則)
d_complex_t ribbon_corr(const seg_t *s1, const seg_t *s2, double kw, int nsub)
{
	double t1[3], w1[3], n1[3], c1[3], t2[3], w2[3], n2[3], c2[3];
	local_frame(s1, t1, w1, n1, c1);
	local_frame(s2, t2, w2, n2, c2);

	double sr = 0, si = 0;
	for (int ia = 0; ia < nsub; ia++) {
	for (int ja = 0; ja < nsub; ja++) {
		for (int i = 0; i < 4; i++) {
			const double a1 = s1->len * (((ia + (0.5 * (1 + xg4[i]))) / nsub) - 0.5);
			for (int j = 0; j < 4; j++) {
				const double b1 = s1->wid * ((0.5 * (1 + xg4[j])) - 0.5);
				double pt[3];
				for (int k = 0; k < 3; k++) {
					pt[k] = c1[k] + (a1 * t1[k]) + (b1 * w1[k]);
				}
				for (int m = 0; m < 4; m++) {
					const double a2 = s2->len * (((ja + (0.5 * (1 + xg4[m]))) / nsub) - 0.5);
					for (int q = 0; q < 4; q++) {
						const double b2 = s2->wid * ((0.5 * (1 + xg4[q])) - 0.5);
						const double dx = c2[0] + (a2 * t2[0]) + (b2 * w2[0]) - pt[0];
						const double dy = c2[1] + (a2 * t2[1]) + (b2 * w2[1]) - pt[1];
						const double dz = c2[2] + (a2 * t2[2]) + (b2 * w2[2]) - pt[2];
						const double r = sqrt((dx * dx) + (dy * dy) + (dz * dz));
						const double w = wg4[i] * wg4[j] * wg4[m] * wg4[q];
						if (r > 1e-300) {
							sr += w * (cos(kw * r) - 1) / r;
							si -= w * sin(kw * r) / r;
						}
						else {
							si -= w * kw;
						}
					}
				}
			}
		}
	}
	}

	// ∬∬ dS dS' = (len1 wid1/4)(len2 wid2/4) sum(...) / nsub^2、これを w1 w2 で割る
	const double fac = s1->len * s2->len / (16 * (double)nsub * nsub);

	return d_complex(fac * sr, fac * si);
}
