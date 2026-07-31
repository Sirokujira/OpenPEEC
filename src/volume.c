/*
volume.c

体積セル (矩形バー) どうしの幾何二重積分 — Hoer-Love の閉形式

  C. Hoer and C. Love, "Exact inductance equations for rectangular
  conductors with applications to more complicated geometries",
  J. Res. NBS 69C(2), pp.127-137, 1965.

平行で断面軸もそろった 2 本の矩形バー
  bar1 : [0,a]x[0,b]x[0,c]、bar2 : [E,E+d]x[P,P+e]x[l3,l3+f]
に対して (z が軸方向、断面は a x b / d x e)

  U = ∬∬ dV dV' / R  =  Σ_{i,j,k=1..4} (-1)^(i+j+k+1) F(q_i, r_j, s_k)
  q = {E-a, E+d-a, E+d, E}
  r = {P-b, P+e-b, P+e, P}
  s = {l3-c, l3+f-c, l3+f, l3}

面セル (surface.c) と同様、断面積で規格化した Î = U/(A1 A2) を返す。
これにより thick -> 0 でリボンの Î に一致し、
  Lp = (mu0/4pi)(t1・t2) Î
が細線・リボンとまったく同じ式で使える。

自己項は E = P = l3 = 0 (2 本が完全に重なる) とした同じ式で厳密に求まる。
重なり・面共有・分離のどの配置でも有効で、体積を分割した部分和が全体に
一致する加法性が機械精度で成り立つ (検証済み)。

数値上の注意 :
- 対数項は log((x+R)/sqrt(y^2+z^2)) ではなく asinh(x/sqrt(y^2+z^2)) と書く
  (等価だが x < 0 で桁落ちしない)。
- ゼロ引数では「係数が消える項だけ」を 0 にする。たとえば F(0,y,z) は
  x を含む 4 項が消えるだけで残りは有限に残る。引数が 0 だからといって
  F 全体を 0 にすると、重なり領域・自己項が不連続に狂う (実際に踏んだ)。
- 遠く離れた対では 64 項の交代和が桁落ちする (F ~ r^5 に対し U ~ V1 V2/r)。
  中心間距離が寸法の HL_RFAR 倍を超えたら呼び出し側でリボン積分に
  切り替える (bar_use_hl)。その距離では厚みの寄与は O((t/r)^2) で無視できる。
*/

#include "peec.h"

#define HL_RFAR 8.0

// Hoer-Love の原始関数 F(x,y,z)
static double hl_F(double x, double y, double z)
{
	const double x2 = x * x;
	const double y2 = y * y;
	const double z2 = z * z;
	const double R = sqrt(x2 + y2 + z2);
	const double yz = y2 + z2;
	const double xz = x2 + z2;
	const double xy = x2 + y2;
	double t = 0;

	if ((x != 0) && (yz > 0)) {
		t += ((y2 * z2 / 4) - (y2 * y2 / 24) - (z2 * z2 / 24)) * x * asinh(x / sqrt(yz));
	}
	if ((y != 0) && (xz > 0)) {
		t += ((x2 * z2 / 4) - (x2 * x2 / 24) - (z2 * z2 / 24)) * y * asinh(y / sqrt(xz));
	}
	if ((z != 0) && (xy > 0)) {
		t += ((x2 * y2 / 4) - (x2 * x2 / 24) - (y2 * y2 / 24)) * z * asinh(z / sqrt(xy));
	}
	t += ((x2 * x2) + (y2 * y2) + (z2 * z2)
	    - 3 * ((x2 * y2) + (y2 * z2) + (z2 * x2))) * R / 60;
	if (z != 0) t -= (x * y * z * z2 / 6) * atan((x * y) / (z * R));
	if (y != 0) t -= (x * y * y2 * z / 6) * atan((x * z) / (y * R));
	if (x != 0) t -= (x * x2 * y * z / 6) * atan((y * z) / (x * R));

	return t;
}

static double dot3(const double *a, const double *b)
{
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

// セルの局所枠 : 軸 t、幅 wv、厚み nv = t x wv、中心 c
static void vol_frame(const seg_t *s, double t[3], double w[3], double nv[3], double c[3])
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

// 閉形式が使える対か : 軸・幅方向がともに (±) 平行で、遠すぎない
int bar_use_hl(const seg_t *s1, const seg_t *s2)
{
	double t1[3], w1[3], n1[3], c1[3], t2[3], w2[3], n2[3], c2[3];
	vol_frame(s1, t1, w1, n1, c1);
	vol_frame(s2, t2, w2, n2, c2);

	if (fabs(dot3(t1, t2)) < 1 - 1e-9) return 0;
	if (fabs(dot3(w1, w2)) < 1 - 1e-9) return 0;

	const double dx = c2[0] - c1[0];
	const double dy = c2[1] - c1[1];
	const double dz = c2[2] - c1[2];
	const double rc = sqrt((dx * dx) + (dy * dy) + (dz * dz));
	double lmax = (s1->len > s2->len) ? s1->len : s2->len;
	if (s1->wid > lmax) lmax = s1->wid;
	if (s2->wid > lmax) lmax = s2->wid;

	return (rc <= HL_RFAR * lmax);
}

// 規格化した幾何二重積分 Î = (1/(A1 A2)) ∬∬ dV dV'/R (静的、常に正)
// s1 == s2 (自己項) もそのまま扱える。
double bar_pair(const seg_t *s1, const seg_t *s2)
{
	double t1[3], w1[3], n1[3], c1[3], t2[3], w2[3], n2[3], c2[3];
	vol_frame(s1, t1, w1, n1, c1);
	vol_frame(s2, t2, w2, n2, c2);

	// bar1 の枠 (x = 厚み n1、y = 幅 w1、z = 軸 t1) で表す。
	// 軸がそろっているので bar2 の寸法はそのまま各軸に沿う
	// (逆向きでも矩形は対称なので寸法・中心だけで決まる)。
	const double a = s1->thick;
	const double b = s1->wid;
	const double c = s1->len;
	const double d = s2->thick;
	const double e = s2->wid;
	const double f = s2->len;

	double rc[3];
	for (int i = 0; i < 3; i++) {
		rc[i] = c2[i] - c1[i];
	}
	// 角どうしのオフセット (Hoer-Love の E, P, l3)
	const double E  = dot3(rc, n1) + ((a - d) / 2);
	const double P  = dot3(rc, w1) + ((b - e) / 2);
	const double l3 = dot3(rc, t1) + ((c - f) / 2);

	const double q[4] = {E - a, E + d - a, E + d, E};
	const double r[4] = {P - b, P + e - b, P + e, P};
	const double s[4] = {l3 - c, l3 + f - c, l3 + f, l3};

	double sum = 0;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			for (int k = 0; k < 4; k++) {
				// 1 起点の (-1)^(i+j+k+1) は 0 起点では (i+j+k) 偶数で +1
				const double sgn = (((i + j + k) % 2) == 0) ? 1.0 : -1.0;
				sum += sgn * hl_F(q[i], r[j], s[k]);
			}
		}
	}

	return sum / (a * b * d * e);
}
