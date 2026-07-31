/*
skin.c

円形断面導体の内部インピーダンス (表皮効果 + 内部インダクタンス)

  Z_int(単位長) = k / (2 pi a sigma) * I0(ka) / I1(ka),   k = (1 + j) / delta
  delta = sqrt(2 / (omega mu0 sigma))                      … 表皮深さ

極限値 (検証に使用) :
  低周波 : Z -> R_dc (1 + q^4/48) + j omega mu0/(8 pi)      (q = a/delta << 1)
  高周波 : R = X -> 1 / (2 pi a delta sigma)                (q >> 1)

I0/I1 は |ka| < 15 で級数、それ以上で漸近展開
  I0(z)/I1(z) = 1 + 1/(2z) + 3/(8z^2) + 3/(8z^3) + O(z^-4)
を用いる (切り替え点での相対誤差は 1e-5 程度)。
*/

#include "peec.h"

#define ZSERIES_MAX 15.0
#define NTERM 60

static d_complex_t i0_over_i1(d_complex_t z)
{
	if (d_abs(z) < ZSERIES_MAX) {
		// 級数 : I0 = sum (z^2/4)^k / (k!)^2, I1 = (z/2) sum (z^2/4)^k / (k! (k+1)!)
		const d_complex_t z2 = d_mul(z, z);

		d_complex_t i0 = d_complex(1, 0);
		d_complex_t t = d_complex(1, 0);
		for (int k = 1; k <= NTERM; k++) {
			t = d_rmul(1.0 / (4.0 * k * k), d_mul(t, z2));
			i0 = d_add(i0, t);
			if (d_abs(t) < 1e-17 * d_abs(i0)) break;
		}

		d_complex_t s = d_complex(1, 0);
		t = d_complex(1, 0);
		for (int k = 1; k <= NTERM; k++) {
			t = d_rmul(1.0 / (4.0 * k * (k + 1)), d_mul(t, z2));
			s = d_add(s, t);
			if (d_abs(t) < 1e-17 * d_abs(s)) break;
		}
		const d_complex_t i1 = d_mul(d_rmul(0.5, z), s);

		return d_div(i0, i1);
	}

	// 漸近展開
	const d_complex_t u = d_inv(z);
	const d_complex_t u2 = d_mul(u, u);
	const d_complex_t u3 = d_mul(u2, u);

	return d_add(d_add(d_complex(1, 0), d_rmul(0.5, u)),
	             d_add(d_rmul(0.375, u2), d_rmul(0.375, u3)));
}

// 長さ len、半径 a、導電率 sigma の丸線 1 本分の内部インピーダンス
d_complex_t zint_round(double len, double a, double sigma, double freq)
{
	if ((sigma <= 0) || (len <= 0) || (a <= 0)) return d_complex(0, 0);

	const double rdc = len / (sigma * PI * a * a);
	const double omega = 2 * PI * freq;
	if (omega <= 0) return d_complex(rdc, 0);

	const double delta = sqrt(2 / (omega * MU0 * sigma));
	const double q = a / delta;

	// ka = (1 + j) a / delta
	const d_complex_t ratio = i0_over_i1(d_complex(q, q));

	// k / (2 pi a sigma) = (1 + j) / (2 pi a delta sigma)
	const double kf = 1 / (2 * PI * a * delta * sigma);

	return d_rmul(len, d_mul(d_complex(kf, kf), ratio));
}

/*
角線 (矩形断面) の内部インピーダンス

丸線のような閉形式が無いため、両極限が厳密になる工学的な合成式を使う :
  低周波 : R -> R_dc = len/(sigma A)、L_int -> mu0 len/(8 pi)
  高周波 : 電流は周長 P の表皮層 (深さ delta) を流れ
           R = X/omega*... -> len/(sigma delta P)  (R = X)
合成 : R = sqrt(R_dc^2 + R_hf^2)、L_int = L_dc L_hf / sqrt(L_dc^2 + L_hf^2)
       (どちらも滑らかで単調、両極限で厳密)
遷移域 (delta ~ 断面寸法) は近似なので、その領域の精度が要る用途では
丸線モデル (SHAPE_ROUND) を使うこと。
*/
d_complex_t zint_bar(double len, double area, double perim, double sigma, double freq)
{
	if ((sigma <= 0) || (len <= 0) || (area <= 0) || (perim <= 0)) return d_complex(0, 0);

	const double rdc = len / (sigma * area);
	const double omega = 2 * PI * freq;
	if (omega <= 0) return d_complex(rdc, 0);

	const double delta = sqrt(2 / (omega * MU0 * sigma));
	const double rhf = len / (sigma * delta * perim);

	const double ldc = MU0 * len / (8 * PI);
	const double lhf = rhf / omega;

	const double r = sqrt((rdc * rdc) + (rhf * rhf));
	const double l = ldc * lhf / sqrt((ldc * ldc) + (lhf * lhf));

	return d_complex(r, omega * l);
}

// 区間の内部インピーダンス (断面形状で振り分ける)
d_complex_t zint_seg(const seg_t *s, double freq)
{
	// 体積セル (plate の厚み分割) : 厚み方向の電流再配分は Lp 行列が
	// 陽に解くので、合成式を使うと二重計上になる。DC 抵抗のみ返す
	// (内部インダクタンスも体積自己項に含まれている)。
	if (s->vol) return d_complex(s->res, 0);

	// 多角形セル (三角形メッシュ) : DC 抵抗は cotangent 重み (FEM) で決まる。
	// 合成式は幾何断面を仮定していて DC 極限が cotangent 抵抗とずれるので
	// 適用しない (面内の電流再配分はメッシュが解く。厚み方向は未分解)。
	if (s->shape == SHAPE_POLY) return d_complex(s->res, 0);

	// 丸線のみ Bessel 閉形式。角線・面セル (plate) は断面 (幅 x 厚さ) の
	// 合成式を使う。以前は SHAPE_BAR 以外がすべて丸線式に落ち、面セルは
	// radius = 0 のため内部インピーダンスが 0 になっていた。
	return (s->shape == SHAPE_ROUND)
		? zint_round(s->len, s->radius, s->sigma, freq)
		: zint_bar(s->len, s->area, s->perim, s->sigma, freq);
}
