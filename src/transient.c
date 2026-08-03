/*
transient.c

過渡応答 (transient) — 周波数掃引の逆フーリエ変換による時間波形

周波数領域ソルバーの結果 H(f) に励振パルスのスペクトル X(f) を掛けて
時間波形へ戻す :

  y(t) = ∫ H(f) X(f) e^{j 2 pi f t} df
       ≈ df [ H(0)X(0) + 2 Re Σ_{k=1..N} H(f_k) X(f_k) e^{j 2 pi f_k t} ]

掃引点は f_k = k df (k = 1..N) でなければならない (DC の整数倍で等間隔)。
`frequency = df  N*df  N-1` と書けばこの条件になる。対数掃引は不可。

励振はガウスパルス
  x(t) = exp(-(t - t0)^2 / (2 sigma^2))        (ピーク振幅 1)
  X(f) = sigma sqrt(2 pi) exp(-2 pi^2 sigma^2 f^2) exp(-j 2 pi f t0)
帯域端 f_max = f_N での減衰が att [dB] になるよう sigma を決める :
  2 pi^2 sigma^2 f_max^2 = att ln(10)/20
これで帯域打ち切りによるリンギング (Gibbs) を抑える。t0 = 4 sigma とすると
t < 0 側の裾は ~1e-7 で、因果性の破れは表示上問題にならない。

時間軸は df と f_max だけで決まる :
  サンプル数 M = 2N、 dt = 1/(2 f_max)、 全長 T = 1/df
応答が T 以内に減衰しない (Q の高い共振) と巻き込み (時間領域エイリアス)
が起きる。その場合は分割数を増やして df を小さくする。

DC 項 H(0) は掃引に含まれない (f > 0 が必須) ので、実数条件 Im H(0) = 0 の
もとで下 2 点から線形補外する (H(0) = 2 H(f_1) - H(f_2))。実数系の伝達
関数は DC で実数になるので、これは物理的にも正しい形。

出力 (tran.csv) は long 形式 :
  X : 励振パルスそのもの (表示の参照用。H = 1 として同じ合成式で作る)
  S : S 行列各成分の時間応答 (i = j なら TDR、i != j なら伝送)
  V : 平面波入射による誘起端子電圧の波形 [V] (planewave 指定時)
X と S を同じ合成式で作るので、周波数依存の無い系 (純抵抗) では
y(t) = S x(t) が打ち切り誤差込みで厳密に成り立つ。
*/

#include "peec.h"

// 掃引が f_k = k df (等間隔かつ DC の整数倍) になっているか
static int tran_check(const peec_t *p, FILE *fp_log)
{
	if (p->flog) {
		printf("%s\n", "*** transient : logarithmic sweep is not supported (use a linear sweep)");
		fprintf(fp_log, "%s\n", "*** transient : logarithmic sweep is not supported");
		return 1;
	}
	if (p->nfreq < 3) {
		printf("%s\n", "*** transient : at least 3 frequency points are required");
		fprintf(fp_log, "%s\n", "*** transient : at least 3 frequency points are required");
		return 1;
	}
	const double df = (p->f1 - p->f0) / (p->nfreq - 1);
	if (fabs(p->f0 - df) > 1e-6 * df) {
		printf("*** transient : sweep must be f_k = k df (start = spacing); f0 = %.6e, df = %.6e\n",
			p->f0, df);
		fprintf(fp_log, "*** transient : sweep must be f_k = k df (start = spacing); f0 = %.6e, df = %.6e\n",
			p->f0, df);
		return 1;
	}

	return 0;
}

/*
時間波形の合成。
  h    : [nfreq] 周波数応答 (NULL なら H = 1 : 励振パルスそのもの)
  xs   : [nfreq] 励振スペクトル、x0 = その DC 値 (実数)
  y    : [m] 出力波形
*/
static void tran_synth(const peec_t *p, const d_complex_t *h, const d_complex_t *xs,
	double x0, double df, double dt, int m, double *y)
{
	// DC : H(0) は実数条件のもとで下 2 点から線形補外
	const double h0 = (h == NULL) ? 1.0 : ((2 * h[0].r) - h[1].r);
	const double g0 = h0 * x0;

	for (int n = 0; n < m; n++) {
		const double t = n * dt;
		double sum = g0;
		for (int k = 0; k < p->nfreq; k++) {
			const d_complex_t g = (h == NULL) ? xs[k] : d_mul(h[k], xs[k]);
			const double ang = 2 * PI * freq_at(p, k) * t;
			sum += 2 * ((g.r * cos(ang)) - (g.i * sin(ang)));
		}
		y[n] = df * sum;
	}
}

int output_tran(const peec_t *p, const char *fn, FILE *fp_log)
{
	if (!p->tran) return 0;
	if (tran_check(p, fp_log)) return 1;

	const int nf = p->nfreq;
	const double df = p->f0;                  // = 掃引間隔 (tran_check で確認済み)
	const double fmax = p->f1;
	const int m = 2 * nf;                     // 時間サンプル数
	const double dt = 1 / (2 * fmax);

	// ガウスパルス : 帯域端で att [dB] 減衰
	const double sigma = sqrt(p->tranatt * log(10.0) / (20 * 2 * PI * PI)) / fmax;
	const double t0 = 4 * sigma;

	d_complex_t *xs = (d_complex_t *)malloc((size_t)nf * sizeof(d_complex_t));
	d_complex_t *hh = (d_complex_t *)malloc((size_t)nf * sizeof(d_complex_t));
	double *y = (double *)malloc((size_t)m * sizeof(double));
	if ((xs == NULL) || (hh == NULL) || (y == NULL)) {
		printf("%s\n", "*** memory allocation error (transient)");
		free(xs); free(hh); free(y);
		return 1;
	}

	// 励振スペクトル X(f) = sigma sqrt(2pi) exp(-2 pi^2 sigma^2 f^2) e^{-j 2 pi f t0}
	const double x0 = sigma * sqrt(2 * PI);
	for (int k = 0; k < nf; k++) {
		const double f = freq_at(p, k);
		const double a = x0 * exp(-2 * PI * PI * sigma * sigma * f * f);
		const double ph = -2 * PI * f * t0;
		xs[k] = d_complex(a * cos(ph), a * sin(ph));
	}

	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fn);
		free(xs); free(hh); free(y);
		return 1;
	}
	fprintf(fp, "type,i,j,time[s],value\n");

	// 励振パルス (参照用。S と同じ合成式なので打ち切り誤差も共通)
	tran_synth(p, NULL, xs, x0, df, dt, m, y);
	for (int n = 0; n < m; n++) {
		fprintf(fp, "X,0,0,%.9e,%.9e\n", n * dt, y[n]);
	}

	// S 行列各成分の時間応答
	if (p->smat != NULL) {
		for (int i = 0; i < p->nport; i++) {
		for (int j = 0; j < p->nport; j++) {
			for (int k = 0; k < nf; k++) {
				hh[k] = p->smat[ZIDX(p, k, i, j)];
			}
			tran_synth(p, hh, xs, x0, df, dt, m, y);
			for (int n = 0; n < m; n++) {
				fprintf(fp, "S,%d,%d,%.9e,%.9e\n", i + 1, j + 1, n * dt, y[n]);
			}
		}
		}
	}

	// 平面波入射による誘起端子電圧の波形 [V]
	if (p->voc != NULL) {
		for (int i = 0; i < p->nport; i++) {
			for (int k = 0; k < nf; k++) {
				hh[k] = p->voc[(size_t)i * nf + k];
			}
			tran_synth(p, hh, xs, x0, df, dt, m, y);
			for (int n = 0; n < m; n++) {
				fprintf(fp, "V,%d,0,%.9e,%.9e\n", i + 1, n * dt, y[n]);
			}
		}
	}

	fprintf(fp_log, "transient : %d samples, dt = %.5e s, span = %.5e s, gaussian sigma = %.5e s (t0 = %.5e s, -%.1f dB at %.5e Hz)\n",
		m, dt, m * dt, sigma, t0, p->tranatt, fmax);
	fflush(fp_log);

	free(xs);
	free(hh);
	free(y);
	fclose(fp);

	return 0;
}
