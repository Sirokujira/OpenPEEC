/*
output.c

入力インピーダンス表 (書式は OpenFDTD sol/outputZin.c に準拠、|Z| 列を追加) と
機械読み取り用 CSV (zin.csv)
*/

#include "peec.h"

static void _output_zin(const peec_t *p, FILE *fp)
{
	fprintf(fp, "=== input impedance ===\n");

	for (int iport = 0; iport < p->nport; iport++) {
		fprintf(fp, "port #%d (Z0[ohm] = %.2f)\n", iport + 1, p->port[iport].z0);
		fprintf(fp, "  %s\n", "frequency[Hz] Rin[ohm]   Xin[ohm]    Gin[mS]    Bin[mS]   |Z|[ohm]");
		for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
			const d_complex_t zin = p->zin[(size_t)iport * p->nfreq + ifreq];
			const d_complex_t yin = d_inv(zin);
			fprintf(fp, "%13.5e%11.3f%11.3f%11.3f%11.3f%11.3f\n",
				freq_at(p, ifreq), zin.r, zin.i, yin.r * 1e3, yin.i * 1e3, d_abs(zin));
		}
	}

	fflush(fp);
}

void output_zin(const peec_t *p, FILE *fp_log)
{
	_output_zin(p, stdout);
	_output_zin(p, fp_log);
}

int output_csv(const peec_t *p, const char *fn)
{
	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fn);
		return 1;
	}

	fprintf(fp, "port,frequency[Hz],Rin[ohm],Xin[ohm],Gin[mS],Bin[mS],Zabs[ohm]\n");
	for (int iport = 0; iport < p->nport; iport++) {
		for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
			const d_complex_t zin = p->zin[(size_t)iport * p->nfreq + ifreq];
			const d_complex_t yin = d_inv(zin);
			fprintf(fp, "%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
				iport + 1, freq_at(p, ifreq), zin.r, zin.i, yin.r * 1e3, yin.i * 1e3, d_abs(zin));
		}
	}

	fclose(fp);
	return 0;
}


// S パラメータ表 (振幅 dB と位相 deg)。1 ポートのときは S11 のみ = 反射係数。
static void _output_spara(const peec_t *p, FILE *fp)
{
	const int np = p->nport;

	fprintf(fp, "=== S parameters ===\n");
	fprintf(fp, "reference impedance [ohm] :");
	for (int i = 0; i < np; i++) {
		fprintf(fp, " port#%d=%.2f", i + 1, p->port[i].z0);
	}
	fprintf(fp, "\n");

	for (int i = 0; i < np; i++) {
	for (int j = 0; j < np; j++) {
		fprintf(fp, "S%d%d\n", i + 1, j + 1);
		fprintf(fp, "  %s\n", "frequency[Hz]      real      imag   |S|[dB]  phase[deg]");
		for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
			const d_complex_t s = p->smat[ZIDX(p, ifreq, i, j)];
			const double mag = d_abs(s);
			// |S| = 0 は -inf dB になるので下限で丸める
			const double db = (mag > 1e-30) ? (20 * log10(mag)) : -600;
			fprintf(fp, "%13.5e%10.5f%10.5f%10.3f%12.3f\n",
				freq_at(p, ifreq), s.r, s.i, db, atan2(s.i, s.r) * 180 / PI);
		}
	}
	}

	fflush(fp);
}

void output_spara(const peec_t *p, FILE *fp_log)
{
	if ((p->nport <= 0) || (p->smat == NULL)) return;
	_output_spara(p, stdout);
	_output_spara(p, fp_log);
}


/*
Touchstone 出力 (peec.sNp, N = ポート数)

書式は Touchstone 1.1 の実部/虚部形式:
    # Hz S RI R <z0>
2 ポートのときだけ列順が S11 S21 S12 S22 と転置になるのが仕様なので分岐する
(3 ポート以上は行優先 S11 S12 ... S1N / S21 ...)。

基準抵抗は Touchstone 1.1 では全ポート共通の 1 個しか書けない。ポートごとに
異なる z0 を与えている場合は port#1 の値を書き、その旨を注記する。
*/
int output_touchstone(const peec_t *p)
{
	const int np = p->nport;
	if ((np <= 0) || (p->smat == NULL)) return 0;

	char fn[64];
	snprintf(fn, sizeof(fn), "peec.s%dp", np);

	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fn);
		return 1;
	}

	int mixed = 0;
	for (int i = 1; i < np; i++) {
		if (p->port[i].z0 != p->port[0].z0) mixed = 1;
	}

	fprintf(fp, "! %s Ver.%d.%d.%d\n", PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
	fprintf(fp, "! %s\n", p->title);
	if (mixed) {
		fprintf(fp, "! NOTE: ports have different reference impedances;");
		fprintf(fp, " Touchstone 1.1 can only record one, so port#1 is used.\n");
		for (int i = 0; i < np; i++) {
			fprintf(fp, "!   port#%d z0 = %.6g ohm\n", i + 1, p->port[i].z0);
		}
	}
	fprintf(fp, "# Hz S RI R %.6g\n", p->port[0].z0);

	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		fprintf(fp, "%.9e", freq_at(p, ifreq));
		if (np == 2) {
			// 2 ポートのみ S11 S21 S12 S22 の順 (Touchstone の仕様)
			const int order[4][2] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};
			for (int k = 0; k < 4; k++) {
				const d_complex_t s = p->smat[ZIDX(p, ifreq, order[k][0], order[k][1])];
				fprintf(fp, " %.9e %.9e", s.r, s.i);
			}
			fprintf(fp, "\n");
		}
		else {
			for (int i = 0; i < np; i++) {
				if (i > 0) fprintf(fp, "\n%*s", 15, "");
				for (int j = 0; j < np; j++) {
					const d_complex_t s = p->smat[ZIDX(p, ifreq, i, j)];
					fprintf(fp, " %.9e %.9e", s.r, s.i);
				}
			}
			fprintf(fp, "\n");
		}
	}

	fclose(fp);
	return 0;
}


/*
平面波入射の応答 (pw.csv) — planewave 指定時のみ

ポート自体は MNA に何もスタンプしない (観測点にすぎない) ので、solve() が
捕捉した端子間電圧はそのまま**開放端電圧 Voc** になる。ただしポート間に
素子 (resistor 等) を置いた場合はその負荷が掛かった端子電圧になるので、
Voc として読めるのはポート間に何も繋いでいないときだけ。あわせて

  実効長 |l_eff| = |Voc| / |E0|                       … アンテナの受信有効長
  利用可能電力 Pav = |Voc|^2 / (8 Re Zin)             … 共役整合負荷への電力

を出力する (どちらも開放時の解釈)。Pav はそのポートの Zin を内部インピー
ダンスとする Thevenin 等価なので、多ポートでは他ポートが開放のときの値。
*/
int output_pw(const peec_t *p, const char *fn, FILE *fp_log)
{
	if (!p->pw || (p->voc == NULL)) return 0;

	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fn);
		return 1;
	}

	fprintf(fp, "port,frequency[Hz],Voc_real[V],Voc_imag[V],Voc_abs[V],leff[m],Pav[W]\n");

	fprintf(fp_log, "=== plane wave response === (arrival theta = %.2f deg, phi = %.2f deg, %s pol., E0 = %.4g V/m)\n",
		p->pwth, p->pwph, (p->pwpol == 2) ? "phi" : "theta", p->pwamp);
	for (int iport = 0; iport < p->nport; iport++) {
		fprintf(fp_log, "port #%d\n", iport + 1);
		fprintf(fp_log, "  %s\n", "frequency[Hz]  |Voc|[V]   leff[m]      Pav[W]");
		for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
			const d_complex_t v = p->voc[(size_t)iport * p->nfreq + ifreq];
			const double va = d_abs(v);
			const double rin = p->zin[(size_t)iport * p->nfreq + ifreq].r;
			const double pav = (rin > 0) ? (va * va / (8 * rin)) : 0;
			fprintf(fp_log, "%13.5e%11.4e%11.4e%12.4e\n",
				freq_at(p, ifreq), va, va / p->pwamp, pav);
			fprintf(fp, "%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
				iport + 1, freq_at(p, ifreq), v.r, v.i, va, va / p->pwamp, pav);
		}
	}
	fflush(fp_log);

	fclose(fp);
	return 0;
}


/*
電流・電荷分布 (dist.csv) — distribution = 1 のときのみ

  - I   : ポート励振時の各導体区間の電流 [A] (区間中点の座標つき)
  - Ipw : 平面波入射 (planewave) で誘起された区間電流 [A] — EMC イミュニティ用
  - Q   : ポート励振時の各容量セルの電荷 [C] (capacitance = 1 のときのみ)
を 1 ファイルにまとめる。type 列で行を区別する。

port 列は**その行を励振したポート番号** : I / Q は「ポート j に 1A を注入し
他ポートは開放」の状態なので、多ポートではポートごとに 1 組ずつ出る
(Z 行列の第 j 列に対応するのでクロストークの経路が追える)。Ipw は平面波が
励振源で特定のポートに属さないので 0。

区間電流は MNA 未知数ベクトルの [offS, offS + nseg) がそのまま枝電流なので
追加計算は不要。セル電荷は q = C v (v は各セルのノード電位) で求める。
*/
int output_dist(const peec_t *p, const char *fn)
{
	if (!p->dist) return 0;
	if ((p->segi == NULL) && (p->cellq == NULL) && (p->segipw == NULL)) return 0;

	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fn);
		return 1;
	}

	fprintf(fp, "type,port,index,frequency[Hz],x[m],y[m],z[m],real,imag,abs\n");

	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		const double f = freq_at(p, ifreq);

		// 区間電流 [A] : 座標は区間の中点。励振ポートごとに 1 組
		if (p->segi != NULL) {
			for (int j = 0; j < p->nport; j++) {
				for (int m = 0; m < p->nseg; m++) {
					const seg_t *s = &p->seg[m];
					const d_complex_t v = p->segi[DIDX(p, ifreq, j, m)];
					fprintf(fp, "I,%d,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
						j + 1, m, f,
						0.5 * (s->x1[0] + s->x2[0]),
						0.5 * (s->x1[1] + s->x2[1]),
						0.5 * (s->x1[2] + s->x2[2]),
						v.r, v.i, d_abs(v));
				}
			}
		}

		// 平面波入射による誘起区間電流 [A] (planewave 指定時のみ、port = 0)
		if (p->segipw != NULL) {
			for (int m = 0; m < p->nseg; m++) {
				const seg_t *s = &p->seg[m];
				const d_complex_t v = p->segipw[(size_t)ifreq * p->nseg + m];
				fprintf(fp, "Ipw,0,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
					m, f,
					0.5 * (s->x1[0] + s->x2[0]),
					0.5 * (s->x1[1] + s->x2[1]),
					0.5 * (s->x1[2] + s->x2[2]),
					v.r, v.i, d_abs(v));
			}
		}

		// セル電荷 [C] : 座標はそのセルに属する電荷セルの中点の平均
		if (p->cellq != NULL) {
			for (int j = 0; j < p->nport; j++) {
				for (int m = 0; m < p->ncell; m++) {
					double xc[3] = {0, 0, 0};
					int cnt = 0;
					for (int h = 0; h < p->nchg; h++) {
						if (p->cellof[h] != m) continue;
						for (int d = 0; d < 3; d++) {
							xc[d] += 0.5 * (p->chg[h].x1[d] + p->chg[h].x2[d]);
						}
						cnt++;
					}
					if (cnt > 0) {
						for (int d = 0; d < 3; d++) xc[d] /= cnt;
					}
					const d_complex_t v = p->cellq[QIDX(p, ifreq, j, m)];
					fprintf(fp, "Q,%d,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e\n",
						j + 1, m, f, xc[0], xc[1], xc[2], v.r, v.i, d_abs(v));
				}
			}
		}
	}

	fclose(fp);
	return 0;
}
