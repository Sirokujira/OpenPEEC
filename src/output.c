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

// Z 行列の表 (2 ポート以上のとき)。Z[i][j] = ポート j に 1A、他ポート開放での
// ポート i の電圧。相互インダクタンス M = Im Z21 / omega など、S に戻す手間なく
// 結合量を読めるようにする (1 ポートは入力インピーダンス表と同じなので出さない)。
static void _output_zmat(const peec_t *p, FILE *fp)
{
	const int np = p->nport;

	fprintf(fp, "=== Z matrix === (port j driven by 1 A, other ports open)\n");
	for (int i = 0; i < np; i++) {
	for (int j = 0; j < np; j++) {
		fprintf(fp, "Z%d%d\n", i + 1, j + 1);
		fprintf(fp, "  %s\n", "frequency[Hz]      real[ohm]      imag[ohm]");
		for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
			const d_complex_t z = p->zmat[ZIDX(p, ifreq, i, j)];
			fprintf(fp, "%13.5e%15.6e%15.6e\n", freq_at(p, ifreq), z.r, z.i);
		}
	}
	}

	fflush(fp);
}

void output_spara(const peec_t *p, FILE *fp_log)
{
	if ((p->nport <= 0) || (p->smat == NULL)) return;
	if ((p->nport >= 2) && (p->zmat != NULL)) {
		_output_zmat(p, stdout);
		_output_zmat(p, fp_log);
	}
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
Touchstone 2.0 出力 (peec.ts) — touchstone = 2 のときだけ (peec.sNp と併せて書く)

1.1 との違いは [Reference] でポートごとの基準抵抗を全部書けること。OpenPEEC の
S は電力波の定義でポートごとの z0 を使う (solve.c の z_to_s) ので、z0 が混在する
多ポートは 2.0 でないと正しく記録できない (1.1 は port#1 の値だけ)。
キーワードは [Version] 2.0 / オプション行 / [Number of Ports] /
[Two-Port Data Order] (2 ポートのみ必須、ここでは 12_21 = 行優先) /
[Number of Frequencies] / [Reference] / [Network Data] ... [End]。データは全ポート
数で行優先 S11 S12 ... S1N S21 ... (2 ポートも 1.1 の 21_12 順にしない)。
オプション行の R はポートごとの基準が [Reference] にあるときは無視される仕様だが、
1.1 の読み手のために port#1 の値を書く。
*/
int output_touchstone2(const peec_t *p, const char *fn)
{
	const int np = p->nport;
	if ((p->tsver != 2) || (np <= 0) || (p->smat == NULL)) return 0;

	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fn);
		return 1;
	}

	fprintf(fp, "! %s Ver.%d.%d.%d\n", PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
	fprintf(fp, "! %s\n", p->title);
	fprintf(fp, "[Version] 2.0\n");
	fprintf(fp, "# Hz S RI R %.16g\n", p->port[0].z0);
	fprintf(fp, "[Number of Ports] %d\n", np);
	if (np == 2) {
		fprintf(fp, "[Two-Port Data Order] 12_21\n");
	}
	fprintf(fp, "[Number of Frequencies] %d\n", p->nfreq);
	fprintf(fp, "[Reference]");
	for (int i = 0; i < np; i++) {
		fprintf(fp, " %.16g", p->port[i].z0);
	}
	fprintf(fp, "\n[Network Data]\n");

	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		fprintf(fp, "%.9e", freq_at(p, ifreq));
		for (int i = 0; i < np; i++) {
			if (i > 0) fprintf(fp, "\n%*s", 15, "");
			for (int j = 0; j < np; j++) {
				const d_complex_t sv = p->smat[ZIDX(p, ifreq, i, j)];
				fprintf(fp, " %.9e %.9e", sv.r, sv.i);
			}
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "[End]\n");

	fclose(fp);
	return 0;
}


/*
混合モード S パラメータ (mixed.csv) — mixedmode 指定時のみ

差動対 k = (正側ポート p_k, 負側ポート n_k) について、電力波を
    a_dk = (a_pk - a_nk)/sqrt(2),  a_ck = (a_pk + a_nk)/sqrt(2)
と組み替える (Bockelman & Eisenstadt 1995)。対に入らないポートは単端のまま
通す。この組替え行列 M は直交 (M^-1 = M^T) なので
    S_mm = M S M^T
が混合モード S 行列で、行と列の並びは d1..dP, c1..cP, 単端ポート (番号順)。
例えば 2 対の 4 ポートでは S(d2,d1) が差動伝送 Sdd21、S(c2,d1) が差動 -> 同相
の変換 Scd21 になる。基準抵抗は対の 2 ポートの z0 が等しいことを入力検査で
保証しており、差動モードは 2 z0、同相モードは z0/2 になる (この関係は z0 が
等しいときにだけ成り立つ)。M は直交なので行列全体の |S|^2 の和 (フロベニウス
ノルム) は組替えで変わらない。

書式 : frequency[Hz],row,col,real,imag,|S|[dB] (row/col は d1 / c1 / s3 など。
s の番号は元の単端ポート番号)。peec.log にも S 表と同じ体裁で出す。
*/
int output_mixed(const peec_t *p, const char *fn, FILE *fp_log)
{
	const int np = p->nport;
	if ((p->nmm <= 0) || (np <= 0) || (p->smat == NULL)) return 0;

	const int nmm = p->nmm;
	double *mt = (double *)calloc((size_t)np * np, sizeof(double));
	char (*lab)[16] = malloc((size_t)np * sizeof(*lab));
	double *zref = (double *)malloc((size_t)np * sizeof(double));
	d_complex_t *tmp = (d_complex_t *)malloc((size_t)np * np * sizeof(d_complex_t));
	d_complex_t *smm = (d_complex_t *)malloc((size_t)p->nfreq * np * np * sizeof(d_complex_t));
	int *inpair = (int *)calloc((size_t)np, sizeof(int));
	if ((mt == NULL) || (lab == NULL) || (zref == NULL) || (tmp == NULL) ||
	    (smm == NULL) || (inpair == NULL)) {
		printf("%s\n", "*** memory allocation error (mixed-mode S)");
		free(mt); free(lab); free(zref); free(tmp); free(smm); free(inpair);
		return 1;
	}

	// 組替え行列 M (行 = 混合モードのポート、列 = 単端ポート)
	const double h = 1 / sqrt(2.0);
	for (int k = 0; k < nmm; k++) {
		const int ip = p->mmport[2 * k];
		const int in = p->mmport[(2 * k) + 1];
		inpair[ip] = inpair[in] = 1;
		mt[((size_t)k * np) + ip] = h;
		mt[((size_t)k * np) + in] = -h;
		mt[((size_t)(nmm + k) * np) + ip] = h;
		mt[((size_t)(nmm + k) * np) + in] = h;
		snprintf(lab[k], sizeof(lab[k]), "d%d", k + 1);
		snprintf(lab[nmm + k], sizeof(lab[nmm + k]), "c%d", k + 1);
		zref[k] = 2 * p->port[ip].z0;
		zref[nmm + k] = p->port[ip].z0 / 2;
	}
	int r = 2 * nmm;
	for (int i = 0; i < np; i++) {
		if (inpair[i]) continue;
		mt[((size_t)r * np) + i] = 1;
		snprintf(lab[r], sizeof(lab[r]), "s%d", i + 1);
		zref[r] = p->port[i].z0;
		r++;
	}

	// S_mm = M S M^T (周波数ごと)
	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		// tmp = S M^T
		for (int a = 0; a < np; a++) {
			for (int j = 0; j < np; j++) {
				d_complex_t v = d_complex(0, 0);
				for (int b = 0; b < np; b++) {
					const double m = mt[((size_t)j * np) + b];
					if (m != 0) v = d_add(v, d_rmul(m, p->smat[ZIDX(p, ifreq, a, b)]));
				}
				tmp[((size_t)a * np) + j] = v;
			}
		}
		for (int i = 0; i < np; i++) {
			for (int j = 0; j < np; j++) {
				d_complex_t v = d_complex(0, 0);
				for (int a = 0; a < np; a++) {
					const double m = mt[((size_t)i * np) + a];
					if (m != 0) v = d_add(v, d_rmul(m, tmp[((size_t)a * np) + j]));
				}
				smm[ZIDX(p, ifreq, i, j)] = v;
			}
		}
	}

	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fn);
		free(mt); free(lab); free(zref); free(tmp); free(smm); free(inpair);
		return 1;
	}
	fprintf(fp, "frequency[Hz],row,col,real,imag,|S|[dB]\n");
	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		for (int i = 0; i < np; i++) {
			for (int j = 0; j < np; j++) {
				const d_complex_t s = smm[ZIDX(p, ifreq, i, j)];
				const double mag = d_abs(s);
				const double db = (mag > 1e-30) ? (20 * log10(mag)) : -600;
				fprintf(fp, "%.9e,%s,%s,%.12e,%.12e,%.6f\n",
					freq_at(p, ifreq), lab[i], lab[j], s.r, s.i, db);
			}
		}
	}
	fclose(fp);

	// peec.log : モードの定義と S 表と同じ体裁の表
	fprintf(fp_log, "=== mixed-mode S parameters === (S_mm = M S M^T)\n");
	for (int k = 0; k < nmm; k++) {
		fprintf(fp_log, "d%d / c%d = port#%d - port#%d (Zdiff = %.2f ohm, Zcomm = %.2f ohm)\n",
			k + 1, k + 1, p->mmport[2 * k] + 1, p->mmport[(2 * k) + 1] + 1, zref[k], zref[nmm + k]);
	}
	for (int i = 2 * nmm; i < np; i++) {
		fprintf(fp_log, "%s = single-ended (Z0 = %.2f ohm)\n", lab[i], zref[i]);
	}
	for (int i = 0; i < np; i++) {
		for (int j = 0; j < np; j++) {
			fprintf(fp_log, "S%s%s\n", lab[i], lab[j]);
			fprintf(fp_log, "  %s\n", "frequency[Hz]      real      imag   |S|[dB]  phase[deg]");
			for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
				const d_complex_t s = smm[ZIDX(p, ifreq, i, j)];
				const double mag = d_abs(s);
				const double db = (mag > 1e-30) ? (20 * log10(mag)) : -600;
				fprintf(fp_log, "%13.5e%10.5f%10.5f%10.3f%12.3f\n",
					freq_at(p, ifreq), s.r, s.i, db, atan2(s.i, s.r) * 180 / PI);
			}
		}
	}
	fflush(fp_log);

	free(mt); free(lab); free(zref); free(tmp); free(smm); free(inpair);
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
	// 複数波の同時入射 : Voc は全波の重ね合わせ、leff は 1 波目の振幅で割った値
	for (int iw = 1; iw < p->npw; iw++) {
		fprintf(fp_log, "  + wave #%d : arrival theta = %.2f deg, phi = %.2f deg, %s pol., E0 = %.4g V/m (%.2f deg)\n",
			iw + 1, p->pwi_th[iw], p->pwi_ph[iw], (p->pwi_pol[iw] == 2) ? "phi" : "theta",
			p->pwi_amp[iw], p->pwi_phase[iw]);
	}
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


/*
励振の重み (excitation) による重ね合わせ

ポート j の解は「ポート j に 1 A、他のポートは開放」なので、線形性から
  x = sum_j w_j x_j
が「全ポートに電流 w_j を同時に注入した」解になる (近似なし)。遠方界・
近傍界はこの電流・電荷で評価し、入力電力は
  Pin = (1/2) Re(sum_i conj(w_i) V_i) = (1/2) Re(w^H Z w)
(V = Z w、相互結合 Z_ij を含む)。excitation の行が無ければ port #1 に 1 A
(base をそのまま返す・Pin = Re(Zin_1)/2) で、従来とビット一致する。
*/
const d_complex_t *exc_mix(const peec_t *p, const d_complex_t *base, int n, d_complex_t *buf)
{
	if (p->nexc <= 0) return base;
	for (int m = 0; m < n; m++) {
		buf[m] = d_complex(0, 0);
	}
	for (int k = 0; k < p->nexc; k++) {
		const d_complex_t w = d_polar_deg(p->exc[k].amp, p->exc[k].phase);
		const d_complex_t *xj = &base[(size_t)p->exc[k].port * n];
		for (int m = 0; m < n; m++) {
			buf[m] = d_add(buf[m], d_mul(w, xj[m]));
		}
	}
	return buf;
}

double exc_pin(const peec_t *p, int ifreq)
{
	if (p->nexc <= 0) return 0.5 * p->zin[ifreq].r;
	d_complex_t s = d_complex(0, 0);
	for (int a = 0; a < p->nexc; a++) {
		const d_complex_t wa = d_polar_deg(p->exc[a].amp, p->exc[a].phase);
		for (int b = 0; b < p->nexc; b++) {
			const d_complex_t wb = d_polar_deg(p->exc[b].amp, p->exc[b].phase);
			const d_complex_t z = p->zmat[ZIDX(p, ifreq, p->exc[a].port, p->exc[b].port)];
			s = d_add(s, d_mul(d_conj(wa), d_mul(z, wb)));
		}
	}
	return 0.5 * s.r;
}
