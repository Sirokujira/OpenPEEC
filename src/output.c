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
