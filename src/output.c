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
