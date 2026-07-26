/*
solve.c

周波数掃引 : 各周波数で MNA 行列を組み立て LU 分解し、
ポートごとに 1A 電流源励振で Zin を求める。
独立電源 (vsource/isource) がある場合は追加で 1 回解いて節点電圧をログ出力する。
*/

#include "peec.h"

static void print_vnode(const peec_t *p, const d_complex_t *x, double f, FILE *fp_log)
{
	fprintf(fp_log, "=== node voltage === (f = %.5e Hz)\n", f);
	fprintf(fp_log, "  %s\n", "node      real[V]      imag[V]       |V|[V]");
	for (int id = 0; id <= p->maxid; id++) {
		const int m = p->nodemap[id];
		if (m == -2) continue;
		const d_complex_t v = (m < 0) ? d_complex(0, 0) : x[m];
		fprintf(fp_log, "%6d %12.4e %12.4e %12.4e\n", id, v.r, v.i, d_abs(v));
	}
}

int solve(peec_t *p, FILE *fp_log)
{
	const int n = p->nunknown;

	d_complex_t *a = (d_complex_t *)malloc((size_t)n * n * sizeof(d_complex_t));
	d_complex_t *b = (d_complex_t *)malloc((size_t)n * sizeof(d_complex_t));
	int *piv = (int *)malloc((size_t)n * sizeof(int));
	p->zin = (d_complex_t *)malloc((size_t)p->nport * p->nfreq * sizeof(d_complex_t));
	if ((a == NULL) || (b == NULL) || (piv == NULL) || (p->zin == NULL)) {
		printf("%s\n", "*** memory allocation error (matrix)");
		return 1;
	}

	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		const double f = freq_at(p, ifreq);

		// 遅延ありのときは部分要素が周波数依存になるので毎回作り直す
		if (p->retardation) {
			lp_fill(p, f, fp_log);
			if (pot_fill(p, f, fp_log)) return 1;
		}

		mna_assemble(p, f, a);
		const int ising = lu_decomp(n, a, piv);
		if (ising >= 0) {
			printf("*** singular matrix at f = %.5e Hz (row %d : floating node or ideal voltage source loop?)\n", f, ising);
			fprintf(fp_log, "*** singular matrix at f = %.5e Hz (row %d)\n", f, ising);
			return 1;
		}

		for (int iport = 0; iport < p->nport; iport++) {
			mna_rhs_port(p, iport, b);
			lu_solve(n, a, piv, b);
			const int i1 = p->nodemap[p->port[iport].n1];
			const int i2 = p->nodemap[p->port[iport].n2];
			const d_complex_t v1 = (i1 < 0) ? d_complex(0, 0) : b[i1];
			const d_complex_t v2 = (i2 < 0) ? d_complex(0, 0) : b[i2];
			p->zin[(size_t)iport * p->nfreq + ifreq] = d_sub(v1, v2);
		}

		if (p->nsrc > 0) {
			mna_rhs_sources(p, b);
			lu_solve(n, a, piv, b);
			print_vnode(p, b, f, fp_log);
		}
	}

	free(a);
	free(b);
	free(piv);

	return 0;
}
