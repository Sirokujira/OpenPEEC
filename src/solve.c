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

/*
Z 行列 -> S 行列 (電力波の定義、Kurokawa)

    S = F (Z - Z0) (Z + Z0)^-1 F^-1,   F = diag(1 / (2 sqrt(z0_i)))

基準抵抗 z0 は実数 (port キーの第 3 引数) なので、成分で書くと

    S_ij = sqrt(z0_j / z0_i) * [(Z - Z0)(Z + Z0)^-1]_ij

全ポートの z0 が等しいときは教科書どおりの (Z - Z0)(Z + Z0)^-1 に一致する。
(Z + Z0) の逆行列は既存の複素 LU で単位ベクトルを解いて作る (ポート数は小さい)。
戻り値 : 0 = 正常、非 0 = (Z + Z0) が特異
*/
static int z_to_s(peec_t *p, int ifreq, d_complex_t *m, d_complex_t *minv,
	d_complex_t *col, int *piv)
{
	const int np = p->nport;

	// M = Z + Z0
	for (int i = 0; i < np; i++) {
		for (int j = 0; j < np; j++) {
			d_complex_t v = p->zmat[ZIDX(p, ifreq, i, j)];
			if (i == j) v = d_add(v, d_complex(p->port[i].z0, 0));
			m[(size_t)i * np + j] = v;
		}
	}

	if (lu_decomp(np, m, piv) >= 0) return 1;

	// M^-1 を列ごとに求める
	for (int k = 0; k < np; k++) {
		memset(col, 0, (size_t)np * sizeof(d_complex_t));
		col[k] = d_complex(1, 0);
		lu_solve(np, m, piv, col);
		for (int i = 0; i < np; i++) {
			minv[(size_t)i * np + k] = col[i];
		}
	}

	// S = sqrt(z0_j / z0_i) * (Z - Z0) M^-1
	for (int i = 0; i < np; i++) {
		for (int j = 0; j < np; j++) {
			d_complex_t s = d_complex(0, 0);
			for (int k = 0; k < np; k++) {
				d_complex_t nk = p->zmat[ZIDX(p, ifreq, i, k)];
				if (i == k) nk = d_sub(nk, d_complex(p->port[i].z0, 0));
				s = d_add(s, d_mul(nk, minv[(size_t)k * np + j]));
			}
			const double f = sqrt(p->port[j].z0 / p->port[i].z0);
			p->smat[ZIDX(p, ifreq, i, j)] = d_rmul(f, s);
		}
	}

	return 0;
}


int solve(peec_t *p, FILE *fp_log)
{
	const int n = p->nunknown;

	d_complex_t *a = (d_complex_t *)malloc((size_t)n * n * sizeof(d_complex_t));
	d_complex_t *b = (d_complex_t *)malloc((size_t)n * sizeof(d_complex_t));
	int *piv = (int *)malloc((size_t)n * sizeof(int));
	p->zin = (d_complex_t *)malloc((size_t)p->nport * p->nfreq * sizeof(d_complex_t));

	// 多ポートの Z / S 行列と、Z -> S 変換の作業領域
	const size_t nmat = (size_t)p->nfreq * p->nport * p->nport;
	const int np = p->nport;
	p->zmat = (d_complex_t *)malloc(nmat * sizeof(d_complex_t));
	p->smat = (d_complex_t *)malloc(nmat * sizeof(d_complex_t));
	d_complex_t *sm   = (d_complex_t *)malloc((size_t)np * np * sizeof(d_complex_t));
	d_complex_t *sinv = (d_complex_t *)malloc((size_t)np * np * sizeof(d_complex_t));
	d_complex_t *scol = (d_complex_t *)malloc((size_t)np * sizeof(d_complex_t));
	int *spiv = (int *)malloc((size_t)np * sizeof(int));

	// 電流・電荷分布 (distribution = 1 のときのみ)
	if (p->dist) {
		if (p->nseg > 0) {
			p->segi = (d_complex_t *)malloc((size_t)p->nfreq * p->nseg * sizeof(d_complex_t));
		}
		if (p->capacitance && (p->ncell > 0)) {
			p->cellq = (d_complex_t *)malloc((size_t)p->nfreq * p->ncell * sizeof(d_complex_t));
		}
	}

	if ((a == NULL) || (b == NULL) || (piv == NULL) || (p->zin == NULL) ||
	    (p->zmat == NULL) || (p->smat == NULL) ||
	    (sm == NULL) || (sinv == NULL) || (scol == NULL) || (spiv == NULL)) {
		printf("%s\n", "*** memory allocation error (matrix)");
		free(a); free(b); free(piv);
		free(sm); free(sinv); free(scol); free(spiv);
		return 1;
	}

	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		const double f = freq_at(p, ifreq);

		// 遅延ありのときは部分要素が周波数依存になるので毎回作り直す
		if (p->retardation) {
			lp_fill(p, f, fp_log);
			if (pot_fill(p, f, fp_log)) {
				free(a); free(b); free(piv);
				free(sm); free(sinv); free(scol); free(spiv);
				return 1;
			}
		}

		mna_assemble(p, f, a);
		const int ising = lu_decomp(n, a, piv);
		if (ising >= 0) {
			printf("*** singular matrix at f = %.5e Hz (row %d : floating node or ideal voltage source loop?)\n", f, ising);
			fprintf(fp_log, "*** singular matrix at f = %.5e Hz (row %d)\n", f, ising);
			free(a); free(b); free(piv);
			free(sm); free(sinv); free(scol); free(spiv);
			return 1;
		}

		// ポート j に 1A を注入し、全ポート i の端子電圧を読む。
		// 他のポートには電流を強制していない = 開放条件なので、これがそのまま
		// Z 行列の第 j 列 (Z[i][j] = V_i / I_j, I_j = 1A) になる。
		// 対角が従来の Zin。
		for (int j = 0; j < p->nport; j++) {
			mna_rhs_port(p, j, b);
			lu_solve(n, a, piv, b);
			for (int i = 0; i < p->nport; i++) {
				const int i1 = p->nodemap[p->port[i].n1];
				const int i2 = p->nodemap[p->port[i].n2];
				const d_complex_t v1 = (i1 < 0) ? d_complex(0, 0) : b[i1];
				const d_complex_t v2 = (i2 < 0) ? d_complex(0, 0) : b[i2];
				p->zmat[ZIDX(p, ifreq, i, j)] = d_sub(v1, v2);
			}
			p->zin[(size_t)j * p->nfreq + ifreq] = p->zmat[ZIDX(p, ifreq, j, j)];

			// 分布は port #1 を 1A で励振したときの解から取る
			if (p->dist && (j == 0)) {
				// 区間電流 : 未知数ベクトルの [offS, offS + nseg) がそのまま枝電流
				if (p->segi != NULL) {
					for (int m = 0; m < p->nseg; m++) {
						p->segi[(size_t)ifreq * p->nseg + m] = b[p->offS + m];
					}
				}
				// セル電荷 : q = C v (v は各セルのノード電位、基準ノードは 0V)
				if (p->cellq != NULL) {
					for (int m = 0; m < p->ncell; m++) {
						d_complex_t q = d_complex(0, 0);
						for (int l = 0; l < p->ncell; l++) {
							const int im = p->nodemap[p->cellid[l]];
							const d_complex_t v = (im < 0) ? d_complex(0, 0) : b[im];
							q = d_add(q, d_mul(p->cmat[(size_t)m * p->ncell + l], v));
						}
						p->cellq[(size_t)ifreq * p->ncell + m] = q;
					}
				}
			}
		}

		// Z -> S (基準抵抗はポートごとの z0)
		if (z_to_s(p, ifreq, sm, sinv, scol, spiv)) {
			printf("*** singular (Z + Z0) at f = %.5e Hz\n", f);
			fprintf(fp_log, "*** singular (Z + Z0) at f = %.5e Hz\n", f);
			free(a); free(b); free(piv);
			free(sm); free(sinv); free(scol); free(spiv);
			return 1;
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
	free(sm);
	free(sinv);
	free(scol);
	free(spiv);

	return 0;
}
