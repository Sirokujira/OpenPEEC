/*
potential.c

容量性 PEEC (電位係数 P と節点容量行列 C = P^-1)

容量セルは幾何ノードに対応させる (標準的な thin-wire PEEC)。
各区間はその半分ずつを両端ノードのセルに与えるので、セル n の電荷は
そのノードに接する半区間上に一様に分布するものとして扱う。

  P(n,m) = 1 / (4 pi eps0 Ln Lm) * sum_{p in n} sum_{q in m} I(p,q)
  C      = P^-1                              … Maxwell の容量行列
  Ctotal = sum_ij C(i,j)                     … 対無限遠の総容量

MNA では節点ブロックに j omega C を加える (基準ノードの行・列は落とす =
基準ノードを無限遠の電位基準に固定することに相当する)。
*/

#include "peec.h"

int pot_fill(peec_t *p, FILE *fp_log)
{
	if (!p->capacitance || (p->nseg <= 0)) return 0;

	// 半区間リスト (区間 k -> 2k : n1 側、2k+1 : n2 側)
	const int nh = 2 * p->nseg;
	seg_t *half = (seg_t *)malloc((size_t)nh * sizeof(seg_t));
	int *hnode = (int *)malloc((size_t)nh * sizeof(int));
	if ((half == NULL) || (hnode == NULL)) {
		printf("%s\n", "*** memory allocation error (potential)");
		return 1;
	}

	for (int k = 0; k < p->nseg; k++) {
		const seg_t *s = &p->seg[k];
		double mid[3];
		for (int c = 0; c < 3; c++) {
			mid[c] = 0.5 * (s->x1[c] + s->x2[c]);
		}
		for (int h = 0; h < 2; h++) {
			seg_t *t = &half[(2 * k) + h];
			for (int c = 0; c < 3; c++) {
				t->x1[c] = (h == 0) ? s->x1[c] : mid[c];
				t->x2[c] = (h == 0) ? mid[c]   : s->x2[c];
			}
			t->len = 0.5 * s->len;
			t->radius = s->radius;
			t->sigma = s->sigma;
			t->res = 0;
			t->n1 = t->n2 = 0;
			hnode[(2 * k) + h] = (h == 0) ? s->n1 : s->n2;
		}
	}

	// セル (幾何ノード) の列挙 : 出現順に番号を振る (決定的)
	p->cellid = (int *)malloc((size_t)nh * sizeof(int));
	int *cellof = (int *)malloc((size_t)nh * sizeof(int));
	if ((p->cellid == NULL) || (cellof == NULL)) {
		printf("%s\n", "*** memory allocation error (potential)");
		return 1;
	}
	p->ncell = 0;
	for (int h = 0; h < nh; h++) {
		int found = -1;
		for (int i = 0; i < p->ncell; i++) {
			if (p->cellid[i] == hnode[h]) {
				found = i;
				break;
			}
		}
		if (found < 0) {
			found = p->ncell;
			p->cellid[p->ncell++] = hnode[h];
		}
		cellof[h] = found;
	}
	const int nc = p->ncell;

	// セル長
	double *clen = (double *)calloc((size_t)nc, sizeof(double));
	for (int h = 0; h < nh; h++) {
		clen[cellof[h]] += half[h].len;
	}

	// 半区間どうしの幾何二重積分 (対称)
	double *ig = (double *)malloc((size_t)nh * nh * sizeof(double));
	if ((clen == NULL) || (ig == NULL)) {
		printf("%s\n", "*** memory allocation error (potential)");
		return 1;
	}

	// MSVC の OpenMP 2.0 は for 文内でのインデックス宣言を許さない (C3015)
	int h1;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (h1 = 0; h1 < nh; h1++) {
		for (int h2 = h1; h2 < nh; h2++) {
			ig[(size_t)h1 * nh + h2] = (h1 == h2)
				? neumann_self(half[h1].len, half[h1].radius)
				: neumann_pair(&half[h1], &half[h2]);
		}
	}
	for (h1 = 0; h1 < nh; h1++) {
		for (int h2 = 0; h2 < h1; h2++) {
			ig[(size_t)h1 * nh + h2] = ig[(size_t)h2 * nh + h1];
		}
	}

	// 電位係数行列 P (セル単位に集約してからセル長で規格化)
	double *pmat = (double *)calloc((size_t)nc * nc, sizeof(double));
	if (pmat == NULL) {
		printf("%s\n", "*** memory allocation error (potential)");
		return 1;
	}
	for (h1 = 0; h1 < nh; h1++) {
		for (int h2 = 0; h2 < nh; h2++) {
			pmat[(size_t)cellof[h1] * nc + cellof[h2]] += ig[(size_t)h1 * nh + h2];
		}
	}
	for (int i = 0; i < nc; i++) {
		for (int j = 0; j < nc; j++) {
			pmat[(size_t)i * nc + j] /= 4 * PI * EPS0 * clen[i] * clen[j];
		}
	}

	// C = P^-1 (複素 LU を流用 : 虚部 0)
	d_complex_t *am = (d_complex_t *)malloc((size_t)nc * nc * sizeof(d_complex_t));
	d_complex_t *b = (d_complex_t *)malloc((size_t)nc * sizeof(d_complex_t));
	int *piv = (int *)malloc((size_t)nc * sizeof(int));
	p->cmat = (double *)malloc((size_t)nc * nc * sizeof(double));
	if ((am == NULL) || (b == NULL) || (piv == NULL) || (p->cmat == NULL)) {
		printf("%s\n", "*** memory allocation error (potential)");
		return 1;
	}
	for (int i = 0; i < nc; i++) {
		for (int j = 0; j < nc; j++) {
			am[(size_t)i * nc + j] = d_complex(pmat[(size_t)i * nc + j], 0);
		}
	}
	const int ising = lu_decomp(nc, am, piv);
	if (ising >= 0) {
		printf("*** singular potential coefficient matrix (row %d)\n", ising);
		fprintf(fp_log, "*** singular potential coefficient matrix (row %d)\n", ising);
		return 1;
	}
	for (int j = 0; j < nc; j++) {
		memset(b, 0, (size_t)nc * sizeof(d_complex_t));
		b[j] = d_complex(1, 0);
		lu_solve(nc, am, piv, b);
		for (int i = 0; i < nc; i++) {
			p->cmat[(size_t)i * nc + j] = b[i].r;
		}
	}

	// 対無限遠の総容量
	p->ctotal = 0;
	for (int i = 0; i < nc; i++) {
		for (int j = 0; j < nc; j++) {
			p->ctotal += p->cmat[(size_t)i * nc + j];
		}
	}

	printf("PEEC: capacitive cells = %d, total capacitance = %.6e F\n", nc, p->ctotal);
	fprintf(fp_log, "PEEC: capacitive cells = %d, total capacitance = %.6e F\n", nc, p->ctotal);
	fflush(fp_log);

	free(half);
	free(hnode);
	free(cellof);
	free(clen);
	free(ig);
	free(pmat);
	free(am);
	free(b);
	free(piv);

	return 0;
}
