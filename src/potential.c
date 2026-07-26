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

retardation = 1 のときは I が複素・周波数依存になるので、周波数ごとに
P を作り直して逆行列を取る。
*/

#include "peec.h"

// 半区間リストと容量セルの構成 (幾何のみ。周波数に依らないので 1 回だけ)
static int build_cells(peec_t *p)
{
	const int nh = 2 * p->nseg;

	p->half = (seg_t *)malloc((size_t)nh * sizeof(seg_t));
	p->cellof = (int *)malloc((size_t)nh * sizeof(int));
	p->cellid = (int *)malloc((size_t)nh * sizeof(int));
	if ((p->half == NULL) || (p->cellof == NULL) || (p->cellid == NULL)) return 1;

	int *hnode = (int *)malloc((size_t)nh * sizeof(int));
	if (hnode == NULL) return 1;

	for (int k = 0; k < p->nseg; k++) {
		const seg_t *s = &p->seg[k];
		double mid[3];
		for (int c = 0; c < 3; c++) {
			mid[c] = 0.5 * (s->x1[c] + s->x2[c]);
		}
		for (int h = 0; h < 2; h++) {
			seg_t *t = &p->half[(2 * k) + h];
			*t = *s;
			for (int c = 0; c < 3; c++) {
				t->x1[c] = (h == 0) ? s->x1[c] : mid[c];
				t->x2[c] = (h == 0) ? mid[c]   : s->x2[c];
			}
			t->len = 0.5 * s->len;
			hnode[(2 * k) + h] = (h == 0) ? s->n1 : s->n2;
		}
	}

	// セル (幾何ノード) の列挙 : 出現順に番号を振る (決定的)
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
		p->cellof[h] = found;
	}

	// セル長
	p->clen = (double *)calloc((size_t)p->ncell, sizeof(double));
	if (p->clen == NULL) return 1;
	for (int h = 0; h < nh; h++) {
		p->clen[p->cellof[h]] += p->half[h].len;
	}

	free(hnode);

	return 0;
}

int pot_fill(peec_t *p, double f, FILE *fp_log)
{
	const char errmem[] = "*** memory allocation error (potential)";

	if (!p->capacitance || (p->nseg <= 0)) return 0;

	if (p->half == NULL) {
		if (build_cells(p)) {
			printf("%s\n", errmem);
			return 1;
		}
	}

	const int nh = 2 * p->nseg;
	const int nc = p->ncell;
	const double kw = p->retardation ? (2 * PI * f / C0) : 0;

	// 半区間どうしの幾何二重積分 (対称)
	d_complex_t *ig = (d_complex_t *)malloc((size_t)nh * nh * sizeof(d_complex_t));
	d_complex_t *pmat = (d_complex_t *)calloc((size_t)nc * nc, sizeof(d_complex_t));
	if ((ig == NULL) || (pmat == NULL)) {
		printf("%s\n", errmem);
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
				? neumann_self_k(&p->half[h1], p->half[h1].aP, kw)
				: neumann_pair_k(&p->half[h1], &p->half[h2],
				                 p->half[h1].aP, p->half[h2].aP, kw);
		}
	}
	for (h1 = 0; h1 < nh; h1++) {
		for (int h2 = 0; h2 < h1; h2++) {
			ig[(size_t)h1 * nh + h2] = ig[(size_t)h2 * nh + h1];
		}
	}

	// 電位係数行列 P (セル単位に集約してからセル長で規格化)
	for (h1 = 0; h1 < nh; h1++) {
		for (int h2 = 0; h2 < nh; h2++) {
			const size_t id = (size_t)p->cellof[h1] * nc + p->cellof[h2];
			pmat[id] = d_add(pmat[id], ig[(size_t)h1 * nh + h2]);
		}
	}
	for (int i = 0; i < nc; i++) {
		for (int j = 0; j < nc; j++) {
			pmat[(size_t)i * nc + j] = d_rmul(
				1 / (4 * PI * EPS0 * p->clen[i] * p->clen[j]), pmat[(size_t)i * nc + j]);
		}
	}

	// C = P^-1
	d_complex_t *b = (d_complex_t *)malloc((size_t)nc * sizeof(d_complex_t));
	int *piv = (int *)malloc((size_t)nc * sizeof(int));
	if (p->cmat == NULL) {
		p->cmat = (d_complex_t *)malloc((size_t)nc * nc * sizeof(d_complex_t));
	}
	if ((b == NULL) || (piv == NULL) || (p->cmat == NULL)) {
		printf("%s\n", errmem);
		return 1;
	}
	const int ising = lu_decomp(nc, pmat, piv);
	if (ising >= 0) {
		printf("*** singular potential coefficient matrix (row %d)\n", ising);
		fprintf(fp_log, "*** singular potential coefficient matrix (row %d)\n", ising);
		return 1;
	}
	for (int j = 0; j < nc; j++) {
		memset(b, 0, (size_t)nc * sizeof(d_complex_t));
		b[j] = d_complex(1, 0);
		lu_solve(nc, pmat, piv, b);
		for (int i = 0; i < nc; i++) {
			p->cmat[(size_t)i * nc + j] = b[i];
		}
	}

	// 対無限遠の総容量 (準静的なら実数。ログは初回のみ)
	if (!p->clogged) {
		p->clogged = 1;
		double ctotal = 0;
		for (int i = 0; i < nc; i++) {
			for (int j = 0; j < nc; j++) {
				ctotal += p->cmat[(size_t)i * nc + j].r;
			}
		}
		printf("PEEC: capacitive cells = %d, total capacitance = %.6e F\n", nc, ctotal);
		fprintf(fp_log, "PEEC: capacitive cells = %d, total capacitance = %.6e F\n", nc, ctotal);
		fflush(fp_log);
	}

	free(ig);
	free(pmat);
	free(b);
	free(piv);

	return 0;
}
