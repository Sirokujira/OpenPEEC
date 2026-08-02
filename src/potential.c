/*
potential.c

容量性 PEEC (電位係数 P と節点容量行列 C = P^-1)

電荷セルは幾何段 (wire.c) が作る :
- 線導体 : 各区間の半分ずつを両端ノードに与えた半区間
- 面導体 : 各格子ノードの双対矩形 (端は半分、隅は 1/4)
これらを同じノードごとにまとめたものが容量セルになる。

  P(n,m) = 1 / (4 pi eps0 An Am) * sum_{p in n} sum_{q in m} I(p,q)
  C      = P^-1                              … Maxwell の容量行列
  Ctotal = sum_ij C(i,j)                     … 対無限遠の総容量

ここで I は幾何二重積分 (細線なら ∬dl dl'/R、面なら幅で規格化した
∬∬dS dS'/(w1 w2 R))、An はセルの長さ (細線) / 面積を幅で割った量。
どちらも同じ式で扱えるよう、セルの「長さ相当量」 len を足し上げる。

MNA では節点ブロックに j omega C を加える (基準ノードの行・列は落とす =
基準ノードを無限遠の電位基準に固定することに相当する)。

retardation = 1 のときは I が複素・周波数依存になるので、周波数ごとに
P を作り直して逆行列を取る。
*/

#include "peec.h"

// 容量セル (ノードごとの集約) の構成。幾何のみなので 1 回だけ。
static int build_cells(peec_t *p)
{
	const int nh = p->nchg;

	p->cellof = (int *)malloc((size_t)nh * sizeof(int));
	p->cellid = (int *)malloc((size_t)nh * sizeof(int));
	if ((p->cellof == NULL) || (p->cellid == NULL)) return 1;

	// セル (幾何ノード) の列挙 : 出現順に番号を振る (決定的)
	p->ncell = 0;
	for (int h = 0; h < nh; h++) {
		int found = -1;
		for (int i = 0; i < p->ncell; i++) {
			if (p->cellid[i] == p->chgnode[h]) {
				found = i;
				break;
			}
		}
		if (found < 0) {
			found = p->ncell;
			p->cellid[p->ncell++] = p->chgnode[h];
		}
		p->cellof[h] = found;
	}

	// セルの長さ相当量 (面導体では面積/幅ではなく長さを足す : 幅は I 側で規格化済み)
	p->carea = (double *)calloc((size_t)p->ncell, sizeof(double));
	if (p->carea == NULL) return 1;
	for (int h = 0; h < nh; h++) {
		p->carea[p->cellof[h]] += p->chg[h].len;
	}

	return 0;
}

int pot_fill(peec_t *p, double f, FILE *fp_log)
{
	const char errmem[] = "*** memory allocation error (potential)";

	if (!p->capacitance || (p->nchg <= 0)) return 0;

	if (p->cellof == NULL) {
		if (build_cells(p)) {
			printf("%s\n", errmem);
			return 1;
		}
	}

	const int nh = p->nchg;
	const int nc = p->ncell;
	const double kw = p->retardation ? (2 * PI * f / C0) : 0;

	// 電荷セルどうしの幾何二重積分 (対称)
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
			d_complex_t v = (h1 == h2)
				? neumann_self_k(&p->chg[h1], p->chg[h1].aP, kw)
				: neumann_pair_k(&p->chg[h1], &p->chg[h2],
				                 p->chg[h1].aP, p->chg[h2].aP, kw);
			if (p->gp) {
				// 地板 : 鏡像電荷は -q なので鏡像との相互項を減算する
				// (自己項 h1 == h2 も自分の鏡像との相互項を引く)
				seg_t m2;
				seg_mirror(&p->chg[h2], p->gpz, &m2);
				v = d_sub(v, neumann_pair_k(&p->chg[h1], &m2,
					p->chg[h1].aP, p->chg[h2].aP, kw));
			}
			ig[(size_t)h1 * nh + h2] = v;
		}
	}
	for (h1 = 0; h1 < nh; h1++) {
		for (int h2 = 0; h2 < h1; h2++) {
			ig[(size_t)h1 * nh + h2] = ig[(size_t)h2 * nh + h1];
		}
	}

	// 電位係数行列 P (セル単位に集約してから規格化)
	for (h1 = 0; h1 < nh; h1++) {
		for (int h2 = 0; h2 < nh; h2++) {
			const size_t id = (size_t)p->cellof[h1] * nc + p->cellof[h2];
			pmat[id] = d_add(pmat[id], ig[(size_t)h1 * nh + h2]);
		}
	}
	for (int i = 0; i < nc; i++) {
		for (int j = 0; j < nc; j++) {
			pmat[(size_t)i * nc + j] = d_rmul(
				1 / (4 * PI * EPS0 * p->carea[i] * p->carea[j]), pmat[(size_t)i * nc + j]);
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
