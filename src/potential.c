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

/*
セル -> 電荷サブセルの索引 (圧縮経路で 1 要素ずつ評価するのに要る)

密経路は nchg x nchg の ig を作ってから集約するが、圧縮経路では
P(n, m) を単独で求めたいので、セルに属するサブセルを引けるようにする。
サブセル数はノードに集まる区間数 (線導体) か 1 (面導体の双対矩形) なので
数個で、1 要素の評価コストは数回の二重積分に収まる。
*/
static int build_cellsub(peec_t *p)
{
	const int nh = p->nchg;
	const int nc = p->ncell;

	p->csoff = (int *)calloc((size_t)(nc + 1), sizeof(int));
	p->csidx = (int *)malloc((size_t)(nh > 0 ? nh : 1) * sizeof(int));
	p->cpt = (double *)malloc((size_t)(nc > 0 ? nc : 1) * 3 * sizeof(double));
	if ((p->csoff == NULL) || (p->csidx == NULL) || (p->cpt == NULL)) return 1;

	for (int h = 0; h < nh; h++) {
		p->csoff[p->cellof[h] + 1]++;
	}
	for (int i = 0; i < nc; i++) {
		p->csoff[i + 1] += p->csoff[i];
	}
	{
		int *fill = (int *)malloc((size_t)(nc + 1) * sizeof(int));
		if (fill == NULL) return 1;
		memcpy(fill, p->csoff, (size_t)(nc + 1) * sizeof(int));
		for (int h = 0; h < nh; h++) {
			p->csidx[fill[p->cellof[h]]++] = h;
		}
		free(fill);
	}

	// セルの代表点 : サブセル中点を長さで重み付けした重心 (クラスタツリー用)
	for (int i = 0; i < nc; i++) {
		double w = 0, x[3] = {0, 0, 0};
		for (int e = p->csoff[i]; e < p->csoff[i + 1]; e++) {
			const seg_t *s = &p->chg[p->csidx[e]];
			for (int c = 0; c < 3; c++) {
				x[c] += s->len * 0.5 * (s->x1[c] + s->x2[c]);
			}
			w += s->len;
		}
		for (int c = 0; c < 3; c++) {
			p->cpt[(3 * i) + c] = (w > 0) ? (x[c] / w) : 0;
		}
	}

	return 0;
}

/*
電位係数行列の 1 要素 P(i, j)

  P(i, j) = 1/(4 pi eps0 Ai Aj) * sum_{p in i} sum_{q in j} I(p, q)

pot_fill() の集約と同じ式で、地板の鏡像減算も同じ規約 (不変条件 8)。
密経路と圧縮経路が同じ関数を通るので構造的にビット一致する。
*/
d_complex_t pot_entry(const peec_t *p, int i, int j, double kw)
{
	d_complex_t s = d_complex(0, 0);

	for (int e = p->csoff[i]; e < p->csoff[i + 1]; e++) {
		const int h1 = p->csidx[e];
		for (int g = p->csoff[j]; g < p->csoff[j + 1]; g++) {
			const int h2 = p->csidx[g];
			// 密経路は上三角を評価して写すので、ここでも順序を正規化する
			const int a = (h1 <= h2) ? h1 : h2;
			const int b = (h1 <= h2) ? h2 : h1;
			d_complex_t v = (a == b)
				? neumann_self_k(&p->chg[a], p->chg[a].aP, kw)
				: neumann_pair_k(&p->chg[a], &p->chg[b],
				                 p->chg[a].aP, p->chg[b].aP, kw);
			if (p->gp) {
				seg_t m2;
				seg_mirror(&p->chg[b], p->gpz, &m2);
				v = d_sub(v, neumann_pair_k(&p->chg[a], &m2,
					p->chg[a].aP, p->chg[b].aP, kw));
			}
			s = d_add(s, v);
		}
	}

	return d_rmul(1 / (4 * PI * EPS0 * p->carea[i] * p->carea[j]), s);
}

static d_complex_t pot_entry_cb(void *ctx, int i, int j, double kw)
{
	return pot_entry((const peec_t *)ctx, i, j, kw);
}

// 電位係数 P を H 行列で圧縮する (compression = 1 かつ capacitance = 1)
struct hmat_t *hmat_build_pot(peec_t *p, double f, FILE *fp_log)
{
	if (pot_cells(p)) return NULL;
	if (p->ncell <= 0) return NULL;
	const double kw = p->retardation ? (2 * PI * f / C0) : 0;

	return hmat_build_gen(p->ncell, p->cpt, pot_entry_cb, p, kw, p->ctol, "P", fp_log);
}

// 容量セルの構成 (幾何のみなので 1 回だけ)。密・圧縮の両経路から呼ぶ。
int pot_cells(peec_t *p)
{
	if (p->cellof != NULL) return 0;
	if (p->nchg <= 0) return 1;
	if (build_cells(p)) return 1;

	return build_cellsub(p);
}

int pot_fill(peec_t *p, double f, FILE *fp_log)
{
	const char errmem[] = "*** memory allocation error (potential)";

	if (!p->capacitance || (p->nchg <= 0)) return 0;

	if (pot_cells(p)) {
		printf("%s\n", errmem);
		return 1;
	}

	/*
	圧縮経路は密な P も C = P^-1 も作らない (電荷が未知数で、P は H 行列)。
	総容量 Ctotal = 1^T P^-1 1 は逆行列を要するのでここでは出さない
	(必要なら compression を外して求める)。
	*/
	if (p->compress) {
		if (!p->clogged) {
			p->clogged = 1;
			printf("PEEC: capacitive cells = %d (compressed : charge unknowns)\n", p->ncell);
			fprintf(fp_log, "PEEC: capacitive cells = %d (compressed : charge unknowns)\n",
				p->ncell);
			fflush(fp_log);
		}
		return 0;
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
