/*
precond.c

行列フリー経路 (compression = 1) の前処理 : 葉ブロック消去 + 回路 Schur 補元

MNA の未知数は [回路 (節点・i_L・i_V) | 圧縮ブロック…] に分かれる。
圧縮ブロックは区間電流 (カーネル Lp) と、容量性 PEEC のときは電荷
(カーネル P) の 2 つで、いずれも回路ブロックとは疎にしか結合しない :

  A = | S11  B_1  B_2 |      n1 = offS (回路)
      | C_1  D_1   0  |      D_b = 圧縮ブロック b の密なカーネル + 対角
      | C_2   0   D_2 |      (B_b / C_b は接続行列 = 1 行 2 要素)

D_b は密だが、前処理では**クラスタツリーの葉の対角ブロックだけ**に落とす
(実測で密行列の 1.6〜3.1%)。するとブロック対角になって厳密に消去でき、

  M^-1 は  D_b の葉ごとの LU (m <= HM_LEAF) と
  Schur 補元 S = S11 - sum_b B_b D_b^-1 C_b (n1 x n1 密) の LU

だけで適用できる。適用は前進消去 → Schur 解 → 後退代入の直列 3 段。
異なる圧縮ブロックどうしは直接結合しない (電流と電荷は回路ブロックを
介してのみつながる) ので、消去は互いに独立に足し込める。

設計の根拠 (プロトタイプでの実測、密 LU を前処理に使った GMRES の反復数) :
  - 近傍ブロック全部     : 18〜29 反復 (Lp の 9.6〜32%)
  - 葉の対角ブロックのみ : 24 反復     (Lp の 1.6〜3.1%)   <- 採用
  - Lp の対角のみ        : 30〜36 反復
  - 近傍パターンの ILU(0): 収束せず (棄却)
  - 近傍場の厳密疎 LU    : fill が n^2 の 75% (密と変わらないので棄却)

メモリは葉 LU が ~HM_LEAF x (未知数)、Schur が n1^2。Schur の n1^2 は
O(N^2) のまま残るが、部分要素の充填 O(N^2) を H 行列が O(kN log N) に
落とすのが主目的なので、この段階ではこれでよい。

スレッド数不変性 : 葉ブロックの LU は互いに独立なので並列化してよい。
Schur の集約と適用は直列 (GMRES 自体の内積が直列なのでボトルネックではない)。
*/

#include "peec.h"

#define PC_MAXM  64                   // 葉の最大セル数 (HM_LEAF) の上限
#define PC_MAXT  (4 * PC_MAXM)        // 1 葉が触れる回路未知数の上限

// 圧縮ブロック 1 個ぶんの消去データ
typedef struct {
	int off, n;                   // 未知数の先頭と個数
	int nleaf;
	int    *lof;                  // [nleaf+1] lidx へのオフセット
	int    *lidx;                 // [n] 葉ごとの (ブロック内) 番号
	size_t *dof;                  // [nleaf+1] dlu へのオフセット (m^2)
	d_complex_t *dlu;
	int    *dpiv;                 // [n]
	// 疎な結合 : B (回路行 x ブロック列)、C (ブロック行 x 回路列)
	int nb, nc;
	int *bri, *bci; d_complex_t *bval;
	int *cri, *cci; d_complex_t *cval;
	int *boff, *bidx;             // ブロック列ごとの B 要素 (CSR 風)
	int *coff, *cidx;             // ブロック行ごとの C 要素
} pcblk_t;

struct precond_t {
	int n, n1;
	int nblk;
	pcblk_t blk[HBLK_MAX];
	d_complex_t *slu;             // Schur 補元 (n1 x n1) の LU
	int *spiv;
	d_complex_t *w1;              // [n1]
	d_complex_t *wb;              // [max n] 作業領域
};

static void pcblk_free(pcblk_t *b)
{
	free(b->lof); free(b->lidx); free(b->dof); free(b->dlu); free(b->dpiv);
	free(b->bri); free(b->bci); free(b->bval);
	free(b->cri); free(b->cci); free(b->cval);
	free(b->boff); free(b->bidx); free(b->coff); free(b->cidx);
	memset(b, 0, sizeof(pcblk_t));
}

void pc_free(struct precond_t *pc)
{
	if (pc == NULL) return;
	for (int b = 0; b < pc->nblk; b++) {
		pcblk_free(&pc->blk[b]);
	}
	free(pc->slu); free(pc->spiv); free(pc->w1); free(pc->wb);
	free(pc);
}

// 葉の取り込みと D_b の初期化 (カーネル x scale)。戻り値 0 = 正常
static int blk_leaves(pcblk_t *b, const hblock_t *src)
{
	b->off = src->off;
	b->n = src->n;
	b->nleaf = hmat_nleaf(src->h);
	b->lof  = (int *)malloc((size_t)(b->nleaf + 1) * sizeof(int));
	b->lidx = (int *)malloc((size_t)(b->n > 0 ? b->n : 1) * sizeof(int));
	b->dof  = (size_t *)malloc((size_t)(b->nleaf + 1) * sizeof(size_t));
	b->dpiv = (int *)malloc((size_t)(b->n > 0 ? b->n : 1) * sizeof(int));
	if ((b->lof == NULL) || (b->lidx == NULL) || (b->dof == NULL)
	 || (b->dpiv == NULL)) return 1;

	b->lof[0] = 0;
	b->dof[0] = 0;
	for (int r = 0; r < b->nleaf; r++) {
		const int *idx;
		int m;
		if (hmat_leaf_block(src->h, r, &idx, &m) == NULL) return 1;
		if (m > PC_MAXM) return 1;
		memcpy(&b->lidx[b->lof[r]], idx, (size_t)m * sizeof(int));
		b->lof[r + 1] = b->lof[r] + m;
		b->dof[r + 1] = b->dof[r] + ((size_t)m * m);
	}
	b->dlu = (d_complex_t *)malloc((b->dof[b->nleaf] > 0 ? b->dof[b->nleaf] : 1)
		* sizeof(d_complex_t));
	if (b->dlu == NULL) return 1;

	// D_b = scale * (葉の対角カーネルブロック)
	int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (r = 0; r < b->nleaf; r++) {
		const int *idx;
		int m;
		const d_complex_t *lb = hmat_leaf_block(src->h, r, &idx, &m);
		d_complex_t *d = &b->dlu[b->dof[r]];
		for (int a = 0; a < m; a++) {
			for (int c = 0; c < m; c++) {
				d[(size_t)a * m + c] = d_mul(src->scale, lb[(size_t)a * m + c]);
			}
		}
	}

	return 0;
}

struct precond_t *pc_build(const peec_t *p, const mna_sparse_t *sp,
	const hblock_t *blk, int nblk)
{
	const int n1 = p->offS;

	struct precond_t *pc = (struct precond_t *)calloc(1, sizeof(struct precond_t));
	if (pc == NULL) return NULL;
	pc->n = p->nunknown;
	pc->n1 = n1;
	pc->nblk = nblk;

	int maxn = 1;
	for (int b = 0; b < nblk; b++) {
		if (blk_leaves(&pc->blk[b], &blk[b])) {
			pc_free(pc);
			return NULL;
		}
		if (blk[b].n > maxn) maxn = blk[b].n;
	}

	// 未知数 -> (ブロック, ブロック内番号)
	int *ublk = (int *)malloc((size_t)pc->n * sizeof(int));
	int *upos = (int *)malloc((size_t)pc->n * sizeof(int));
	if ((ublk == NULL) || (upos == NULL)) {
		free(ublk); free(upos);
		pc_free(pc);
		return NULL;
	}
	for (int i = 0; i < pc->n; i++) {
		ublk[i] = -1;
		upos[i] = -1;
	}
	for (int b = 0; b < nblk; b++) {
		for (int i = 0; i < blk[b].n; i++) {
			ublk[blk[b].off + i] = b;
			upos[blk[b].off + i] = i;
		}
	}
	// ブロック内番号 -> (葉, 葉内位置)
	int *sleaf[HBLK_MAX], *spos[HBLK_MAX];
	for (int b = 0; b < nblk; b++) {
		pcblk_t *pb = &pc->blk[b];
		sleaf[b] = (int *)malloc((size_t)(pb->n > 0 ? pb->n : 1) * sizeof(int));
		spos[b]  = (int *)malloc((size_t)(pb->n > 0 ? pb->n : 1) * sizeof(int));
		if ((sleaf[b] == NULL) || (spos[b] == NULL)) {
			for (int c = 0; c <= b; c++) { free(sleaf[c]); free(spos[c]); }
			free(ublk); free(upos);
			pc_free(pc);
			return NULL;
		}
		for (int r = 0; r < pb->nleaf; r++) {
			for (int e = pb->lof[r]; e < pb->lof[r + 1]; e++) {
				sleaf[b][pb->lidx[e]] = r;
				spos[b][pb->lidx[e]] = e - pb->lof[r];
			}
		}
	}

	// ── 疎部の振り分け : S11 / B_b / C_b / D_b への追加分 ──────────
	pc->slu = (d_complex_t *)calloc((size_t)(n1 > 0 ? n1 : 1) * (n1 > 0 ? n1 : 1),
		sizeof(d_complex_t));
	pc->spiv = (int *)malloc((size_t)(n1 > 0 ? n1 : 1) * sizeof(int));
	if ((pc->slu == NULL) || (pc->spiv == NULL)) goto fail;

	for (int b = 0; b < nblk; b++) {
		pc->blk[b].nb = 0;
		pc->blk[b].nc = 0;
	}
	for (int e = 0; e < sp->nnz; e++) {
		const int i = sp->ri[e], j = sp->ci[e];
		if ((i < n1) && (j >= n1)) pc->blk[ublk[j]].nb++;
		else if ((i >= n1) && (j < n1)) pc->blk[ublk[i]].nc++;
	}
	for (int b = 0; b < nblk; b++) {
		pcblk_t *pb = &pc->blk[b];
		pb->bri = (int *)malloc((size_t)(pb->nb > 0 ? pb->nb : 1) * sizeof(int));
		pb->bci = (int *)malloc((size_t)(pb->nb > 0 ? pb->nb : 1) * sizeof(int));
		pb->bval = (d_complex_t *)malloc((size_t)(pb->nb > 0 ? pb->nb : 1) * sizeof(d_complex_t));
		pb->cri = (int *)malloc((size_t)(pb->nc > 0 ? pb->nc : 1) * sizeof(int));
		pb->cci = (int *)malloc((size_t)(pb->nc > 0 ? pb->nc : 1) * sizeof(int));
		pb->cval = (d_complex_t *)malloc((size_t)(pb->nc > 0 ? pb->nc : 1) * sizeof(d_complex_t));
		if ((pb->bri == NULL) || (pb->bci == NULL) || (pb->bval == NULL)
		 || (pb->cri == NULL) || (pb->cci == NULL) || (pb->cval == NULL)) goto fail;
		pb->nb = 0;
		pb->nc = 0;
	}

	for (int e = 0; e < sp->nnz; e++) {
		const int i = sp->ri[e], j = sp->ci[e];
		if ((i < n1) && (j < n1)) {
			pc->slu[(size_t)i * n1 + j] = d_add(pc->slu[(size_t)i * n1 + j], sp->val[e]);
		}
		else if (i < n1) {
			pcblk_t *pb = &pc->blk[ublk[j]];
			pb->bri[pb->nb] = i;
			pb->bci[pb->nb] = upos[j];
			pb->bval[pb->nb] = sp->val[e];
			pb->nb++;
		}
		else if (j < n1) {
			pcblk_t *pb = &pc->blk[ublk[i]];
			pb->cri[pb->nc] = upos[i];
			pb->cci[pb->nc] = j;
			pb->cval[pb->nc] = sp->val[e];
			pb->nc++;
		}
		else if (ublk[i] == ublk[j]) {
			// 同じブロック内の疎な要素 (Z_int の対角など)。同じ葉なら D に足す
			const int b = ublk[i];
			const int k = upos[i], l = upos[j];
			if (sleaf[b][k] == sleaf[b][l]) {
				pcblk_t *pb = &pc->blk[b];
				const int lr = sleaf[b][k];
				const int m = pb->lof[lr + 1] - pb->lof[lr];
				const size_t id = pb->dof[lr] + ((size_t)spos[b][k] * m) + spos[b][l];
				pb->dlu[id] = d_add(pb->dlu[id], sp->val[e]);
			}
		}
	}

	// ── 葉ブロックの LU (互いに独立 = 並列化してよい) ────────────
	for (int b = 0; b < nblk; b++) {
		pcblk_t *pb = &pc->blk[b];
		int fail = 0;
		int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
		for (r = 0; r < pb->nleaf; r++) {
			const int m = pb->lof[r + 1] - pb->lof[r];
			if (lu_decomp(m, &pb->dlu[pb->dof[r]], &pb->dpiv[pb->lof[r]]) >= 0) {
				fail = 1;
			}
		}
		if (fail) goto fail;
	}

	// ── Schur 補元 S = S11 - sum_b B_b D_b^-1 C_b ─────────────────
	// 葉 r ごとに : 触れる回路未知数の集合 t、C_r (m x t) を D_r で解き、
	// B の該当列から S のクリークに引く。集約は直列 (共有ノードの競合回避)。
	{
		int *tmap = (int *)malloc((size_t)(n1 > 0 ? n1 : 1) * sizeof(int));
		int *tset = (int *)malloc((size_t)PC_MAXT * sizeof(int));
		d_complex_t *tblk = (d_complex_t *)malloc((size_t)PC_MAXM * PC_MAXT
			* sizeof(d_complex_t));
		if ((tmap == NULL) || (tset == NULL) || (tblk == NULL)) {
			free(tmap); free(tset); free(tblk);
			goto fail;
		}
		for (int j = 0; j < n1; j++) tmap[j] = -1;

		int ok = 1;
		for (int b = 0; (b < nblk) && ok; b++) {
			pcblk_t *pb = &pc->blk[b];
			// ブロック行ごとの C 要素、ブロック列ごとの B 要素 (CSR 風)
			pb->coff = (int *)calloc((size_t)(pb->n + 1), sizeof(int));
			pb->cidx = (int *)malloc((size_t)(pb->nc > 0 ? pb->nc : 1) * sizeof(int));
			pb->boff = (int *)calloc((size_t)(pb->n + 1), sizeof(int));
			pb->bidx = (int *)malloc((size_t)(pb->nb > 0 ? pb->nb : 1) * sizeof(int));
			if ((pb->coff == NULL) || (pb->cidx == NULL)
			 || (pb->boff == NULL) || (pb->bidx == NULL)) {
				ok = 0;
				break;
			}
			for (int e = 0; e < pb->nc; e++) pb->coff[pb->cri[e] + 1]++;
			for (int k = 0; k < pb->n; k++) pb->coff[k + 1] += pb->coff[k];
			for (int e = 0; e < pb->nc; e++) pb->cidx[pb->coff[pb->cri[e]]++] = e;
			for (int k = pb->n; k > 0; k--) pb->coff[k] = pb->coff[k - 1];
			pb->coff[0] = 0;
			for (int e = 0; e < pb->nb; e++) pb->boff[pb->bci[e] + 1]++;
			for (int k = 0; k < pb->n; k++) pb->boff[k + 1] += pb->boff[k];
			for (int e = 0; e < pb->nb; e++) pb->bidx[pb->boff[pb->bci[e]]++] = e;
			for (int k = pb->n; k > 0; k--) pb->boff[k] = pb->boff[k - 1];
			pb->boff[0] = 0;

			for (int r = 0; r < pb->nleaf; r++) {
				const int m = pb->lof[r + 1] - pb->lof[r];
				const int *idx = &pb->lidx[pb->lof[r]];
				// 触れる回路未知数の集合 (順序は走査順で決定的)
				int t = 0;
				for (int a = 0; (a < m) && ok; a++) {
					for (int q = pb->coff[idx[a]]; q < pb->coff[idx[a] + 1]; q++) {
						const int j = pb->cci[pb->cidx[q]];
						if (tmap[j] < 0) {
							if (t >= PC_MAXT) { ok = 0; break; }
							tmap[j] = t;
							tset[t++] = j;
						}
					}
				}
				if (!ok) break;
				if (t == 0) continue;
				// T = C_r (m x t) を組み、D_r^-1 T を列ごとに解く
				memset(tblk, 0, (size_t)m * t * sizeof(d_complex_t));
				for (int a = 0; a < m; a++) {
					for (int q = pb->coff[idx[a]]; q < pb->coff[idx[a] + 1]; q++) {
						const int e = pb->cidx[q];
						const size_t id = (size_t)a * t + tmap[pb->cci[e]];
						tblk[id] = d_add(tblk[id], pb->cval[e]);
					}
				}
				for (int c = 0; c < t; c++) {
					d_complex_t col[PC_MAXM];
					for (int a = 0; a < m; a++) col[a] = tblk[(size_t)a * t + c];
					lu_solve(m, &pb->dlu[pb->dof[r]], &pb->dpiv[pb->lof[r]], col);
					for (int a = 0; a < m; a++) tblk[(size_t)a * t + c] = col[a];
				}
				// S[i][tset[c]] -= B[i][idx[a]] * X[a][c]
				for (int a = 0; a < m; a++) {
					for (int q = pb->boff[idx[a]]; q < pb->boff[idx[a] + 1]; q++) {
						const int e = pb->bidx[q];
						const int i = pb->bri[e];
						for (int c = 0; c < t; c++) {
							const size_t id = (size_t)i * n1 + tset[c];
							pc->slu[id] = d_sub(pc->slu[id],
								d_mul(pb->bval[e], tblk[(size_t)a * t + c]));
						}
					}
				}
				for (int c = 0; c < t; c++) tmap[tset[c]] = -1;
			}
		}
		free(tmap);
		free(tset);
		free(tblk);
		if (!ok) goto fail;
	}

	for (int b = 0; b < nblk; b++) { free(sleaf[b]); free(spos[b]); }
	free(ublk);
	free(upos);

	if ((n1 > 0) && (lu_decomp(n1, pc->slu, pc->spiv) >= 0)) {
		pc_free(pc);
		return NULL;
	}

	pc->w1 = (d_complex_t *)malloc((size_t)(n1 > 0 ? n1 : 1) * sizeof(d_complex_t));
	pc->wb = (d_complex_t *)malloc((size_t)maxn * sizeof(d_complex_t));
	if ((pc->w1 == NULL) || (pc->wb == NULL)) {
		pc_free(pc);
		return NULL;
	}

	return pc;

fail:
	for (int b = 0; b < nblk; b++) { free(sleaf[b]); free(spos[b]); }
	free(ublk);
	free(upos);
	pc_free(pc);

	return NULL;
}

// z = D_b^-1 u (葉ごとに小さな LU で解く)
static void blk_dsolve(const pcblk_t *b, const d_complex_t *u, d_complex_t *z)
{
	for (int r = 0; r < b->nleaf; r++) {
		const int m = b->lof[r + 1] - b->lof[r];
		const int *idx = &b->lidx[b->lof[r]];
		d_complex_t col[PC_MAXM];
		for (int a = 0; a < m; a++) col[a] = u[idx[a]];
		lu_solve(m, &b->dlu[b->dof[r]], &b->dpiv[b->lof[r]], col);
		for (int a = 0; a < m; a++) z[idx[a]] = col[a];
	}
}

// b <- M^-1 b (直列 : GMRES の中から呼ばれ、順序を固定する)
void pc_apply(const struct precond_t *pc, d_complex_t *b)
{
	const int n1 = pc->n1;

	// z1 = S^-1 (u1 - sum_b B_b D_b^-1 u_b)
	for (int i = 0; i < n1; i++) pc->w1[i] = b[i];
	for (int k = 0; k < pc->nblk; k++) {
		const pcblk_t *pb = &pc->blk[k];
		blk_dsolve(pb, &b[pb->off], pc->wb);
		for (int e = 0; e < pb->nb; e++) {
			pc->w1[pb->bri[e]] = d_sub(pc->w1[pb->bri[e]],
				d_mul(pb->bval[e], pc->wb[pb->bci[e]]));
		}
	}
	if (n1 > 0) lu_solve(n1, pc->slu, pc->spiv, pc->w1);

	// z_b = D_b^-1 (u_b - C_b z1)
	for (int k = 0; k < pc->nblk; k++) {
		const pcblk_t *pb = &pc->blk[k];
		d_complex_t *u = &b[pb->off];
		for (int e = 0; e < pb->nc; e++) {
			u[pb->cri[e]] = d_sub(u[pb->cri[e]],
				d_mul(pb->cval[e], pc->w1[pb->cci[e]]));
		}
		blk_dsolve(pb, u, pc->wb);
		memcpy(u, pc->wb, (size_t)pb->n * sizeof(d_complex_t));
	}
	for (int i = 0; i < n1; i++) b[i] = pc->w1[i];
}

double pc_memory_mb(const struct precond_t *pc)
{
	if (pc == NULL) return 0;
	size_t nz = (size_t)pc->n1 * pc->n1;
	for (int b = 0; b < pc->nblk; b++) {
		nz += pc->blk[b].dof[pc->blk[b].nleaf]
		    + (size_t)pc->blk[b].nb + pc->blk[b].nc;
	}

	return (double)nz * sizeof(d_complex_t) / (1024.0 * 1024.0);
}
