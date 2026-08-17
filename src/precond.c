/*
precond.c

行列フリー経路 (compression = 1) の前処理 : 葉ブロック消去 + 回路 Schur 補元

MNA の未知数は [回路 (節点・i_L・i_V) | 区間電流] の 2 ブロックで、

  A = | S11  B  |      n1 = offS (回路)、n2 = nseg (区間)
      | C    D  |

D は -(Z_int + jw Lp) で、Lp が密なため D も密。前処理では Lp を
クラスタツリーの**葉の対角ブロックだけ**に落とす (実測で密行列の 1.6〜3.1%)。
すると D はブロック対角になり厳密に消去できて、

  M^-1 は  D_r の LU (葉ごと、m <= HM_LEAF)  と
  Schur 補元 S = S11 - B D^-1 C (n1 x n1 密) の LU

だけで適用できる。適用は前進消去 → Schur 解 → 後退代入の直列 3 段。

設計の根拠 (プロトタイプでの実測、密 LU を前処理に使った GMRES の反復数) :
  - 近傍ブロック全部     : 18〜29 反復 (Lp の 9.6〜32%)
  - 葉の対角ブロックのみ : 24 反復     (Lp の 1.6〜3.1%)   <- 採用
  - Lp の対角のみ        : 30〜36 反復
  - 近傍パターンの ILU(0): 収束せず (棄却)
  - 近傍場の厳密疎 LU    : fill が n^2 の 75% (密と変わらないので棄却)

メモリは葉 LU が ~HM_LEAF x nseg、Schur が n1^2。面格子では n1 ~ nseg/2
なので密 MNA ((offS+nseg)^2) + 密 Lp (nseg^2) に比べ ~1/13 になる。
Schur の n1^2 は O(N^2) のまま残るが、部分要素の充填 O(N^2) を H 行列が
O(kN log N) に落とすのが主目的なので、段階 1 ではこれでよい
(節点数 << 区間数 の体積セル分割では利得はさらに大きい)。

スレッド数不変性 : 葉ブロックの LU は互いに独立なので並列化してよい。
Schur の集約と適用は直列 (GMRES 自体の内積が直列なのでボトルネックではない)。
*/

#include "peec.h"

struct precond_t {
	int n, n1, nseg, offS;
	// 葉ブロック LU (フラット格納)
	int     nleaf;
	int    *lof;                  // [nleaf+1] lidx へのオフセット
	int    *lidx;                 // [nseg] 葉ごとの区間番号 (元の番号)
	size_t *dof;                  // [nleaf+1] dlu へのオフセット (m_r^2)
	d_complex_t *dlu;
	int    *dpiv;                 // [nseg]
	// 疎な結合 (三つ組) : B (回路行 x 区間列)、C (区間行 x 回路列)
	int nb, nc;
	int *bri, *bci; d_complex_t *bval;
	int *cri, *cci; d_complex_t *cval;
	// Schur 補元の LU (n1 x n1)
	d_complex_t *slu;
	int *spiv;
	// 適用の作業領域
	d_complex_t *w2;              // [nseg]
	d_complex_t *w1;              // [n1]
};

void pc_free(struct precond_t *pc)
{
	if (pc == NULL) return;
	free(pc->lof); free(pc->lidx); free(pc->dof); free(pc->dlu); free(pc->dpiv);
	free(pc->bri); free(pc->bci); free(pc->bval);
	free(pc->cri); free(pc->cci); free(pc->cval);
	free(pc->slu); free(pc->spiv);
	free(pc->w2); free(pc->w1);
	free(pc);
}

struct precond_t *pc_build(const peec_t *p, const mna_sparse_t *sp,
	const struct hmat_t *h, double f)
{
	const int n1 = p->offS;
	const int ns = p->nseg;
	const double omega = 2 * PI * f;

	struct precond_t *pc = (struct precond_t *)calloc(1, sizeof(struct precond_t));
	if (pc == NULL) return NULL;
	pc->n = p->nunknown;
	pc->n1 = n1;
	pc->nseg = ns;
	pc->offS = n1;

	// ── 葉ブロックの取り込み ──────────────────────────────────
	pc->nleaf = hmat_nleaf(h);
	pc->lof  = (int *)malloc((size_t)(pc->nleaf + 1) * sizeof(int));
	pc->lidx = (int *)malloc((size_t)(ns > 0 ? ns : 1) * sizeof(int));
	pc->dof  = (size_t *)malloc((size_t)(pc->nleaf + 1) * sizeof(size_t));
	pc->dpiv = (int *)malloc((size_t)(ns > 0 ? ns : 1) * sizeof(int));
	if ((pc->lof == NULL) || (pc->lidx == NULL) || (pc->dof == NULL)
	 || (pc->dpiv == NULL)) {
		pc_free(pc);
		return NULL;
	}
	pc->lof[0] = 0;
	pc->dof[0] = 0;
	for (int r = 0; r < pc->nleaf; r++) {
		const int *idx;
		int m;
		if (hmat_leaf_block(h, r, &idx, &m) == NULL) {
			pc_free(pc);
			return NULL;
		}
		memcpy(&pc->lidx[pc->lof[r]], idx, (size_t)m * sizeof(int));
		pc->lof[r + 1] = pc->lof[r] + m;
		pc->dof[r + 1] = pc->dof[r] + ((size_t)m * m);
	}
	pc->dlu = (d_complex_t *)malloc((pc->dof[pc->nleaf] > 0 ? pc->dof[pc->nleaf] : 1)
		* sizeof(d_complex_t));
	if (pc->dlu == NULL) {
		pc_free(pc);
		return NULL;
	}

	// 区間番号 -> (葉, 葉内位置)
	int *segpos = (int *)malloc((size_t)(ns > 0 ? ns : 1) * sizeof(int));
	int *segleaf = (int *)malloc((size_t)(ns > 0 ? ns : 1) * sizeof(int));
	if ((segpos == NULL) || (segleaf == NULL)) {
		free(segpos); free(segleaf);
		pc_free(pc);
		return NULL;
	}
	for (int r = 0; r < pc->nleaf; r++) {
		for (int e = pc->lof[r]; e < pc->lof[r + 1]; e++) {
			segleaf[pc->lidx[e]] = r;
			segpos[pc->lidx[e]] = e - pc->lof[r];
		}
	}

	// ── 疎部の振り分け : S11 / B / C / D への追加分 ───────────────
	pc->slu = (d_complex_t *)calloc((size_t)(n1 > 0 ? n1 : 1) * (n1 > 0 ? n1 : 1),
		sizeof(d_complex_t));
	pc->spiv = (int *)malloc((size_t)(n1 > 0 ? n1 : 1) * sizeof(int));
	if ((pc->slu == NULL) || (pc->spiv == NULL)) {
		free(segpos); free(segleaf);
		pc_free(pc);
		return NULL;
	}
	int nb = 0, nc = 0;
	for (int e = 0; e < sp->nnz; e++) {
		nb += ((sp->ri[e] < n1) && (sp->ci[e] >= n1));
		nc += ((sp->ri[e] >= n1) && (sp->ci[e] < n1));
	}
	pc->nb = nb;
	pc->nc = nc;
	pc->bri = (int *)malloc((size_t)(nb > 0 ? nb : 1) * sizeof(int));
	pc->bci = (int *)malloc((size_t)(nb > 0 ? nb : 1) * sizeof(int));
	pc->bval = (d_complex_t *)malloc((size_t)(nb > 0 ? nb : 1) * sizeof(d_complex_t));
	pc->cri = (int *)malloc((size_t)(nc > 0 ? nc : 1) * sizeof(int));
	pc->cci = (int *)malloc((size_t)(nc > 0 ? nc : 1) * sizeof(int));
	pc->cval = (d_complex_t *)malloc((size_t)(nc > 0 ? nc : 1) * sizeof(d_complex_t));
	if ((pc->bri == NULL) || (pc->bci == NULL) || (pc->bval == NULL)
	 || (pc->cri == NULL) || (pc->cci == NULL) || (pc->cval == NULL)) {
		free(segpos); free(segleaf);
		pc_free(pc);
		return NULL;
	}

	// 葉ブロック D_r = -jw Lp_rr から初期化し、疎部の (区間, 区間) 要素を足す
	// MSVC の OpenMP 2.0 のためループ変数は事前宣言
	int r;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (r = 0; r < pc->nleaf; r++) {
		const int *idx;
		int m;
		const d_complex_t *lb = hmat_leaf_block(h, r, &idx, &m);
		d_complex_t *d = &pc->dlu[pc->dof[r]];
		const d_complex_t jw = d_complex(0, -omega);
		for (int a = 0; a < m; a++) {
			for (int b = 0; b < m; b++) {
				d[(size_t)a * m + b] = d_mul(jw, lb[(size_t)a * m + b]);
			}
		}
	}

	nb = 0;
	nc = 0;
	for (int e = 0; e < sp->nnz; e++) {
		const int i = sp->ri[e], j = sp->ci[e];
		if ((i < n1) && (j < n1)) {
			pc->slu[(size_t)i * n1 + j] = d_add(pc->slu[(size_t)i * n1 + j], sp->val[e]);
		}
		else if (i < n1) {
			pc->bri[nb] = i;
			pc->bci[nb] = j - n1;
			pc->bval[nb] = sp->val[e];
			nb++;
		}
		else if (j < n1) {
			pc->cri[nc] = i - n1;
			pc->cci[nc] = j;
			pc->cval[nc] = sp->val[e];
			nc++;
		}
		else {
			// (区間, 区間) : Z_int の対角など。同じ葉なら D に足す
			// (葉をまたぐ要素は疎部には存在しない : 対角は必ず同一葉)
			const int k = i - n1, l = j - n1;
			if (segleaf[k] == segleaf[l]) {
				const int lr = segleaf[k];
				const int m = pc->lof[lr + 1] - pc->lof[lr];
				pc->dlu[pc->dof[lr] + ((size_t)segpos[k] * m) + segpos[l]] =
					d_add(pc->dlu[pc->dof[lr] + ((size_t)segpos[k] * m) + segpos[l]],
						sp->val[e]);
			}
		}
	}

	// ── 葉ブロックの LU (互いに独立 = 並列化してよい) ────────────
	int fail = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (r = 0; r < pc->nleaf; r++) {
		const int m = pc->lof[r + 1] - pc->lof[r];
		if (lu_decomp(m, &pc->dlu[pc->dof[r]], &pc->dpiv[pc->lof[r]]) >= 0) {
			fail = 1;
		}
	}
	if (fail) {
		free(segpos); free(segleaf);
		pc_free(pc);
		return NULL;
	}

	// ── Schur 補元 S = S11 - B D^-1 C ─────────────────────────────
	// 葉 r ごとに : 触れる回路未知数の集合 t、C_r (m x t) を D_r で解き、
	// B の該当列から S のクリークに引く。集約は直列 (共有ノードの競合回避)。
	{
		// 区間行ごとの C 要素 (CSR 風)
		int *coff = (int *)calloc((size_t)(ns + 1), sizeof(int));
		int *cidx = (int *)malloc((size_t)(nc > 0 ? nc : 1) * sizeof(int));
		// 区間列ごとの B 要素
		int *boff = (int *)calloc((size_t)(ns + 1), sizeof(int));
		int *bidx = (int *)malloc((size_t)(nb > 0 ? nb : 1) * sizeof(int));
		int *tmap = (int *)malloc((size_t)(n1 > 0 ? n1 : 1) * sizeof(int));
		const int HM = 64;            // 葉の最大セル数 (HM_LEAF) の余裕をみた上限
		int *tset = (int *)malloc((size_t)(4 * HM) * sizeof(int));
		d_complex_t *tblk = (d_complex_t *)malloc((size_t)HM * 4 * HM * sizeof(d_complex_t));
		if ((coff == NULL) || (cidx == NULL) || (boff == NULL) || (bidx == NULL)
		 || (tmap == NULL) || (tset == NULL) || (tblk == NULL)) {
			free(coff); free(cidx); free(boff); free(bidx);
			free(tmap); free(tset); free(tblk);
			free(segpos); free(segleaf);
			pc_free(pc);
			return NULL;
		}
		for (int e = 0; e < nc; e++) coff[pc->cri[e] + 1]++;
		for (int k = 0; k < ns; k++) coff[k + 1] += coff[k];
		for (int e = 0; e < nc; e++) cidx[coff[pc->cri[e]]++] = e;
		for (int k = ns; k > 0; k--) coff[k] = coff[k - 1];
		coff[0] = 0;
		for (int e = 0; e < nb; e++) boff[pc->bci[e] + 1]++;
		for (int k = 0; k < ns; k++) boff[k + 1] += boff[k];
		for (int e = 0; e < nb; e++) bidx[boff[pc->bci[e]]++] = e;
		for (int k = ns; k > 0; k--) boff[k] = boff[k - 1];
		boff[0] = 0;
		for (int j = 0; j < n1; j++) tmap[j] = -1;

		int ok = 1;
		for (r = 0; r < pc->nleaf; r++) {
			const int m = pc->lof[r + 1] - pc->lof[r];
			const int *idx = &pc->lidx[pc->lof[r]];
			// 触れる回路未知数の集合 (順序は走査順で決定的)
			int t = 0;
			for (int a = 0; a < m; a++) {
				for (int q = coff[idx[a]]; q < coff[idx[a] + 1]; q++) {
					const int j = pc->cci[cidx[q]];
					if (tmap[j] < 0) {
						if (t >= 4 * HM) { ok = 0; break; }
						tmap[j] = t;
						tset[t++] = j;
					}
				}
				if (!ok) break;
			}
			if (!ok) break;
			if (t == 0) continue;
			// T = C_r (m x t) を組み、D_r^-1 T を列ごとに解く
			memset(tblk, 0, (size_t)m * t * sizeof(d_complex_t));
			for (int a = 0; a < m; a++) {
				for (int q = coff[idx[a]]; q < coff[idx[a] + 1]; q++) {
					const int e = cidx[q];
					tblk[(size_t)a * t + tmap[pc->cci[e]]] =
						d_add(tblk[(size_t)a * t + tmap[pc->cci[e]]], pc->cval[e]);
				}
			}
			// 列ごとに解く (lu_solve は 1 本ずつなので列を抜き出す)
			for (int c = 0; c < t; c++) {
				d_complex_t col[64];
				for (int a = 0; a < m; a++) col[a] = tblk[(size_t)a * t + c];
				lu_solve(m, &pc->dlu[pc->dof[r]], &pc->dpiv[pc->lof[r]], col);
				for (int a = 0; a < m; a++) tblk[(size_t)a * t + c] = col[a];
			}
			// S[i][tset[c]] -= B[i][idx[a]] * X[a][c]
			for (int a = 0; a < m; a++) {
				for (int q = boff[idx[a]]; q < boff[idx[a] + 1]; q++) {
					const int e = bidx[q];
					const int i = pc->bri[e];
					for (int c = 0; c < t; c++) {
						pc->slu[(size_t)i * n1 + tset[c]] =
							d_sub(pc->slu[(size_t)i * n1 + tset[c]],
								d_mul(pc->bval[e], tblk[(size_t)a * t + c]));
					}
				}
			}
			for (int c = 0; c < t; c++) tmap[tset[c]] = -1;
		}
		free(coff); free(cidx); free(boff); free(bidx);
		free(tmap); free(tset); free(tblk);
		if (!ok) {
			free(segpos); free(segleaf);
			pc_free(pc);
			return NULL;
		}
	}
	free(segpos);
	free(segleaf);

	if ((n1 > 0) && (lu_decomp(n1, pc->slu, pc->spiv) >= 0)) {
		pc_free(pc);
		return NULL;
	}

	pc->w2 = (d_complex_t *)malloc((size_t)(ns > 0 ? ns : 1) * sizeof(d_complex_t));
	pc->w1 = (d_complex_t *)malloc((size_t)(n1 > 0 ? n1 : 1) * sizeof(d_complex_t));
	if ((pc->w2 == NULL) || (pc->w1 == NULL)) {
		pc_free(pc);
		return NULL;
	}

	return pc;
}

// b <- M^-1 b (直列 : GMRES の中から呼ばれ、順序を固定する)
void pc_apply(const struct precond_t *pc, d_complex_t *b)
{
	const int n1 = pc->n1;
	const int ns = pc->nseg;
	d_complex_t *u1 = b;
	d_complex_t *u2 = &b[n1];

	// w2 = D^-1 u2
	for (int r = 0; r < pc->nleaf; r++) {
		const int m = pc->lof[r + 1] - pc->lof[r];
		const int *idx = &pc->lidx[pc->lof[r]];
		d_complex_t col[64];
		for (int a = 0; a < m; a++) col[a] = u2[idx[a]];
		lu_solve(m, &pc->dlu[pc->dof[r]], &pc->dpiv[pc->lof[r]], col);
		for (int a = 0; a < m; a++) pc->w2[idx[a]] = col[a];
	}

	// z1 = S^-1 (u1 - B w2)
	for (int i = 0; i < n1; i++) pc->w1[i] = u1[i];
	for (int e = 0; e < pc->nb; e++) {
		pc->w1[pc->bri[e]] = d_sub(pc->w1[pc->bri[e]],
			d_mul(pc->bval[e], pc->w2[pc->bci[e]]));
	}
	if (n1 > 0) lu_solve(n1, pc->slu, pc->spiv, pc->w1);

	// z2 = D^-1 (u2 - C z1)
	for (int e = 0; e < pc->nc; e++) {
		u2[pc->cri[e]] = d_sub(u2[pc->cri[e]],
			d_mul(pc->cval[e], pc->w1[pc->cci[e]]));
	}
	for (int r = 0; r < pc->nleaf; r++) {
		const int m = pc->lof[r + 1] - pc->lof[r];
		const int *idx = &pc->lidx[pc->lof[r]];
		d_complex_t col[64];
		for (int a = 0; a < m; a++) col[a] = u2[idx[a]];
		lu_solve(m, &pc->dlu[pc->dof[r]], &pc->dpiv[pc->lof[r]], col);
		for (int a = 0; a < m; a++) u2[idx[a]] = col[a];
	}
	for (int i = 0; i < n1; i++) u1[i] = pc->w1[i];
	(void)ns;
}

double pc_memory_mb(const struct precond_t *pc)
{
	if (pc == NULL) return 0;
	size_t nz = pc->dof[pc->nleaf] + ((size_t)pc->n1 * pc->n1)
		+ (size_t)pc->nb + pc->nc;

	return (double)nz * sizeof(d_complex_t) / (1024.0 * 1024.0);
}
