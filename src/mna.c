/*
mna.c

修正節点解析 (MNA)
未知数の並び : [節点電位 (基準除く) | i_L (L 素子) | i_V (電圧源) | i_S (ワイヤ区間)]
ワイヤ区間の枝方程式 : v_a - v_b - R_k i_k - jw Σ_m Lp[k][m] i_m = 0
*/

#include "peec.h"

static void markused(int *used, int id)
{
	used[id] = 1;
}

int mna_numbering(peec_t *p, FILE *fp_log)
{
	const int nid = p->maxid + 1;
	int *used = (int *)malloc((size_t)nid * sizeof(int));
	memset(used, 0, (size_t)nid * sizeof(int));

	for (int i = 0; i < p->nres; i++) {
		markused(used, p->res[i].n1);
		markused(used, p->res[i].n2);
	}
	for (int i = 0; i < p->ncap; i++) {
		markused(used, p->cap[i].n1);
		markused(used, p->cap[i].n2);
	}
	for (int i = 0; i < p->nind; i++) {
		markused(used, p->ind[i].n1);
		markused(used, p->ind[i].n2);
	}
	for (int i = 0; i < p->nsrc; i++) {
		markused(used, p->src[i].n1);
		markused(used, p->src[i].n2);
	}
	for (int i = 0; i < p->nport; i++) {
		markused(used, p->port[i].n1);
		markused(used, p->port[i].n2);
	}
	for (int i = 0; i < p->nseg; i++) {
		markused(used, p->seg[i].n1);
		markused(used, p->seg[i].n2);
	}

	// 基準ノード
	// - 容量性 PEEC を使うときは無限遠が電位の基準になる。実ノードを接地すると
	//   そこから無限遠へ電荷が逃げる経路ができ、構造の総電荷が保存しなくなる。
	//   遅延ありのときこれは放射抵抗の打ち消し (sum q = 0) を壊し、Rin が負に
	//   なるので、node 0 (= 無限遠) 以外は基準にしない。
	//   node 0 が未使用ならどの実ノードも消去しない (容量が電位を確定させる)。
	// - それ以外は従来どおり : ノード 0、無ければ port #1 の n2 (浮遊回路の特異回避)
	const int usecap = (p->capacitance && (p->ncell > 0));
	p->refnode = (usecap || used[0]) ? 0 : p->port[0].n2;

	p->nodemap = (int *)malloc((size_t)nid * sizeof(int));
	int idx = 0;
	for (int id = 0; id < nid; id++) {
		if (id == p->refnode) {
			p->nodemap[id] = -1;
		}
		else if (!used[id]) {
			p->nodemap[id] = -2;
		}
		else {
			p->nodemap[id] = idx++;
		}
	}
	p->nnode = idx;
	free(used);

	// 電圧源の枝番号
	int nvsrc = 0;
	for (int i = 0; i < p->nsrc; i++) {
		if (p->src[i].isvsrc) {
			p->src[i].ibr = nvsrc++;
		}
	}

	p->offL = p->nnode;
	p->offV = p->offL + p->nind;
	p->offS = p->offV + nvsrc;
	p->nunknown = p->offS + p->nseg;

	fprintf(fp_log, "MNA : %d nodes (reference = node %d), %d unknowns\n",
		p->nnode + 1, p->refnode, p->nunknown);

	return 0;
}

// a[i*n+j] += v (i, j < 0 は基準ノード : スキップ)
static void stamp(d_complex_t *a, int n, int i, int j, d_complex_t v)
{
	if ((i >= 0) && (j >= 0)) {
		a[(size_t)i * n + j] = d_add(a[(size_t)i * n + j], v);
	}
}

void mna_assemble(const peec_t *p, double f, d_complex_t *a)
{
	const int n = p->nunknown;
	const double omega = 2 * PI * f;
	const int *map = p->nodemap;

	memset(a, 0, (size_t)n * n * sizeof(d_complex_t));

	// gmin (指定時のみ : 全ノード対地コンダクタンス)
	if (p->gmin > 0) {
		for (int i = 0; i < p->nnode; i++) {
			stamp(a, n, i, i, d_complex(p->gmin, 0));
		}
	}

	// resistor
	for (int i = 0; i < p->nres; i++) {
		const int i1 = map[p->res[i].n1];
		const int i2 = map[p->res[i].n2];
		const d_complex_t y = d_complex(1 / p->res[i].val, 0);
		stamp(a, n, i1, i1, y);
		stamp(a, n, i2, i2, y);
		stamp(a, n, i1, i2, d_rmul(-1, y));
		stamp(a, n, i2, i1, d_rmul(-1, y));
	}

	// capacitor
	for (int i = 0; i < p->ncap; i++) {
		const int i1 = map[p->cap[i].n1];
		const int i2 = map[p->cap[i].n2];
		const d_complex_t y = d_complex(0, omega * p->cap[i].val);
		stamp(a, n, i1, i1, y);
		stamp(a, n, i2, i2, y);
		stamp(a, n, i1, i2, d_rmul(-1, y));
		stamp(a, n, i2, i1, d_rmul(-1, y));
	}

	// inductor (枝電流)
	for (int i = 0; i < p->nind; i++) {
		const int i1 = map[p->ind[i].n1];
		const int i2 = map[p->ind[i].n2];
		const int r = p->offL + i;
		stamp(a, n, i1, r, d_complex(1, 0));
		stamp(a, n, r, i1, d_complex(1, 0));
		stamp(a, n, i2, r, d_complex(-1, 0));
		stamp(a, n, r, i2, d_complex(-1, 0));
		stamp(a, n, r, r, d_complex(0, -omega * p->ind[i].val));
	}

	// mutual (M = k sqrt(L1 L2)、ドットは各 inductor の n1 側)
	for (int i = 0; i < p->nmut; i++) {
		const int r1 = p->offL + p->mut[i].l1;
		const int r2 = p->offL + p->mut[i].l2;
		const double m = p->mut[i].k * sqrt(p->ind[p->mut[i].l1].val * p->ind[p->mut[i].l2].val);
		stamp(a, n, r1, r2, d_complex(0, -omega * m));
		stamp(a, n, r2, r1, d_complex(0, -omega * m));
	}

	// 電圧源 (枝電流、右辺は mna_rhs_sources で設定。ポート解析時は短絡として働く)
	for (int i = 0; i < p->nsrc; i++) {
		if (!p->src[i].isvsrc) continue;
		const int i1 = map[p->src[i].n1];
		const int i2 = map[p->src[i].n2];
		const int r = p->offV + p->src[i].ibr;
		stamp(a, n, i1, r, d_complex(1, 0));
		stamp(a, n, r, i1, d_complex(1, 0));
		stamp(a, n, i2, r, d_complex(-1, 0));
		stamp(a, n, r, i2, d_complex(-1, 0));
	}

	// ワイヤ区間 (枝電流 + 密な相互インダクタンスブロック)
	for (int k = 0; k < p->nseg; k++) {
		const int i1 = map[p->seg[k].n1];
		const int i2 = map[p->seg[k].n2];
		const int rk = p->offS + k;
		stamp(a, n, i1, rk, d_complex(1, 0));
		stamp(a, n, rk, i1, d_complex(1, 0));
		stamp(a, n, i2, rk, d_complex(-1, 0));
		stamp(a, n, rk, i2, d_complex(-1, 0));
		// 内部インピーダンス : 既定は DC 抵抗、skineffect = 1 で表皮効果 + 内部 L
		const d_complex_t zint = p->skin
			? zint_seg(&p->seg[k], f)
			: d_complex(p->seg[k].res, 0);
		stamp(a, n, rk, rk, d_rmul(-1, zint));
		for (int m = 0; m < p->nseg; m++) {
			// -j omega Lp (retardation = 1 のとき Lp は複素)
			stamp(a, n, rk, p->offS + m,
				d_mul(d_complex(0, -omega), p->lp[(size_t)k * p->nseg + m]));
		}
	}

	// 容量性 PEEC : 節点ブロックに j omega C を加える
	// (基準ノードは電位 0 なので stamp() 側で行・列とも落ちる)
	if (p->capacitance && (p->ncell > 0)) {
		for (int i = 0; i < p->ncell; i++) {
			const int mi = map[p->cellid[i]];
			for (int j = 0; j < p->ncell; j++) {
				const int mj = map[p->cellid[j]];
				stamp(a, n, mi, mj,
					d_mul(d_complex(0, omega), p->cmat[(size_t)i * p->ncell + j]));
			}
		}
	}
}

// ポート励振 : n1 に +1A 注入、n2 から -1A (Zin = V(n1) - V(n2))
void mna_rhs_port(const peec_t *p, int iport, d_complex_t *b)
{
	const int n = p->nunknown;
	const int i1 = p->nodemap[p->port[iport].n1];
	const int i2 = p->nodemap[p->port[iport].n2];

	memset(b, 0, (size_t)n * sizeof(d_complex_t));
	if (i1 >= 0) b[i1] = d_add(b[i1], d_complex(1, 0));
	if (i2 >= 0) b[i2] = d_sub(b[i2], d_complex(1, 0));
}

// 独立電源励振 (isource : 電流は電源内部を n1 -> n2 に流れる = SPICE 規約)
void mna_rhs_sources(const peec_t *p, d_complex_t *b)
{
	const int n = p->nunknown;

	memset(b, 0, (size_t)n * sizeof(d_complex_t));
	for (int i = 0; i < p->nsrc; i++) {
		const src_t *s = &p->src[i];
		const d_complex_t v = d_polar_deg(s->amp, s->phase);
		if (s->isvsrc) {
			b[p->offV + s->ibr] = v;
		}
		else {
			const int i1 = p->nodemap[s->n1];
			const int i2 = p->nodemap[s->n2];
			if (i1 >= 0) b[i1] = d_sub(b[i1], v);
			if (i2 >= 0) b[i2] = d_add(b[i2], v);
		}
	}
}
