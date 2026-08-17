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

	/*
	圧縮経路 (compression = 1) の容量性 PEEC は電荷 q を陽な未知数にする。

	密経路は節点ブロックに jw C (C = P^-1) を積むが、P の逆行列は密で
	O(ncell^2) のメモリと O(ncell^3) の時間を要し、圧縮の意味が無くなる。
	そこで逆行列を作らず、電荷を未知数に足して

	  KCL 行   : ... + jw q_i = 0
	  電荷行 i : sum_j P(i,j) q_j - v_i = 0

	と組む。電荷行から q = P^-1 v なので KCL には jw (C v)_i が入り、
	密経路と同じ系になる (基準ノードは v = 0 として列が落ちるだけ)。
	P は H 行列で持てるので密な逆行列は現れない。
	*/
	p->offQ = p->nunknown;
	if (p->compress && p->capacitance && (p->nchg > 0)) {
		if (pot_cells(p)) {
			printf("%s\n", "*** memory allocation error (potential cells)");
			return 1;
		}
		p->nunknown = p->offQ + p->ncell;
	}

	fprintf(fp_log, "MNA : %d nodes (reference = node %d), %d unknowns\n",
		p->nnode + 1, p->refnode, p->nunknown);

	return 0;
}

/*
誘電体セルの直列インピーダンス (分極電流の枝)

  Z = 1/Y,  Y = jw eps0 (epsr*(f) - 1) A/len

非分散 (frelax = 0) : epsr* = epsr (1 - j tand) は周波数に依らないので
事前計算済みの cexc / gexc から Y = w (gexc + j cexc)。
分散 (frelax > 0) : 単極 Debye
  epsr*(f) = epsinf + (epss - epsinf)/(1 + j f/frelax)
(因果的で Kramers-Kronig を満たす)。epsr* = a + jb (b <= 0) とすると
  Y = jw eps0 gfac (a - 1 + jb) = w eps0 gfac (-b + j(a - 1))
で、-b >= 0 が緩和損失のコンダクタンスになる。
*/
static d_complex_t zint_diel(const seg_t *s, double omega, double f)
{
	if (s->frelax > 0) {
		const d_complex_t er = d_add(d_complex(s->epsinf, 0),
			d_div(d_complex(s->epss - s->epsinf, 0), d_complex(1, f / s->frelax)));
		return d_inv(d_complex(omega * s->gfac * (-er.i),
		                       omega * s->gfac * (er.r - 1)));
	}

	return d_rmul(1 / omega, d_inv(d_complex(s->gexc, s->cexc)));
}

/*
スタンプの行き先

密行列 (従来) と三つ組リスト (compression = 1 の行列フリー経路) の両方に
同じ組み立てコードから書けるようにする。sp != NULL のときは Lp ブロックを
除いた疎な部分だけを三つ組で集める (Lp は H 行列が持つ)。
*/
typedef struct {
	d_complex_t *a;               // 密行列 (NULL なら疎リスト)
	int          n;
	mna_sparse_t *sp;
} sink_t;

// a[i*n+j] += v (i, j < 0 は基準ノード : スキップ)
static void stamp(sink_t *s, int i, int j, d_complex_t v)
{
	if ((i < 0) || (j < 0)) return;
	if (s->a != NULL) {
		s->a[(size_t)i * s->n + j] = d_add(s->a[(size_t)i * s->n + j], v);
		return;
	}
	mna_sparse_t *sp = s->sp;
	if (sp->nnz >= sp->cap) {
		sp->cap = sp->cap ? (2 * sp->cap) : 1024;
		sp->ri  = (int *)realloc(sp->ri, (size_t)sp->cap * sizeof(int));
		sp->ci  = (int *)realloc(sp->ci, (size_t)sp->cap * sizeof(int));
		sp->val = (d_complex_t *)realloc(sp->val, (size_t)sp->cap * sizeof(d_complex_t));
	}
	sp->ri[sp->nnz] = i;
	sp->ci[sp->nnz] = j;
	sp->val[sp->nnz] = v;
	sp->nnz++;
}

// 共通の組み立て。skiplp = 1 なら密なカーネルブロック (Lp / P) を飛ばす
static void assemble(const peec_t *p, double f, sink_t *s, int skiplp)
{
	const double omega = 2 * PI * f;
	const double kw = p->retardation ? (2 * PI * f / C0) : 0;
	const int *map = p->nodemap;

	// gmin (指定時のみ : 全ノード対地コンダクタンス)
	if (p->gmin > 0) {
		for (int i = 0; i < p->nnode; i++) {
			stamp(s, i, i, d_complex(p->gmin, 0));
		}
	}

	// resistor
	for (int i = 0; i < p->nres; i++) {
		const int i1 = map[p->res[i].n1];
		const int i2 = map[p->res[i].n2];
		const d_complex_t y = d_complex(1 / p->res[i].val, 0);
		stamp(s, i1, i1, y);
		stamp(s, i2, i2, y);
		stamp(s, i1, i2, d_rmul(-1, y));
		stamp(s, i2, i1, d_rmul(-1, y));
	}

	// capacitor
	for (int i = 0; i < p->ncap; i++) {
		const int i1 = map[p->cap[i].n1];
		const int i2 = map[p->cap[i].n2];
		const d_complex_t y = d_complex(0, omega * p->cap[i].val);
		stamp(s, i1, i1, y);
		stamp(s, i2, i2, y);
		stamp(s, i1, i2, d_rmul(-1, y));
		stamp(s, i2, i1, d_rmul(-1, y));
	}

	// inductor (枝電流)
	for (int i = 0; i < p->nind; i++) {
		const int i1 = map[p->ind[i].n1];
		const int i2 = map[p->ind[i].n2];
		const int r = p->offL + i;
		stamp(s, i1, r, d_complex(1, 0));
		stamp(s, r, i1, d_complex(1, 0));
		stamp(s, i2, r, d_complex(-1, 0));
		stamp(s, r, i2, d_complex(-1, 0));
		stamp(s, r, r, d_complex(0, -omega * p->ind[i].val));
	}

	// mutual (M = k sqrt(L1 L2)、ドットは各 inductor の n1 側)
	for (int i = 0; i < p->nmut; i++) {
		const int r1 = p->offL + p->mut[i].l1;
		const int r2 = p->offL + p->mut[i].l2;
		const double m = p->mut[i].k * sqrt(p->ind[p->mut[i].l1].val * p->ind[p->mut[i].l2].val);
		stamp(s, r1, r2, d_complex(0, -omega * m));
		stamp(s, r2, r1, d_complex(0, -omega * m));
	}

	// 電圧源 (枝電流、右辺は mna_rhs_sources で設定。ポート解析時は短絡として働く)
	for (int i = 0; i < p->nsrc; i++) {
		if (!p->src[i].isvsrc) continue;
		const int i1 = map[p->src[i].n1];
		const int i2 = map[p->src[i].n2];
		const int r = p->offV + p->src[i].ibr;
		stamp(s, i1, r, d_complex(1, 0));
		stamp(s, r, i1, d_complex(1, 0));
		stamp(s, i2, r, d_complex(-1, 0));
		stamp(s, r, i2, d_complex(-1, 0));
	}

	// ワイヤ区間 (枝電流 + 密な相互インダクタンスブロック)
	for (int k = 0; k < p->nseg; k++) {
		const int i1 = map[p->seg[k].n1];
		const int i2 = map[p->seg[k].n2];
		const int rk = p->offS + k;
		stamp(s, i1, rk, d_complex(1, 0));
		stamp(s, rk, i1, d_complex(1, 0));
		stamp(s, i2, rk, d_complex(-1, 0));
		stamp(s, rk, i2, d_complex(-1, 0));
		// 内部インピーダンス : 既定は DC 抵抗、skineffect = 1 で表皮効果 + 内部 L。
		// 誘電体セルは過剰分 (epsr* - 1) の直列インピーダンス (zint_diel)
		const d_complex_t zint = p->seg[k].diel
			? zint_diel(&p->seg[k], omega, f)
			: p->skin
			? zint_seg(&p->seg[k], f)
			: d_complex(p->seg[k].res, 0);
		stamp(s, rk, rk, d_rmul(-1, zint));
		if (skiplp) continue;
		for (int m = 0; m < p->nseg; m++) {
			// -j omega Lp (retardation = 1 のとき Lp は複素)
			stamp(s, rk, p->offS + m,
				d_mul(d_complex(0, -omega), p->lp[(size_t)k * p->nseg + m]));
		}
	}

	// 容量性 PEEC
	// (基準ノードは電位 0 なので stamp() 側で行・列とも落ちる)
	if (p->capacitance && (p->ncell > 0)) {
		if (p->offQ < p->nunknown) {
			// 電荷を陽な未知数にする (圧縮経路)。P ブロックは密なので
			// skiplp のときは H 行列に任せて飛ばす。
			for (int i = 0; i < p->ncell; i++) {
				const int mi = map[p->cellid[i]];
				const int ri = p->offQ + i;
				stamp(s, mi, ri, d_complex(0, omega));      // KCL に + jw q
				stamp(s, ri, mi, d_complex(-1, 0));         // 電荷行の - v
				if (skiplp) continue;
				for (int j = 0; j < p->ncell; j++) {
					stamp(s, ri, p->offQ + j, pot_entry(p, i, j, kw));
				}
			}
		}
		else {
			// 密経路 : 節点ブロックに j omega C (C = P^-1)
			for (int i = 0; i < p->ncell; i++) {
				const int mi = map[p->cellid[i]];
				for (int j = 0; j < p->ncell; j++) {
					const int mj = map[p->cellid[j]];
					stamp(s, mi, mj,
						d_mul(d_complex(0, omega), p->cmat[(size_t)i * p->ncell + j]));
				}
			}
		}
	}
}

// 密行列に組み立てる (従来経路)
void mna_assemble(const peec_t *p, double f, d_complex_t *a)
{
	sink_t s;
	s.a = a;
	s.n = p->nunknown;
	s.sp = NULL;
	memset(a, 0, (size_t)p->nunknown * p->nunknown * sizeof(d_complex_t));
	assemble(p, f, &s, 0);
}

/*
Lp ブロックを除いた疎な部分を三つ組で集める (compression = 1)

Lp を除くと非零は要素あたり数個しか無い (接続行列 + 素子スタンプ + 内部
インピーダンスの対角) ので、O(N) で持てる。同じ (i, j) が複数回現れても
mna_apply() が加算するだけなので、重複は問題にならない。
*/
void mna_sparse(const peec_t *p, double f, mna_sparse_t *sp)
{
	sink_t s;
	s.a = NULL;
	s.n = p->nunknown;
	s.sp = sp;
	sp->nnz = 0;
	assemble(p, f, &s, 1);
}

void mna_sparse_free(mna_sparse_t *sp)
{
	free(sp->ri);
	free(sp->ci);
	free(sp->val);
	memset(sp, 0, sizeof(mna_sparse_t));
}

/*
行列フリーの行列ベクトル積 y = A x (compression = 1)

疎な部分は三つ組を足し込み、密なカーネルブロックは H 行列に任せる :
  区間電流 : y[offS + k] += -j omega sum_m Lp[k][m] x[offS + m]
  電荷     : y[offQ + i] +=          sum_j P[i][j]  x[offQ + j]

blk[] が「未知数の塊 (先頭 off、個数 n) + そのカーネルの H 行列 + 係数」を
表す。work は max(blk[].n) 以上の作業領域 (呼び出し側が使い回す)。
三つ組の加算順は登録順で固定、H 行列側もスレッド数不変なので、結果は
スレッド数によらずビット単位で一致する。
*/
void mna_apply(const peec_t *p, const mna_sparse_t *sp,
	const hblock_t *blk, int nblk,
	const d_complex_t *x, d_complex_t *y, d_complex_t *work)
{
	const int n = p->nunknown;

	memset(y, 0, (size_t)n * sizeof(d_complex_t));
	for (int e = 0; e < sp->nnz; e++) {
		y[sp->ri[e]] = d_add(y[sp->ri[e]], d_mul(sp->val[e], x[sp->ci[e]]));
	}
	for (int b = 0; b < nblk; b++) {
		if ((blk[b].h == NULL) || (blk[b].n <= 0)) continue;
		hmat_matvec(blk[b].h, &x[blk[b].off], work);
		for (int k = 0; k < blk[b].n; k++) {
			y[blk[b].off + k] = d_add(y[blk[b].off + k], d_mul(blk[b].scale, work[k]));
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

/*
平面波入射 (planewave) の右辺 — 外部界励振

導体表面の電界積分方程式 E_inc + E_scat = J/sigma を区間に沿って線積分すると

  v_a - v_b - Z_int i_k - jw sum_m Lp[k][m] i_m = -∫_k E_inc・dl

となる (左辺は mna_assemble() が組んだ枝方程式そのもの)。したがって区間 k の
行 (offS + k) に誘起起電力 -∫E_inc・dl を立てればよい。

平面波は到来方向 r^ = (theta, phi) から来るものとする (伝搬方向は -r^)。
exp(jwt) 系では E(r) = E0 e^ exp(+j k r^・r) で、区間内は直線かつ一定方向
なので線積分は厳密に評価できる :

  ∫E・dl = (E0 e^・t^) len exp(j k r^・rc) sinc(k len (r^・t^)/2)   (rc = 中点)

面セル (幅 wid、横方向 wv) では幅方向の平均にも同じ形の sinc が掛かる
(幅で規格化した Lp と同じ扱い)。多角形セルは wv を持たないので軸方向のみ
(セル寸法 << 波長 では差は O((kw)^2))。

地板 (groundplane) があるときは PEC 面での反射波を足す : 到来方向を z 鏡像に
し、偏波の水平成分を反転・垂直成分を保存する。z = gpz で接線成分が打ち消され
(partial.c の鏡像規約と同一)、位相は面上で一致するよう 2 k r^z gpz だけずらす。
*/
static double pw_sinc(double x)
{
	return (fabs(x) < 1e-8) ? 1.0 : (sin(x) / x);
}

// 1 波ぶんの起電力 ∫E・dl
static d_complex_t pw_emf(const seg_t *s, const double *rh, const double *ev,
	d_complex_t e0, double kw)
{
	double tv[3], rc[3];
	for (int c = 0; c < 3; c++) {
		tv[c] = (s->x2[c] - s->x1[c]) / s->len;
		rc[c] = 0.5 * (s->x1[c] + s->x2[c]);
	}
	const double et = (ev[0] * tv[0]) + (ev[1] * tv[1]) + (ev[2] * tv[2]);
	const double rt = (rh[0] * tv[0]) + (rh[1] * tv[1]) + (rh[2] * tv[2]);
	const double rr = (rh[0] * rc[0]) + (rh[1] * rc[1]) + (rh[2] * rc[2]);

	double amp = et * s->len * pw_sinc(0.5 * kw * s->len * rt);

	// 幅方向の平均 (リボン・体積セル : wv が単位ベクトル)
	const double wl = sqrt((s->wv[0] * s->wv[0]) + (s->wv[1] * s->wv[1])
	                     + (s->wv[2] * s->wv[2]));
	if ((s->wid > 0) && (wl > 0)) {
		const double rw = ((rh[0] * s->wv[0]) + (rh[1] * s->wv[1])
		                 + (rh[2] * s->wv[2])) / wl;
		amp *= pw_sinc(0.5 * kw * s->wid * rw);
	}

	return d_rmul(amp, d_mul(e0, d_complex(cos(kw * rr), sin(kw * rr))));
}

void mna_rhs_planewave(const peec_t *p, double f, d_complex_t *b)
{
	const int n = p->nunknown;
	const double kw = 2 * PI * f / C0;

	memset(b, 0, (size_t)n * sizeof(d_complex_t));
	if (!p->pw) return;

	const double th = p->pwth * PI / 180;
	const double ph = p->pwph * PI / 180;
	const double st = sin(th), ct = cos(th), cp = cos(ph), sp = sin(ph);
	const double rh[3] = {st * cp, st * sp, ct};
	const double thh[3] = {ct * cp, ct * sp, -st};
	const double phh[3] = {-sp, cp, 0};
	double ev[3];
	for (int c = 0; c < 3; c++) {
		ev[c] = (p->pwpol == 2) ? phh[c] : thh[c];
	}
	const d_complex_t e0 = d_polar_deg(p->pwamp, p->pwphase);

	// 反射波 (地板) : 到来方向を z 鏡像、偏波の水平成分を反転、
	// 位相は z = gpz で入射波と一致するようずらす
	const double rh2[3] = {rh[0], rh[1], -rh[2]};
	const double ev2[3] = {-ev[0], -ev[1], ev[2]};
	const double psi = 2 * kw * rh[2] * p->gpz;
	const d_complex_t e0r = d_mul(e0, d_complex(cos(psi), sin(psi)));

	for (int k = 0; k < p->nseg; k++) {
		d_complex_t emf = pw_emf(&p->seg[k], rh, ev, e0, kw);
		if (p->gp) {
			emf = d_add(emf, pw_emf(&p->seg[k], rh2, ev2, e0r, kw));
		}
		b[p->offS + k] = d_rmul(-1, emf);
	}
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
