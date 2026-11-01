/*
spice.c — 準静的な PEEC 等価回路の SPICE ネットリスト出力 (spice = 1 → peec.cir)

部分要素等価回路をそのまま回路素子に書き出す (Ruehli の PEEC の古典的な使い方)。

  導体セル k (n1 -> n2)      : R_k (DC 抵抗) と部分インダクタンス Lp_kk の直列
                               n1 —R— X_k —L— n2 (内部節点 X_k、R = 0 なら L だけ)
  誘電体セル                 : 過剰容量 C_e = eps0 (epsr - 1) A/len と Lp の直列
  相互部分インダクタンス     : K_kl = Lp_kl / sqrt(Lp_kk Lp_ll)  (|K| <= kmin は省く)
  容量性 PEEC (C = P^-1)     : Maxwell の容量行列を枝に分ける
                               節点 - 基準 : Σ_j C_ij、節点 i - j : -C_ij
  集中定数・電源・線路       : そのまま (R / C / L / K / V / I / T)

基準節点 0 は容量性 PEEC なら無限遠 (地板があれば地板)。L の向きは区間の向き
(n1 側に dot) で、Lp の符号 (逆向きの区間で負) は K の符号に入る。部分要素は
準静的 (f = 0) の値で、遅延・表皮効果の周波数依存は表せない (ログとファイルに
注記)。誤差の無い書き出しなので、ネットリストを OpenPEEC の集中定数に読み戻すと
元の PEEC と同じ Zin になる (peec_check.sh が判定する)。
数値は %.16e (17 桁) で書き、double の往復で値が変わらないようにする。
*/

#include "peec.h"

// SPICE の節点名 : 0 はそのまま、それ以外は N<id>
static void nodename(int id, char *buf, size_t n)
{
	if (id == 0) snprintf(buf, n, "0");
	else snprintf(buf, n, "N%d", id);
}

int output_spice(const peec_t *p, const char *fn, FILE *fp_log)
{
	if (!p->spice) return 0;

	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** %s open error\n", fn);
		return 1;
	}

	const int ns = p->nseg;
	const int nc = (p->capacitance && (p->cmat != NULL)) ? p->ncell : 0;
	char a[32], b[32], c[32], d[32];

	fprintf(fp, "* OpenPEEC quasi-static PEEC equivalent circuit\n");
	fprintf(fp, "* title : %s\n", p->title);
	fprintf(fp, "* reference node 0 = %s\n",
		p->gp ? "ground plane (image method)" : (nc > 0 ? "infinity (capacitive PEEC)" : "circuit ground"));
	fprintf(fp, "* partial elements at f = 0 : R = DC resistance, L = partial inductance Lp, K = Lp_kl/sqrt(Lp_kk Lp_ll)\n");
	if (p->retardation) {
		fprintf(fp, "* note : retardation = 1 in the input; the netlist holds the quasi-static (f = 0) elements\n");
	}
	if (p->skin) {
		fprintf(fp, "* note : skineffect = 1 in the input; the netlist holds the DC resistance only\n");
	}
	for (int i = 0; i < p->nport; i++) {
		nodename(p->port[i].n1, a, sizeof(a));
		nodename(p->port[i].n2, b, sizeof(b));
		fprintf(fp, "* port %d : %s %s (z0 = %g ohm)\n", i + 1, a, b, p->port[i].z0);
	}
	fprintf(fp, "* if the DC operating point fails on floating conductors, add a large shunt :\n");
	fprintf(fp, "* .options rshunt=1e12\n");

	// 導体・誘電体セル (R / C_e と Lp の直列)
	fprintf(fp, "*\n* cells : %d (series R or C_e + partial inductance)\n", ns);
	for (int k = 0; k < ns; k++) {
		const seg_t *s = &p->seg[k];
		nodename(s->n1, a, sizeof(a));
		nodename(s->n2, b, sizeof(b));
		const double lkk = p->lp[(size_t)k * ns + k].r;
		if (s->diel) {
			fprintf(fp, "CE%d %s X%d %.16e\n", k + 1, a, k + 1, s->cexc);
			fprintf(fp, "LS%d X%d %s %.16e\n", k + 1, k + 1, b, lkk);
		}
		else if (s->res > 0) {
			fprintf(fp, "RS%d %s X%d %.16e\n", k + 1, a, k + 1, s->res);
			fprintf(fp, "LS%d X%d %s %.16e\n", k + 1, k + 1, b, lkk);
		}
		else {
			fprintf(fp, "LS%d %s %s %.16e\n", k + 1, a, b, lkk);
		}
	}

	// 相互部分インダクタンス
	int nk = 0, ndrop = 0;
	fprintf(fp, "*\n* mutual partial inductances (|k| > %g)\n", p->spkmin);
	for (int k = 0; k < ns; k++) {
		const double lkk = p->lp[(size_t)k * ns + k].r;
		for (int l = k + 1; l < ns; l++) {
			const double lll = p->lp[(size_t)l * ns + l].r;
			const double kv = p->lp[(size_t)k * ns + l].r / sqrt(lkk * lll);
			if (kv == 0) continue;
			if (fabs(kv) <= p->spkmin) {
				ndrop++;
				continue;
			}
			fprintf(fp, "KS%d_%d LS%d LS%d %.16e\n", k + 1, l + 1, k + 1, l + 1, kv);
			nk++;
		}
	}

	// 容量性 PEEC : Maxwell の容量行列を枝に分ける
	int ncap = 0;
	if (nc > 0) {
		fprintf(fp, "*\n* Maxwell capacitance (C = P^-1) : node-to-reference (row sum) and node-to-node (-C_ij)\n");
		for (int i = 0; i < nc; i++) {
			nodename(p->cellid[i], a, sizeof(a));
			double rs = 0;
			for (int j = 0; j < nc; j++) {
				rs += p->cmat[(size_t)i * nc + j].r;
			}
			if (p->cellid[i] != 0) {
				fprintf(fp, "CQ%d_0 %s 0 %.16e\n", i + 1, a, rs);
				ncap++;
			}
			for (int j = i + 1; j < nc; j++) {
				const double v = -p->cmat[(size_t)i * nc + j].r;
				if ((v == 0) || ((p->cellid[i] == 0) && (p->cellid[j] == 0))) continue;
				nodename(p->cellid[j], b, sizeof(b));
				fprintf(fp, "CQ%d_%d %s %s %.16e\n", i + 1, j + 1, a, b, v);
				ncap++;
			}
		}
	}

	// 入力の集中定数・電源・線路
	if ((p->nres + p->ncap + p->nind + p->nmut + p->nsrc + p->ntl) > 0) {
		fprintf(fp, "*\n* lumped elements of the input\n");
	}
	for (int i = 0; i < p->nres; i++) {
		nodename(p->res[i].n1, a, sizeof(a));
		nodename(p->res[i].n2, b, sizeof(b));
		fprintf(fp, "RU%d %s %s %.16e\n", i + 1, a, b, p->res[i].val);
	}
	for (int i = 0; i < p->ncap; i++) {
		nodename(p->cap[i].n1, a, sizeof(a));
		nodename(p->cap[i].n2, b, sizeof(b));
		fprintf(fp, "CU%d %s %s %.16e\n", i + 1, a, b, p->cap[i].val);
	}
	for (int i = 0; i < p->nind; i++) {
		nodename(p->ind[i].n1, a, sizeof(a));
		nodename(p->ind[i].n2, b, sizeof(b));
		fprintf(fp, "LU_%s %s %s %.16e\n", p->ind[i].name, a, b, p->ind[i].val);
	}
	for (int i = 0; i < p->nmut; i++) {
		fprintf(fp, "KU%d LU_%s LU_%s %.16e\n", i + 1,
			p->ind[p->mut[i].l1].name, p->ind[p->mut[i].l2].name, p->mut[i].k);
	}
	for (int i = 0; i < p->nsrc; i++) {
		nodename(p->src[i].n1, a, sizeof(a));
		nodename(p->src[i].n2, b, sizeof(b));
		fprintf(fp, "%cU%d %s %s AC %.16e %.16e\n", p->src[i].isvsrc ? 'V' : 'I', i + 1,
			a, b, p->src[i].amp, p->src[i].phase);
	}
	for (int i = 0; i < p->ntl; i++) {
		const tline_t *t = &p->tl[i];
		nodename(t->n1, a, sizeof(a));
		nodename(t->n2, b, sizeof(b));
		nodename(t->n3, c, sizeof(c));
		nodename(t->n4, d, sizeof(d));
		fprintf(fp, "TU%d %s %s %s %s Z0=%.16e TD=%.16e\n", i + 1, a, b, c, d, t->z0, t->td);
	}
	fprintf(fp, ".end\n");
	fclose(fp);

	fprintf(fp_log, "spice : %s (%d cells, %d couplings, %d dropped with |k| <= %g, %d capacitors)\n",
		fn, ns, nk, ndrop, p->spkmin, ncap);
	if (p->retardation || p->skin) {
		fprintf(fp_log, "*** warning : spice netlist holds quasi-static DC elements (retardation / skin effect not represented)\n");
	}
	return 0;
}
