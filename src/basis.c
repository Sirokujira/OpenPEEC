/*
basis.c — 区分線形の電流基底 (basis = 1、細線のみ)

従来の PEEC (basis = 0) は電流が区間ごとに一定 (パルス p_k)、電荷が節点
セルごとに一定 (平均電位法)。basis = 1 は電流に節点セルの傾斜関数を足す :

  J = Σ_k I_k p_k + Σ_n c_n χ_n,  c_n = jω q_n

χ_n は節点 n の電荷セルを作る半区間 h (長さ L_h、節点から外向きの ô_h) 上で
  χ_n = (L_h - s) ô_h / L_n    (s = 節点からの距離、L_n = Σ L_h = carea)
KCL (Σ_out I + jω q = 注入) と合わせると電流は連続な区分線形 (節点と区間中点で
折れ、中点の値が I_k) になり、セル内で一様な電荷密度と連続の式を満たす。
試験関数も {p_k, χ_n} (Galerkin)。χ_n で試験すると
  <χ_n, -∇φ> = φ(節点) - φ̄_n      (φ̄ はセル平均 = Σ P q)
なので、MNA の節点電位を「節点の点電位」と読めば未知数と式の数が揃い、
集中定数・ポートはそのまま点電位に繋がる。電荷を陽な未知数にした MNA
(圧縮経路と同じ並び、mna_numbering) に足すのは

  KVL 行 k  : - jω Σ_n K_kn c_n                  ( = + ω² Σ K_kn q_n )
  電荷行 n  : + jω (Σ_l K_ln I_l + Σ_m X_nm c_m)  ( = jω Kᵀ I - ω² X q )

K_kn = <p_k, A[χ_n]>、X_nm = <χ_n, A[χ_m]> は (mu0/4pi) ∬ e^{-jkR}/R の重み付き
積分 (細線の縮約カーネル、地板の鏡像込み)。ω -> 0 で従来の PEEC に戻る。
試験と展開が同じなので Z 行列は対称 (相反) のまま。平面波の右辺 (電荷行の
<χ_n, E_inc>) と遠方界の放射ベクトル (傾斜関数の寄与) も同じ関数で足すので、
送信 (farfield) と受信 (planewave) の相反も保たれる。

傾斜関数は「節点から長さ t のパルス」の重ね合わせ
  L_h - s = ∫_0^{L_h} [0 <= s < t] dt
なので、重み付き積分は既存の厳密なパルス積分 neumann_pair_k (平行区間は
解析式、遅延込み) を t について Gauss-Legendre 6 点で積分して求める
(t -> 0 で t log t 型。6 点で相対 1e-5 程度 : tools/basis_study.py で確認)。

収束の性質 (tools/basis_study.py の検討) : 粗い格子では従来の基底より正確
(半波長ダイポールの Rin 誤差が片腕 4 分割で 1/3、8 分割で 1/2) だが、電荷が
区分一定のままなので漸近的には従来の互い違い格子 (2 次) より遅い。
*/

#include "peec.h"

// Gauss-Legendre 6 点 ([0, 1] 上)
#define BN 6
static const double gx[BN] = {
	0.0337652428984240, 0.1693953067668677, 0.3806904069584015,
	0.6193095930415985, 0.8306046932331323, 0.9662347571015760};
static const double gw[BN] = {
	0.0856622461895852, 0.1803807865240693, 0.2339569672863455,
	0.2339569672863455, 0.1803807865240693, 0.0856622461895852};

// 電荷サブセル h の節点側の端点と外向きの単位ベクトル
static void sub_frame(const seg_t *h, double node[3], double out[3])
{
	const double *a = h->nend ? h->x2 : h->x1;
	const double *b = h->nend ? h->x1 : h->x2;
	for (int c = 0; c < 3; c++) {
		node[c] = a[c];
		out[c] = (b[c] - a[c]) / h->len;
	}
}

// 節点から外向きに長さ t のパルス
static void sub_piece(const seg_t *h, double t, seg_t *o)
{
	double node[3], out[3];
	sub_frame(h, node, out);
	*o = *h;
	for (int c = 0; c < 3; c++) {
		o->x1[c] = node[c];
		o->x2[c] = node[c] + (t * out[c]);
	}
	o->len = t;
	o->nend = 0;
}

static double dot3b(const double *a, const double *b)
{
	return (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
}

int basis_fill(peec_t *p, double f, FILE *fp_log)
{
	if (!p->basis) return 0;
	if (pot_cells(p)) return 1;

	const int ns = p->nseg;
	const int nc = p->ncell;
	const int nh = p->nchg;
	const double kw = p->retardation ? (2 * PI * f / C0) : 0;
	const double coef = MU0 / (4 * PI);
	const int nimg = gp_nimg(p);

	if (p->bk == NULL) {
		p->bk = (d_complex_t *)malloc((size_t)(ns > 0 ? ns : 1) * (nc > 0 ? nc : 1) * sizeof(d_complex_t));
		p->bx = (d_complex_t *)malloc((size_t)(nc > 0 ? nc : 1) * (nc > 0 ? nc : 1) * sizeof(d_complex_t));
	}
	// 各サブセルの Gauss 点ごとのパルス (長さ L_h x_q)
	seg_t *pc = (seg_t *)malloc((size_t)(nh > 0 ? nh : 1) * BN * sizeof(seg_t));
	if ((p->bk == NULL) || (p->bx == NULL) || (pc == NULL)) {
		printf("%s\n", "*** memory allocation error (basis)");
		free(pc);
		return 1;
	}
	for (int h = 0; h < nh; h++) {
		for (int q = 0; q < BN; q++) {
			sub_piece(&p->chg[h], p->chg[h].len * gx[q], &pc[(h * BN) + q]);
		}
	}

	// K_kn = (mu0/4pi)(1/L_n) Σ_h Σ_img sgn (t_k・M ô_h) L_h Σ_q w_q I(p_k, M piece_h(L_h x_q))
	for (int k = 0; k < ns; k++) {
		const seg_t *s = &p->seg[k];
		double tk[3];
		for (int c = 0; c < 3; c++) {
			tk[c] = (s->x2[c] - s->x1[c]) / s->len;
		}
		for (int n = 0; n < nc; n++) {
			d_complex_t sum = d_complex(0, 0);
			for (int e = p->csoff[n]; e < p->csoff[n + 1]; e++) {
				const int h = p->csidx[e];
				double node[3], out[3];
				sub_frame(&p->chg[h], node, out);
				for (int m = 0; m <= nimg; m++) {
					double om[3] = {out[0], out[1], out[2]};
					if (m > 0) gp_vec(p, m, om);
					const double dt = dot3b(tk, om);
					if (dt == 0) continue;
					d_complex_t acc = d_complex(0, 0);
					for (int q = 0; q < BN; q++) {
						seg_t pm;
						const seg_t *pq = &pc[(h * BN) + q];
						if (m > 0) {
							seg_image(p, pq, m, &pm);
							pq = &pm;
						}
						acc = d_add(acc, d_rmul(gw[q], neumann_pair_k(s, pq, s->aL, pq->aL, kw)));
					}
					const double sg = (m > 0) ? gp_sign(m) : 1.0;
					sum = d_add(sum, d_rmul(sg * dt * p->chg[h].len, acc));
				}
			}
			p->bk[((size_t)k * nc) + n] = d_rmul(coef / p->carea[n], sum);
		}
	}

	// X_nm = (mu0/4pi)/(L_n L_m) Σ_h Σ_g Σ_img sgn (ô_h・M ô_g) L_h L_g
	//        Σ_q Σ_r w_q w_r I(piece_h, M piece_g)   (対称なので上三角だけ計算)
	for (int n = 0; n < nc; n++) {
		for (int mc = n; mc < nc; mc++) {
			d_complex_t sum = d_complex(0, 0);
			for (int e = p->csoff[n]; e < p->csoff[n + 1]; e++) {
				const int h = p->csidx[e];
				double nh3[3], oh[3];
				sub_frame(&p->chg[h], nh3, oh);
				for (int e2 = p->csoff[mc]; e2 < p->csoff[mc + 1]; e2++) {
					const int g = p->csidx[e2];
					double ng3[3], og[3];
					sub_frame(&p->chg[g], ng3, og);
					for (int m = 0; m <= nimg; m++) {
						double om[3] = {og[0], og[1], og[2]};
						if (m > 0) gp_vec(p, m, om);
						const double dt = dot3b(oh, om);
						if (dt == 0) continue;
						d_complex_t acc = d_complex(0, 0);
						for (int q = 0; q < BN; q++) {
							const seg_t *ph = &pc[(h * BN) + q];
							for (int r = 0; r < BN; r++) {
								seg_t pm;
								const seg_t *pg = &pc[(g * BN) + r];
								if (m > 0) {
									seg_image(p, pg, m, &pm);
									pg = &pm;
								}
								acc = d_add(acc, d_rmul(gw[q] * gw[r],
									neumann_pair_k(ph, pg, ph->aL, pg->aL, kw)));
							}
						}
						const double sg = (m > 0) ? gp_sign(m) : 1.0;
						sum = d_add(sum, d_rmul(sg * dt * p->chg[h].len * p->chg[g].len, acc));
					}
				}
			}
			const d_complex_t v = d_rmul(coef / (p->carea[n] * p->carea[mc]), sum);
			p->bx[((size_t)n * nc) + mc] = v;
			p->bx[((size_t)mc * nc) + n] = v;
		}
	}
	free(pc);

	if (!p->blogged) {
		p->blogged = 1;
		fprintf(fp_log, "basis : piecewise-linear current (basis = 1), ramp couplings K = %d x %d, X = %d x %d\n",
			ns, nc, nc, nc);
		fflush(fp_log);
	}
	return 0;
}

/*
∫_0^L (L - s) e^{j β s} ds = L² (1 + jx - e^{jx}) / x²,  x = βL
実部 = 2 sin²(x/2)/x²、虚部 = (x - sin x)/x² (小さい x は級数で桁落ちを避ける)
*/
static d_complex_t ramp_phase(double beta, double len)
{
	const double x = beta * len;
	const double l2 = len * len;
	if (fabs(x) < 1e-3) {
		const double x2 = x * x;
		return d_complex(l2 * (0.5 - (x2 / 24)), l2 * ((x / 6) - (x * x2 / 120)));
	}
	const double sh = sin(0.5 * x);
	const double im = (fabs(x) < 0.1)
		? ((x / 6) - (x * x * x / 120) + (x * x * x * x * x / 5040))
		: ((x - sin(x)) / (x * x));
	return d_complex(l2 * (2 * sh * sh / (x * x)), l2 * im);
}

void basis_rhs_planewave(const peec_t *p, double kw, const pw_wave_t *w, int nw, d_complex_t *b)
{
	if (!p->basis || (p->offQ >= p->nunknown)) return;
	for (int n = 0; n < p->ncell; n++) {
		d_complex_t v = d_complex(0, 0);
		for (int e = p->csoff[n]; e < p->csoff[n + 1]; e++) {
			const seg_t *h = &p->chg[p->csidx[e]];
			double node[3], out[3];
			sub_frame(h, node, out);
			for (int iw = 0; iw < nw; iw++) {
				const double eo = dot3b(w[iw].ev, out);
				if (eo == 0) continue;
				const double rr = kw * dot3b(w[iw].rh, node);
				const d_complex_t ph = d_mul(w[iw].e0, d_complex(cos(rr), sin(rr)));
				const d_complex_t f = ramp_phase(kw * dot3b(w[iw].rh, out), h->len);
				v = d_add(v, d_rmul(eo, d_mul(ph, f)));
			}
		}
		// 電荷行の右辺は + <χ_n, E_inc> (KVL 行の - <p_k, E_inc> と同じ規約)
		b[p->offQ + n] = d_add(b[p->offQ + n], d_rmul(1 / p->carea[n], v));
	}
}

void basis_ff_add(const peec_t *p, const d_complex_t *q, double omega, double kw,
	const double rh[3], d_complex_t nv[3])
{
	if (!p->basis || (q == NULL)) return;
	const int nimg = gp_nimg(p);
	for (int n = 0; n < p->ncell; n++) {
		// 傾斜関数の振幅 c_n = j omega q_n
		const d_complex_t cn = d_mul(d_complex(0, omega), q[n]);
		for (int e = p->csoff[n]; e < p->csoff[n + 1]; e++) {
			const seg_t *h = &p->chg[p->csidx[e]];
			double node[3], out[3];
			sub_frame(h, node, out);
			for (int m = 0; m <= nimg; m++) {
				double nd[3] = {node[0], node[1], node[2]};
				double om[3] = {out[0], out[1], out[2]};
				double sg = 1;
				if (m > 0) {
					gp_pt(p, m, nd);
					gp_vec(p, m, om);
					sg = gp_sign(m);
				}
				const double rr = kw * dot3b(rh, nd);
				const d_complex_t f = ramp_phase(kw * dot3b(rh, om), h->len);
				const d_complex_t wv = d_mul(d_rmul(sg / p->carea[n], cn),
					d_mul(f, d_complex(cos(rr), sin(rr))));
				for (int c = 0; c < 3; c++) {
					nv[c] = d_add(nv[c], d_rmul(om[c], wv));
				}
			}
		}
	}
}
