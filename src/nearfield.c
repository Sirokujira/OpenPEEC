/*
nearfield.c

近傍界の後処理 (nearfield = x1 y1 z1 x2 y2 z2 n1 n2 n3)

port #1 を 1 A で励振したときの区間電流 I_m (p->segi) とセル電荷 q_i
(p->cellq) から、観測点の電界 E と磁界 H を遅延ポテンシャルで求める。
遠方界 (farfield.c) が r -> inf の漸近形だけを見るのに対し、こちらは
1/R, 1/R^2, 1/R^3 の全項を残すので、EMC の近傍界スキャンやプローブ位置での
結合量の見積りに使える。

  A(r)   = (mu0/4pi) sum_m I_m t^_m <g>_m
  phi(r) = (1/4pi eps0) sum_i q_i <g>_i
  E = -j omega A - grad phi,  H = (1/mu0) rot A

exp(+j omega t) 規約なので Green 関数は g(R) = exp(-jkR)/R で、
勾配は grad_r g = -(1 + jkR) exp(-jkR)/R^3 (r - r')。したがって

  E(r) = -j omega mu0/(4pi) sum_m I_m t^_m <g>_m
       +      1/(4pi eps0) sum_i q_i <(1 + jkR) exp(-jkR) (r - r')/R^3>_i
  H(r) =      1/(4pi)      sum_m I_m <(1 + jkR) exp(-jkR) t^_m x (r - r')/R^3>_m

<.>_m はセル上の積分で、電流セルは幅で規格化した線積分 (不変条件 4 と同じ
規約 : リボン / 多角形は (1/wid) x 面積分 が線積分に一致する)、電荷セルは
セル内一様密度の平均 (電位係数 pot_entry と同じ 1/carea 重み) を取る。
求積点はセルの形に依らず polygon.c の cell_qpts() から取る (セル形状と
求積の対応を二重に持たないため)。

k = 0 (retardation = 0) では静的な Coulomb / Biot-Savart にそのまま落ちる。
部分要素も静的なので、場の側だけ遅延を入れない (解と評価を揃える)。
capacitance = 0 では電荷が無いので E は -j omega A の項だけになる。

地板 (groundplane) があるときは鏡像を加える。規約は partial.c / farfield.c と
同一 : 幾何鏡像セル (z -> 2 gpz - z、t^ の z 成分反転) に電流 -I / 電荷 -q。
これで地板上の接線電界が恒等的に 0 になる。地板より下は導体内部なので
観測点を置いても物理的な意味は無い (ログに注意を出す)。

平面波入射 (planewave) があるときは、同じ観測点で平面波による誘起電流・電荷
(p->segipw / p->cellqpw) の散乱界に入射界 (地板があれば反射波も) を足した
全界を nearpw.csv に出す。入射界は右辺 (mna_rhs_planewave) と同じ pw_waves()
から作るので、到来方向・偏波・反射の規約は 1 箇所に閉じている :

  E_inc(r) = sum_w e0_w ev_w exp(+j k rh_w・r)
  H_inc(r) = sum_w -(1/eta0) rh_w x E_w(r)        (伝搬方向は -rh_w)

入射界は物理的な平面波なので retardation に依らず k = omega/c で評価する
(右辺の起電力も同じ k)。散乱界は上の規約どおり retardation = 0 なら静的。

観測点ごとの評価は互いに独立なので OpenMP 並列にしてもリダクションが無く、
スレッド数によらずビット一致する (不変条件の「並列化」の節)。
*/

#include <limits.h>

#include "peec.h"

// セルごとの求積点列 (CSR 風に詰める : セル k の点は [off[k], off[k+1]))
typedef struct {
	int    *off;
	double *px;
	double *wt;
} nfq_t;

static void nfq_free(nfq_t *q)
{
	free(q->off);
	free(q->px);
	free(q->wt);
	q->off = NULL;
	q->px = NULL;
	q->wt = NULL;
}

// n 個のセルの求積点を詰める (点数はセル形状で決まるので 2 パス)
static int nfq_build(const seg_t *cell, int n, nfq_t *q)
{
	double *sx = (double *)malloc((size_t)POLY_NQMAX * 3 * sizeof(double));
	double *sw = (double *)malloc((size_t)POLY_NQMAX * sizeof(double));
	q->off = (int *)calloc((size_t)n + 1, sizeof(int));
	q->px = NULL;
	q->wt = NULL;
	if ((sx == NULL) || (sw == NULL) || (q->off == NULL)) {
		free(sx);
		free(sw);
		return 1;
	}

	for (int k = 0; k < n; k++) {
		q->off[k + 1] = q->off[k] + cell_qpts(&cell[k], 1, 1, sx, sw);
	}
	const int tot = q->off[n];

	q->px = (double *)malloc((size_t)((tot > 0) ? (3 * tot) : 3) * sizeof(double));
	q->wt = (double *)malloc((size_t)((tot > 0) ? tot : 1) * sizeof(double));
	if ((q->px == NULL) || (q->wt == NULL)) {
		free(sx);
		free(sw);
		return 1;
	}

	for (int k = 0; k < n; k++) {
		const int np = cell_qpts(&cell[k], 1, 1, sx, sw);
		memcpy(&q->px[3 * q->off[k]], sx, (size_t)(3 * np) * sizeof(double));
		memcpy(&q->wt[q->off[k]], sw, (size_t)np * sizeof(double));
	}

	free(sx);
	free(sw);
	return 0;
}

/*
観測点 pt におけるセル k の積分核

  gs   = sum_q w_q exp(-jkR)/R
  gv[] = sum_q w_q (1 + jkR) exp(-jkR) (r - r')/R^3     (= -grad gs)

mirror != 0 なら求積点を z -> zmir - z の幾何鏡像に置いて評価する。
*/
static void nf_kernel(const nfq_t *q, int k, const double *pt, double kw,
	int mirror, double zmir, d_complex_t *gs, d_complex_t gv[3])
{
	*gs = d_complex(0, 0);
	for (int c = 0; c < 3; c++) {
		gv[c] = d_complex(0, 0);
	}

	for (int i = q->off[k]; i < q->off[k + 1]; i++) {
		const double sz = q->px[(3 * i) + 2];
		double d[3];
		d[0] = pt[0] - q->px[3 * i];
		d[1] = pt[1] - q->px[(3 * i) + 1];
		d[2] = pt[2] - (mirror ? (zmir - sz) : sz);
		const double r2 = (d[0] * d[0]) + (d[1] * d[1]) + (d[2] * d[2]);
		if (r2 <= 0) continue;                // 求積点に重なった観測点は落とす
		const double r = sqrt(r2);
		const double kr = kw * r;
		const d_complex_t e = d_complex(cos(kr), -sin(kr));    // exp(-jkR)
		const double w = q->wt[i];
		*gs = d_add(*gs, d_rmul(w / r, e));
		const d_complex_t h = d_rmul(w / (r2 * r), d_mul(d_complex(1, kr), e));
		for (int c = 0; c < 3; c++) {
			gv[c] = d_add(gv[c], d_rmul(d[c], h));
		}
	}
}

// 観測点 1 点の E / H (散乱界)。cur[nseg] は区間電流、chg[ncell] はセル電荷
// (NULL なら電荷項なし)。port #1 励振と平面波励振で共通。
static void nf_point(const peec_t *p, const nfq_t *qi, const nfq_t *qc,
	const d_complex_t *cur0, const d_complex_t *chg, const double *pt,
	double f, double kw, d_complex_t *ev, d_complex_t *hv)
{
	const double ca = 2 * PI * f * MU0 / (4 * PI);
	const double cq = 1 / (4 * PI * EPS0);
	const double ch = 1 / (4 * PI);
	const double zmir = 2 * p->gpz;
	const int npass = p->gp ? 2 : 1;

	for (int c = 0; c < 3; c++) {
		ev[c] = d_complex(0, 0);
		hv[c] = d_complex(0, 0);
	}

	for (int pass = 0; pass < npass; pass++) {
		const int mir = (pass == 1);

		// 電流セル : E の -j omega A 項と H
		for (int k = 0; k < p->nseg; k++) {
			const seg_t *s = &p->seg[k];
			d_complex_t cur = cur0[k];
			double tv[3];
			for (int c = 0; c < 3; c++) {
				tv[c] = (s->x2[c] - s->x1[c]) / s->len;
			}
			if (mir) {
				// 鏡像 : 幾何方向の z 反転 + 電流 -I (水平反転・垂直保存)
				tv[2] = -tv[2];
				cur = d_rmul(-1, cur);
			}
			// 幅の規格化 (不変条件 4) : 面セルは (1/wid) x 面積分 = 線積分
			const d_complex_t ic = d_rmul((s->wid > 0) ? (1 / s->wid) : 1, cur);
			d_complex_t gs, gv[3];
			nf_kernel(qi, k, pt, kw, mir, zmir, &gs, gv);

			const d_complex_t w = d_mul(ic, gs);
			for (int c = 0; c < 3; c++) {
				const d_complex_t u = d_rmul(ca * tv[c], w);
				ev[c] = d_add(ev[c], d_complex(u.i, -u.r));    // x (-j)
			}
			// t^ x <(1+jkR) e (r - r')/R^3>
			d_complex_t cx[3];
			cx[0] = d_sub(d_rmul(tv[1], gv[2]), d_rmul(tv[2], gv[1]));
			cx[1] = d_sub(d_rmul(tv[2], gv[0]), d_rmul(tv[0], gv[2]));
			cx[2] = d_sub(d_rmul(tv[0], gv[1]), d_rmul(tv[1], gv[0]));
			for (int c = 0; c < 3; c++) {
				hv[c] = d_add(hv[c], d_rmul(ch, d_mul(ic, cx[c])));
			}
		}

		// 電荷セル : E の -grad phi 項 (鏡像電荷は -q)
		if ((qc == NULL) || (chg == NULL)) continue;
		for (int i = 0; i < p->ncell; i++) {
			if (p->carea[i] <= 0) continue;
			const d_complex_t qcell = chg[i];
			const double scl = (mir ? -cq : cq) / p->carea[i];
			for (int e = p->csoff[i]; e < p->csoff[i + 1]; e++) {
				const int h = p->csidx[e];
				const seg_t *s = &p->chg[h];
				// サブセル電荷 = q_cell x len_h / carea (セル内一様密度)
				const double wn = (s->wid > 0) ? (1 / s->wid) : 1;
				d_complex_t gs, gv[3];
				nf_kernel(qc, h, pt, kw, mir, zmir, &gs, gv);
				const d_complex_t qs = d_rmul(scl * wn, qcell);
				for (int c = 0; c < 3; c++) {
					ev[c] = d_add(ev[c], d_mul(qs, gv[c]));
				}
			}
		}
	}
}

// 平面波の入射界 (地板があれば反射波も) を ev / hv に加える
static void nf_incident(const pw_wave_t *w, int nw, double kinc, const double *pt,
	d_complex_t *ev, d_complex_t *hv)
{
	const double eta0 = MU0 * C0;
	for (int iw = 0; iw < nw; iw++) {
		const double *rh = w[iw].rh;
		const double ph = kinc * ((rh[0] * pt[0]) + (rh[1] * pt[1]) + (rh[2] * pt[2]));
		const d_complex_t a = d_mul(w[iw].e0, d_complex(cos(ph), sin(ph)));
		d_complex_t e[3];
		for (int c = 0; c < 3; c++) {
			e[c] = d_rmul(w[iw].ev[c], a);
			ev[c] = d_add(ev[c], e[c]);
		}
		// H = -(1/eta0) rh x E
		d_complex_t h[3];
		h[0] = d_sub(d_rmul(rh[1], e[2]), d_rmul(rh[2], e[1]));
		h[1] = d_sub(d_rmul(rh[2], e[0]), d_rmul(rh[0], e[2]));
		h[2] = d_sub(d_rmul(rh[0], e[1]), d_rmul(rh[1], e[0]));
		for (int c = 0; c < 3; c++) {
			hv[c] = d_sub(hv[c], d_rmul(1 / eta0, h[c]));
		}
	}
}

static void nf_header(FILE *fp)
{
	fprintf(fp, "frequency[Hz],x[m],y[m],z[m],"
		"Ex_real[V/m],Ex_imag[V/m],Ey_real[V/m],Ey_imag[V/m],Ez_real[V/m],Ez_imag[V/m],"
		"Hx_real[A/m],Hx_imag[A/m],Hy_real[A/m],Hy_imag[A/m],Hz_real[A/m],Hz_imag[A/m],"
		"absE[V/m],absH[A/m]\n");
}

static void nf_write(FILE *fp, double f, const double *pts, const d_complex_t *ev,
	const d_complex_t *hv, int np)
{
	for (int i = 0; i < np; i++) {
		double ae = 0, ah = 0;
		for (int c = 0; c < 3; c++) {
			ae += d_norm(ev[(3 * i) + c]);
			ah += d_norm(hv[(3 * i) + c]);
		}
		fprintf(fp, "%.9e,%.9e,%.9e,%.9e",
			f, pts[3 * i], pts[(3 * i) + 1], pts[(3 * i) + 2]);
		for (int c = 0; c < 3; c++) {
			fprintf(fp, ",%.9e,%.9e", ev[(3 * i) + c].r, ev[(3 * i) + c].i);
		}
		for (int c = 0; c < 3; c++) {
			fprintf(fp, ",%.9e,%.9e", hv[(3 * i) + c].r, hv[(3 * i) + c].i);
		}
		fprintf(fp, ",%.9e,%.9e\n", sqrt(ae), sqrt(ah));
	}
}

int output_near(const peec_t *p, const char *fn, const char *fnpw, FILE *fp_log)
{
	if (p->nnf <= 0) return 0;
	if ((p->segi == NULL) || (p->nseg <= 0) || (p->nport <= 0)) return 0;

	// 観測点を平坦化する (周波数に依らないので一度だけ作る)
	size_t npt = 0;
	for (int ig = 0; ig < p->nnf; ig++) {
		npt += (size_t)(p->nf[ig].n[0] + 1) * (size_t)(p->nf[ig].n[1] + 1)
		     * (size_t)(p->nf[ig].n[2] + 1);
	}
	if ((npt == 0) || (npt > (size_t)INT_MAX)) return 0;

	double *pts = (double *)malloc(npt * 3 * sizeof(double));
	d_complex_t *ev = (d_complex_t *)malloc(npt * 3 * sizeof(d_complex_t));
	d_complex_t *hv = (d_complex_t *)malloc(npt * 3 * sizeof(d_complex_t));
	if ((pts == NULL) || (ev == NULL) || (hv == NULL)) {
		printf("%s\n", "*** memory allocation error (nearfield)");
		free(pts); free(ev); free(hv);
		return 1;
	}

	size_t ip = 0;
	int below = 0;
	for (int ig = 0; ig < p->nnf; ig++) {
		const nf_t *g = &p->nf[ig];
		for (int i0 = 0; i0 <= g->n[0]; i0++) {
		for (int i1 = 0; i1 <= g->n[1]; i1++) {
		for (int i2 = 0; i2 <= g->n[2]; i2++) {
			const int idx[3] = {i0, i1, i2};
			for (int c = 0; c < 3; c++) {
				pts[(3 * ip) + c] = (g->n[c] > 0)
					? (g->p1[c] + ((g->p2[c] - g->p1[c]) * idx[c] / g->n[c]))
					: g->p1[c];
			}
			if (p->gp && (pts[(3 * ip) + 2] < p->gpz)) below = 1;
			ip++;
		}
		}
		}
	}

	// セルの求積点 (幾何なので周波数に依らない)
	nfq_t qi = {NULL, NULL, NULL};
	nfq_t qc = {NULL, NULL, NULL};
	const int hascharge = (p->cellq != NULL) && (p->nchg > 0) && (p->ncell > 0);
	if (nfq_build(p->seg, p->nseg, &qi) ||
	    (hascharge && nfq_build(p->chg, p->nchg, &qc))) {
		printf("%s\n", "*** memory allocation error (nearfield)");
		nfq_free(&qi); nfq_free(&qc);
		free(pts); free(ev); free(hv);
		return 1;
	}
	if (!hascharge) {
		fprintf(fp_log, "*** warning : nearfield without capacitance = 1 (no charge; E keeps only the -j omega A term)\n");
	}
	if (below) {
		fprintf(fp_log, "*** warning : nearfield point below the ground plane (z < %.5e) is inside the conductor\n", p->gpz);
	}

	// 平面波 (planewave) の全界 : 誘起電流があるときだけ
	const int dopw = p->pw && (p->segipw != NULL) && (fnpw != NULL);

	FILE *fp = fopen(fn, "w");
	FILE *fpw = dopw ? fopen(fnpw, "w") : NULL;
	if ((fp == NULL) || (dopw && (fpw == NULL))) {
		printf("*** file %s open error.\n", (fp == NULL) ? fn : fnpw);
		if (fp != NULL) fclose(fp);
		if (fpw != NULL) fclose(fpw);
		nfq_free(&qi); nfq_free(&qc);
		free(pts); free(ev); free(hv);
		return 1;
	}
	nf_header(fp);
	if (fpw != NULL) nf_header(fpw);

	const int np = (int)npt;
	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		const double f = freq_at(p, ifreq);
		// retardation = 0 の解は静的な部分要素から出ているので場も静的に評価する
		const double kw = p->retardation ? (2 * PI * f / C0) : 0;
		const d_complex_t *cur = &p->segi[DIDX(p, ifreq, 0, 0)];
		const d_complex_t *chg = hascharge ? &p->cellq[QIDX(p, ifreq, 0, 0)] : NULL;
		int i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0; i < np; i++) {
			nf_point(p, &qi, hascharge ? &qc : NULL, cur, chg, &pts[3 * i], f, kw,
				&ev[3 * i], &hv[3 * i]);
		}
		nf_write(fp, f, pts, ev, hv, np);

		if (fpw != NULL) {
			// 散乱界 (誘起電流・電荷) + 入射界 (実際の波数 k = omega/c)
			const double kinc = 2 * PI * f / C0;
			pw_wave_t w[2];
			const int nw = pw_waves(p, kinc, w);
			const d_complex_t *curpw = &p->segipw[(size_t)ifreq * p->nseg];
			const d_complex_t *chgpw = (hascharge && (p->cellqpw != NULL))
				? &p->cellqpw[(size_t)ifreq * p->ncell] : NULL;
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (i = 0; i < np; i++) {
				nf_point(p, &qi, (chgpw != NULL) ? &qc : NULL, curpw, chgpw, &pts[3 * i],
					f, kw, &ev[3 * i], &hv[3 * i]);
				nf_incident(w, nw, kinc, &pts[3 * i], &ev[3 * i], &hv[3 * i]);
			}
			nf_write(fpw, f, pts, ev, hv, np);
		}
	}

	fprintf(fp_log, "nearfield : %d grid(s), %d point(s) x %d frequency%s\n",
		p->nnf, np, p->nfreq, dopw ? " (+ plane-wave total field)" : "");
	fflush(fp_log);

	nfq_free(&qi);
	nfq_free(&qc);
	free(pts);
	free(ev);
	free(hv);
	fclose(fp);
	if (fpw != NULL) fclose(fpw);

	return 0;
}
