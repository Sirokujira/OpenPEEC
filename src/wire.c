/*
wire.c

導体形状の分割 (幾何段)
- 線導体 (wire / bar) : 区間に分割する。電荷セルは各区間の半分。
- 面導体 (plate)      : 格子点にノードを置き、隣接ノード間を電流セル
                        (幅 = 双対格子の横幅をもつリボン)、各ノードの
                        双対矩形を電荷セルとする (標準的な面 PEEC)。
- node = で束縛された座標を種として、端点/格子点を nodetol 以内で照合し、
  一致しない点には maxid+1 から決定的に自動採番する (ファイル記載順)。
*/

#include "peec.h"

static double dist3(const double *a, const double *b)
{
	const double dx = a[0] - b[0];
	const double dy = a[1] - b[1];
	const double dz = a[2] - b[2];
	return sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

static double norm3(const double *a)
{
	return sqrt((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]));
}

// ノード表から nodetol 以内の点を探す (先に登録されたものが優先)
static int find_gnode(const peec_t *p, const double *x)
{
	for (int i = 0; i < p->ngnode; i++) {
		if (dist3(&p->gxyz[3 * i], x) <= p->nodetol) return i;
	}
	return -1;
}

static void add_gnode(peec_t *p, int id, const double *x)
{
	p->gid[p->ngnode] = id;
	for (int k = 0; k < 3; k++) {
		p->gxyz[3 * p->ngnode + k] = x[k];
	}
	p->ngnode++;
}

// 座標 x のノード id を得る (無ければ自動採番)
static int node_at(peec_t *p, const double *x, int *autoid)
{
	const int found = find_gnode(p, x);
	if (found >= 0) return p->gid[found];

	const int id = ++(*autoid);
	add_gnode(p, id, x);

	return id;
}

// 電荷セルを追加する (軸 x1->x2、幅 wid、横方向 wv)
static void add_chg(peec_t *p, int node, const double *x1, const double *x2,
	double wid, const double *wv, const seg_t *proto)
{
	seg_t *c = &p->chg[p->nchg];
	*c = *proto;
	for (int k = 0; k < 3; k++) {
		c->x1[k] = x1[k];
		c->x2[k] = x2[k];
		c->wv[k] = (wv != NULL) ? wv[k] : 0;
	}
	c->len = dist3(x1, x2);
	c->wid = wid;
	p->chgnode[p->nchg] = node;
	p->nchg++;
}

int wire_build(peec_t *p, FILE *fp_log)
{
	if ((p->nwire <= 0) && (p->nplate <= 0)) return 0;

	// 上限の見積り
	int npoint = p->nnodexyz;
	int nseg = 0;
	int nchg = 0;
	for (int i = 0; i < p->nwire; i++) {
		npoint += p->wire[i].ndiv + 1;
		nseg += p->wire[i].ndiv;
		nchg += 2 * p->wire[i].ndiv;
	}
	for (int i = 0; i < p->nplate; i++) {
		const int na = p->plate[i].ndiva;
		const int nb = p->plate[i].ndivb;
		const int nt = (p->plate[i].ndivt > 1) ? p->plate[i].ndivt : 1;
		npoint += (na + 1) * (nb + 1);
		nseg += nt * ((na * (nb + 1)) + ((na + 1) * nb));
		nchg += (na + 1) * (nb + 1);
	}
	p->gid = (int *)malloc((size_t)npoint * sizeof(int));
	p->gxyz = (double *)malloc((size_t)npoint * 3 * sizeof(double));
	p->seg = (seg_t *)malloc((size_t)nseg * sizeof(seg_t));
	p->chg = (seg_t *)malloc((size_t)nchg * sizeof(seg_t));
	p->chgnode = (int *)malloc((size_t)nchg * sizeof(int));
	if ((p->gid == NULL) || (p->gxyz == NULL) || (p->seg == NULL)
	 || (p->chg == NULL) || (p->chgnode == NULL)) {
		printf("%s\n", "*** memory allocation error (geometry)");
		return 1;
	}
	p->ngnode = 0;
	p->nseg = 0;
	p->nchg = 0;

	// 種 : node = の束縛座標 (重複束縛は警告して先勝ち)
	for (int i = 0; i < p->nnodexyz; i++) {
		const int found = find_gnode(p, &p->ncxyz[3 * i]);
		if (found >= 0) {
			fprintf(fp_log, "*** warning : node %d and node %d are within nodetol, node %d is used\n",
				p->gid[found], p->ncid[i], p->gid[found]);
			continue;
		}
		add_gnode(p, p->ncid[i], &p->ncxyz[3 * i]);
	}

	int autoid = p->maxid;

	// ── 線導体 (wire / bar) ──────────────────────────────────────
	for (int i = 0; i < p->nwire; i++) {
		const wire_t *w = &p->wire[i];
		const int ndiv = w->ndiv;
		int idprev = -1;
		for (int k = 0; k <= ndiv; k++) {
			double x[3];
			for (int c = 0; c < 3; c++) {
				x[c] = w->x1[c] + (w->x2[c] - w->x1[c]) * k / ndiv;
			}
			const int id = node_at(p, x, &autoid);
			if (k > 0) {
				seg_t *s = &p->seg[p->nseg];
				memset(s, 0, sizeof(seg_t));
				s->n1 = idprev;
				s->n2 = id;
				for (int c = 0; c < 3; c++) {
					s->x1[c] = w->x1[c] + (w->x2[c] - w->x1[c]) * (k - 1) / ndiv;
					s->x2[c] = x[c];
				}
				s->len = dist3(s->x1, s->x2);
				s->shape = w->shape;
				s->radius = w->radius;
				s->width = w->width;
				s->thick = w->thick;
				s->sigma = w->sigma;
				s->wid = 0;
				if (w->shape == SHAPE_BAR) {
					// 角線 : 等価半径は幾何平均距離 (GMD) 近似
					//   インダクタンス : GMD = 0.2235 (w + t)  (Grover)
					//     -> Lp = (mu0 l/2pi)[ln(2l/(w+t)) + 0.5] を再現する
					//   容量 : 薄板 (t -> 0) の等価半径 w/4 を (w+t)/4 に拡張
					const double wt = w->width + w->thick;
					s->aL = 0.2235 * wt;
					s->aP = 0.25 * wt;
					s->area = w->width * w->thick;
					s->perim = 2 * wt;
				}
				else {
					s->aL = w->radius;
					s->aP = w->radius;
					s->area = PI * w->radius * w->radius;
					s->perim = 2 * PI * w->radius;
				}
				s->res = (w->sigma > 0) ? s->len / (w->sigma * s->area) : 0;

				// 電荷セル : 区間の半分ずつを両端ノードに与える
				double mid[3];
				for (int c = 0; c < 3; c++) {
					mid[c] = 0.5 * (s->x1[c] + s->x2[c]);
				}
				add_chg(p, s->n1, s->x1, mid, 0, NULL, s);
				add_chg(p, s->n2, mid, s->x2, 0, NULL, s);

				p->nseg++;
			}
			idprev = id;
		}
	}

	// ── 面導体 (plate) ──────────────────────────────────────────
	for (int i = 0; i < p->nplate; i++) {
		const plate_t *pl = &p->plate[i];
		const int na = pl->ndiva;
		const int nb = pl->ndivb;
		const double la = norm3(pl->ea);
		const double lb = norm3(pl->eb);
		const double ha = la / na;
		const double hb = lb / nb;
		double ta[3], tb[3];
		for (int c = 0; c < 3; c++) {
			ta[c] = pl->ea[c] / la;
			tb[c] = pl->eb[c] / lb;
		}
		// 格子点のノード id (行優先 [ia*(nb+1)+ib])
		int *nid = (int *)malloc((size_t)(na + 1) * (nb + 1) * sizeof(int));
		if (nid == NULL) {
			printf("%s\n", "*** memory allocation error (plate)");
			return 1;
		}
		for (int ia = 0; ia <= na; ia++) {
			for (int ib = 0; ib <= nb; ib++) {
				double x[3];
				for (int c = 0; c < 3; c++) {
					x[c] = pl->org[c] + (ia * ha * ta[c]) + (ib * hb * tb[c]);
				}
				nid[(ia * (nb + 1)) + ib] = node_at(p, x, &autoid);
			}
		}

		// 電流セルの雛形 (面導体は厚さ thick、幅はセルごとに決まる)
		seg_t proto;
		memset(&proto, 0, sizeof(seg_t));
		proto.shape = SHAPE_PLATE;
		proto.thick = pl->thick;
		proto.sigma = pl->sigma;

		// 厚み方向の分割 (ndivt >= 2 で体積セル化)。
		// 各層は同じ格子ノードを共有する (層間は各格子点で並列接続)。
		// ノードは元の平面に留まるので node = の束縛はそのまま届く。
		// 層のスタックは元の平面を中心に対称に置く : 板の占有体積が
		// ndivt = 1 のリボン (中立面 = 元の平面) と幾何的に一致する。
		const int nt = (pl->ndivt > 1) ? pl->ndivt : 1;
		const double tl = pl->thick / nt;
		double nv[3];                             // 板の法線 = ta x tb
		nv[0] = (ta[1] * tb[2]) - (ta[2] * tb[1]);
		nv[1] = (ta[2] * tb[0]) - (ta[0] * tb[2]);
		nv[2] = (ta[0] * tb[1]) - (ta[1] * tb[0]);

		// 電流セル : 格子線に沿って隣接ノードを結ぶ。幅は双対格子の横幅
		// (内部の行/列は h、端は h/2)。
		for (int it = 0; it < nt; it++) {
			const double zoff = (it - (0.5 * (nt - 1))) * tl;   // nt = 1 なら 0
		for (int dir = 0; dir < 2; dir++) {
			const int n1 = dir ? nb : na;             // 進行方向の分割数
			const int n2 = dir ? na : nb;             // 横方向の分割数
			const double h1 = dir ? hb : ha;
			const double h2 = dir ? ha : hb;
			const double *t1 = dir ? tb : ta;
			const double *t2 = dir ? ta : tb;
			for (int j = 0; j <= n2; j++) {
				const double wid = ((j == 0) || (j == n2)) ? (0.5 * h2) : h2;
				// 端の行は幅が半分なので、セルの中心も半セル分内側にずらす
				const double off = (j == 0) ? (0.25 * h2)
				                 : (j == n2) ? (-0.25 * h2) : 0;
				for (int k = 0; k < n1; k++) {
					seg_t *s = &p->seg[p->nseg];
					proto.wid = wid;
					*s = proto;
					const int ia1 = dir ? j : k;
					const int ib1 = dir ? k : j;
					const int ia2 = dir ? j : (k + 1);
					const int ib2 = dir ? (k + 1) : j;
					s->n1 = nid[(ia1 * (nb + 1)) + ib1];
					s->n2 = nid[(ia2 * (nb + 1)) + ib2];
					for (int c = 0; c < 3; c++) {
						const double base = pl->org[c] + (j * h2 * t2[c]) + (off * t2[c])
						                  + (zoff * nv[c]);
						s->x1[c] = base + (k * h1 * t1[c]);
						s->x2[c] = base + ((k + 1) * h1 * t1[c]);
						s->wv[c] = t2[c];
					}
					s->len = h1;
					s->vol = (nt > 1);
					s->thick = (nt > 1) ? tl : pl->thick;
					s->area = wid * s->thick;
					s->perim = 2 * (wid + s->thick);
					s->width = wid;
					s->aL = 0.2235 * (wid + s->thick);
					s->aP = 0.25 * (wid + s->thick);
					s->res = (pl->sigma > 0) ? s->len / (pl->sigma * s->area) : 0;
					p->nseg++;
				}
			}
		}
		}

		// 電荷セル : 各格子点の双対矩形 (端は半分、隅は 1/4)
		for (int ia = 0; ia <= na; ia++) {
			for (int ib = 0; ib <= nb; ib++) {
				const double wa = ((ia == 0) || (ia == na)) ? (0.5 * ha) : ha;
				const double wb = ((ib == 0) || (ib == nb)) ? (0.5 * hb) : hb;
				const double oa = (ia == 0) ? (0.25 * ha) : (ia == na) ? (-0.25 * ha) : 0;
				const double ob = (ib == 0) ? (0.25 * hb) : (ib == nb) ? (-0.25 * hb) : 0;
				double c0[3], x1[3], x2[3];
				for (int c = 0; c < 3; c++) {
					c0[c] = pl->org[c] + ((ia * ha) + oa) * ta[c]
					                   + ((ib * hb) + ob) * tb[c];
					x1[c] = c0[c] - (0.5 * wa * ta[c]);
					x2[c] = c0[c] + (0.5 * wa * ta[c]);
				}
				add_chg(p, nid[(ia * (nb + 1)) + ib], x1, x2, wb, tb, &proto);
			}
		}

		free(nid);
	}

	p->maxid = autoid;

	// 近接ノード警告 (マージされなかったが nodetol の 10 倍未満 : 意図しない分離の可能性)
	int nwarn = 0;
	for (int i = 0; (i < p->ngnode) && (nwarn < 10); i++) {
		for (int j = i + 1; (j < p->ngnode) && (nwarn < 10); j++) {
			const double d = dist3(&p->gxyz[3 * i], &p->gxyz[3 * j]);
			if (d < 10 * p->nodetol) {
				fprintf(fp_log, "*** warning : node %d and node %d are separated by only %.3e m\n",
					p->gid[i], p->gid[j], d);
				nwarn++;
			}
		}
	}

	fprintf(fp_log, "geometry : %d wires + %d plates -> %d cells, %d charge cells, %d nodes\n",
		p->nwire, p->nplate, p->nseg, p->nchg, p->ngnode);

	return 0;
}
