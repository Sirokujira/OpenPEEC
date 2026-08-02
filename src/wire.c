/*
wire.c

導体形状の分割 (幾何段)
- 線導体 (wire / bar) : 区間に分割する。電荷セルは各区間の半分。
- 面導体 (plate)      : 格子点にノードを置き、隣接ノード間を電流セル
                        (幅 = 双対格子の横幅をもつリボン)、各ノードの
                        双対矩形を電荷セルとする (標準的な面 PEEC)。
- パネル (quad / disk) : 構造格子 (双一次四辺形 / 極座標) で分割する。
                        plate と同じ双対構成の一般化 : 格子点にノードを
                        置き、格子線に沿う枝の電流セルは隣の格子中線まで
                        の帯 (一般四辺形)、電荷セルは格子点まわりの双対
                        多角形 (辺中点とセル中心を巡る環)。どちらも面を
                        重複なく張り、矩形格子では plate のセルに一致する。
- node = で束縛された座標を種として、端点/格子点/メッシュ頂点を nodetol
  以内で照合し、一致しない点には maxid+1 から決定的に自動採番する
  (ファイル記載順)。
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

// 誘電体の体積セル (枝) を追加する。断面 wid x thk、横方向 wv。
// 枝の直列インピーダンスは過剰容量 1/(jw C_e) (mna.c が diel を見て入れる)。
static void diel_seg(peec_t *p, const double *x1, const double *x2, int n1, int n2,
	const double *wv, double wid, double thk, double epsr)
{
	seg_t *s = &p->seg[p->nseg];
	memset(s, 0, sizeof(seg_t));
	s->shape = SHAPE_PLATE;
	s->vol = 1;
	s->diel = 1;
	s->n1 = n1;
	s->n2 = n2;
	for (int c = 0; c < 3; c++) {
		s->x1[c] = x1[c];
		s->x2[c] = x2[c];
		s->wv[c] = wv[c];
	}
	s->len = dist3(x1, x2);
	s->wid = wid;
	s->thick = thk;
	s->width = wid;
	s->area = wid * thk;
	s->perim = 2 * (wid + thk);
	s->aL = 0.2235 * (wid + thk);
	s->aP = 0.25 * (wid + thk);
	s->sigma = 0;
	s->res = 0;
	s->cexc = EPS0 * (epsr - 1) * s->area / s->len;
	p->nseg++;
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

// ── パネル (quad / disk) : 構造格子のヘルパ ─────────────────────
// 格子点座標。quad は双一次写像、disk は極格子 (i = 0 は中心の 1 点)。
static void panel_node(const panel_t *pa, int i, int j, double *x)
{
	if (pa->kind == PANEL_QUAD) {
		const double u = (double)i / pa->ndiva;
		const double v = (double)j / pa->ndivb;
		for (int c = 0; c < 3; c++) {
			x[c] = ((1 - u) * (1 - v) * pa->v[c])
			     + (u * (1 - v) * pa->v[3 + c])
			     + (u * v * pa->v[6 + c])
			     + ((1 - u) * v * pa->v[9 + c]);
		}
		return;
	}

	// disk : i = リング (0 = 中心)、j = セクタ (mod nsec)
	if (i == 0) {
		for (int c = 0; c < 3; c++) {
			x[c] = pa->org[c];
		}
		return;
	}
	// 面内直交軸 : 法線と最も遠い座標軸から作る (決定的)
	double nn[3], ax[3] = {0, 0, 0}, u[3], v[3];
	const double nl = norm3(pa->nrm);
	for (int c = 0; c < 3; c++) {
		nn[c] = pa->nrm[c] / nl;
	}
	int imin = 0;
	if (fabs(nn[1]) < fabs(nn[imin])) imin = 1;
	if (fabs(nn[2]) < fabs(nn[imin])) imin = 2;
	ax[imin] = 1;
	u[0] = (ax[1] * nn[2]) - (ax[2] * nn[1]);
	u[1] = (ax[2] * nn[0]) - (ax[0] * nn[2]);
	u[2] = (ax[0] * nn[1]) - (ax[1] * nn[0]);
	const double ul = norm3(u);
	v[0] = ((nn[1] * u[2]) - (nn[2] * u[1])) / ul;
	v[1] = ((nn[2] * u[0]) - (nn[0] * u[2])) / ul;
	v[2] = ((nn[0] * u[1]) - (nn[1] * u[0])) / ul;
	const double r = pa->radius * i / pa->ndiva;
	const int jm = ((j % pa->ndivb) + pa->ndivb) % pa->ndivb;
	const double th = 2 * PI * jm / pa->ndivb;
	for (int c = 0; c < 3; c++) {
		x[c] = pa->org[c] + (r * cos(th) * u[c] / ul) + (r * sin(th) * v[c]);
	}
}

// 格子点 (i,j) が存在するか (disk の j は周期的)
static int panel_has(const panel_t *pa, int i, int j)
{
	if (i < 0) return 0;
	if (pa->kind == PANEL_QUAD) {
		return (i <= pa->ndiva) && (j >= 0) && (j <= pa->ndivb);
	}
	return (i <= pa->ndiva);                  // j は mod nsec で常に存在
}

// 格子セル (i+1/2, j+1/2) の中心 (4 隅の平均。disk の i = 0 は中心が重複)
static void panel_cellcen(const panel_t *pa, int i, int j, double *x)
{
	double a[3], b[3], c[3], d[3];
	panel_node(pa, i, j, a);
	panel_node(pa, i + 1, j, b);
	panel_node(pa, i + 1, j + 1, c);
	panel_node(pa, i, j + 1, d);
	for (int k = 0; k < 3; k++) {
		x[k] = 0.25 * (a[k] + b[k] + c[k] + d[k]);
	}
}

static void midp(const double *a, const double *b, double *m)
{
	for (int k = 0; k < 3; k++) {
		m[k] = 0.5 * (a[k] + b[k]);
	}
}

// 電流セル (格子線 (i1,j1)-(i2,j2) の枝) を追加する。
// 帯は隣の格子中線まで : 4 隅は隣接ノードとの中点 (境界はノード自身)。
static void panel_seg(peec_t *p, const panel_t *pa, int n1, int n2,
	int i1, int j1, int i2, int j2)
{
	seg_t *s = &p->seg[p->nseg];
	memset(s, 0, sizeof(seg_t));
	s->shape = SHAPE_POLY;
	s->n1 = n1;
	s->n2 = n2;
	panel_node(pa, i1, j1, s->x1);
	panel_node(pa, i2, j2, s->x2);
	s->len = dist3(s->x1, s->x2);

	// 横方向の隣接格子点 (進行方向と直交する側)
	const int di = j2 - j1;                   // 横方向 = (dj, di) を入れ替えた向き
	const int dj = i2 - i1;
	double q[3];
	int n = 0;
	// マイナス側 : [mid(x1,横-), mid(x2,横-)]、無ければ [x1, x2]
	if (panel_has(pa, i1 - di, j1 - dj) && panel_has(pa, i2 - di, j2 - dj)) {
		panel_node(pa, i1 - di, j1 - dj, q);
		midp(s->x1, q, &s->pv[3 * n]);
		n++;
		panel_node(pa, i2 - di, j2 - dj, q);
		midp(s->x2, q, &s->pv[3 * n]);
		n++;
	}
	else {
		memcpy(&s->pv[3 * n], s->x1, 3 * sizeof(double));
		n++;
		memcpy(&s->pv[3 * n], s->x2, 3 * sizeof(double));
		n++;
	}
	// プラス側 (逆順で環を閉じる)
	if (panel_has(pa, i2 + di, j2 + dj) && panel_has(pa, i1 + di, j1 + dj)) {
		panel_node(pa, i2 + di, j2 + dj, q);
		midp(s->x2, q, &s->pv[3 * n]);
		n++;
		panel_node(pa, i1 + di, j1 + dj, q);
		midp(s->x1, q, &s->pv[3 * n]);
		n++;
	}
	else {
		memcpy(&s->pv[3 * n], s->x2, 3 * sizeof(double));
		n++;
		memcpy(&s->pv[3 * n], s->x1, 3 * sizeof(double));
		n++;
	}
	s->npv = n;
	memcpy(s->papex, &s->pv[0], 3 * sizeof(double));

	const double area = poly_area(s);
	s->wid = area / s->len;                   // 平均幅 (面積/長さ)
	s->thick = pa->thick;
	s->sigma = pa->sigma;
	s->area = s->wid * pa->thick;
	s->perim = 2 * (s->wid + pa->thick);
	s->width = s->wid;
	s->aL = 0.2235 * (s->wid + pa->thick);
	s->aP = 0.25 * (s->wid + pa->thick);
	s->res = (pa->sigma > 0) ? s->len / (pa->sigma * s->area) : 0;
	p->nseg++;
}

// 電荷セル (格子点 (i,j) の双対多角形) を追加する。
// 環は E -> NE -> N -> NW -> W -> SW -> S -> SE (辺中点とセル中心が交互)。
// 境界で欠ける区間には格子点自身を 1 度だけ挿入する (plate の端/隅の
// 半分/1/4 双対矩形と同じ扱い)。
static void panel_chg(peec_t *p, const panel_t *pa, int node, int i, int j)
{
	seg_t *c = &p->chg[p->nchg];
	memset(c, 0, sizeof(seg_t));
	c->shape = SHAPE_POLY;
	double v0[3], q[3];
	panel_node(pa, i, j, v0);

	// 巡回順の 8 要素 : (向き, セルか中点か)
	const int dirs[8][2] = {
		{0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}};
	int n = 0;
	int gap = 0;
	for (int k = 0; k < 8; k++) {
		const int di = dirs[k][0];
		const int dj = dirs[k][1];
		int ok;
		if ((di != 0) && (dj != 0)) {
			// 斜め = 格子セル (i..i+1, j..j+1) の中心
			const int ci = (di > 0) ? i : (i - 1);
			const int cj = (dj > 0) ? j : (j - 1);
			ok = panel_has(pa, ci, cj) && panel_has(pa, ci + 1, cj + 1);
			if (ok) {
				panel_cellcen(pa, ci, cj, q);
			}
		}
		else {
			ok = panel_has(pa, i + di, j + dj);
			if (ok) {
				panel_node(pa, i + di, j + dj, q);
				midp(v0, q, q);
			}
		}
		if (ok) {
			memcpy(&c->pv[3 * n], q, 3 * sizeof(double));
			n++;
		}
		else if (!gap) {
			// 欠けた区間 (境界) : 格子点自身を挿入
			memcpy(&c->pv[3 * n], v0, 3 * sizeof(double));
			n++;
			gap = 1;
		}
	}
	c->npv = n;
	memcpy(c->papex, v0, 3 * sizeof(double));
	memcpy(c->x1, v0, 3 * sizeof(double));
	memcpy(c->x2, v0, 3 * sizeof(double));

	const double area = poly_area(c);
	c->len = sqrt(area);                      // len x wid = 面積
	c->wid = sqrt(area);
	c->thick = pa->thick;
	c->sigma = pa->sigma;
	c->aP = 0.25 * (c->wid + pa->thick);
	p->chgnode[p->nchg] = node;
	p->nchg++;
}

// 双対格子の幅とオフセット (端の行/列は幅が半分で中心が内側へずれる)
static double dualw(int i, int n, double h)
{
	return ((i == 0) || (i == n)) ? (0.5 * h) : h;
}

static double dualo(int i, int n, double h)
{
	return (i == 0) ? (0.25 * h) : (i == n) ? (-0.25 * h) : 0;
}

int wire_build(peec_t *p, FILE *fp_log)
{
	if ((p->nwire <= 0) && (p->nplate <= 0) && (p->npanel <= 0) && (p->ndiel <= 0)) return 0;

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
	for (int i = 0; i < p->npanel; i++) {
		const int na = p->panel[i].ndiva;
		const int nb = p->panel[i].ndivb;
		if (p->panel[i].kind == PANEL_QUAD) {
			npoint += (na + 1) * (nb + 1);
			nseg += (na * (nb + 1)) + ((na + 1) * nb);
			nchg += (na + 1) * (nb + 1);
		}
		else {
			npoint += 1 + (na * nb);          // 中心 + nring x nsec
			nseg += 2 * na * nb;              // 半径方向 + 周方向
			nchg += 1 + (na * nb);
		}
	}
	for (int i = 0; i < p->ndiel; i++) {
		const int na = p->diel[i].ndiva;
		const int nb = p->diel[i].ndivb;
		const int nt = p->diel[i].ndivt;
		npoint += (na + 1) * (nb + 1) * (nt + 1);
		nseg += (na * (nb + 1) * (nt + 1)) + ((na + 1) * nb * (nt + 1))
		      + ((na + 1) * (nb + 1) * nt);
		nchg += (na + 1) * (nb + 1) * (nt + 1);
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

	// ── パネル (quad / disk) : 構造格子の双対セル ─────────────────
	for (int ip = 0; ip < p->npanel; ip++) {
		const panel_t *pa = &p->panel[ip];
		const int na = pa->ndiva;
		const int nb = pa->ndivb;

		if (pa->kind == PANEL_QUAD) {
			int *nid = (int *)malloc((size_t)(na + 1) * (nb + 1) * sizeof(int));
			if (nid == NULL) {
				printf("%s\n", "*** memory allocation error (panel)");
				return 1;
			}
			for (int i = 0; i <= na; i++) {
				for (int j = 0; j <= nb; j++) {
					double x[3];
					panel_node(pa, i, j, x);
					nid[(i * (nb + 1)) + j] = node_at(p, x, &autoid);
				}
			}
			// 電流セル : 格子線に沿う枝 (a 方向 / b 方向)
			for (int j = 0; j <= nb; j++) {
				for (int i = 0; i < na; i++) {
					panel_seg(p, pa, nid[(i * (nb + 1)) + j],
						nid[((i + 1) * (nb + 1)) + j], i, j, i + 1, j);
				}
			}
			for (int i = 0; i <= na; i++) {
				for (int j = 0; j < nb; j++) {
					panel_seg(p, pa, nid[(i * (nb + 1)) + j],
						nid[(i * (nb + 1)) + j + 1], i, j, i, j + 1);
				}
			}
			// 電荷セル : 各格子点の双対多角形
			for (int i = 0; i <= na; i++) {
				for (int j = 0; j <= nb; j++) {
					panel_chg(p, pa, nid[(i * (nb + 1)) + j], i, j);
				}
			}
			free(nid);
		}
		else {
			// disk : ノード index は [0] = 中心、[1 + (i-1)*nsec + j]
			int *nid = (int *)malloc((size_t)(1 + (na * nb)) * sizeof(int));
			if (nid == NULL) {
				printf("%s\n", "*** memory allocation error (panel)");
				return 1;
			}
			{
				double x[3];
				panel_node(pa, 0, 0, x);
				nid[0] = node_at(p, x, &autoid);
			}
			for (int i = 1; i <= na; i++) {
				for (int j = 0; j < nb; j++) {
					double x[3];
					panel_node(pa, i, j, x);
					nid[1 + ((i - 1) * nb) + j] = node_at(p, x, &autoid);
				}
			}
			// 電流セル : 半径方向 (中心からのスポークを含む) と周方向
			for (int j = 0; j < nb; j++) {
				for (int i = 0; i < na; i++) {
					const int m1 = (i == 0) ? 0 : (1 + ((i - 1) * nb) + j);
					const int m2 = 1 + (i * nb) + j;
					panel_seg(p, pa, nid[m1], nid[m2], i, j, i + 1, j);
				}
			}
			for (int i = 1; i <= na; i++) {
				for (int j = 0; j < nb; j++) {
					const int m1 = 1 + ((i - 1) * nb) + j;
					const int m2 = 1 + ((i - 1) * nb) + ((j + 1) % nb);
					panel_seg(p, pa, nid[m1], nid[m2], i, j, i, j + 1);
				}
			}
			// 電荷セル : 中心は特別扱い (スポーク中点とセル中心を巡る環)
			{
				seg_t *c = &p->chg[p->nchg];
				memset(c, 0, sizeof(seg_t));
				c->shape = SHAPE_POLY;
				double v0[3], q[3];
				panel_node(pa, 0, 0, v0);
				int n = 0;
				for (int j = 0; j < nb; j++) {
					panel_node(pa, 1, j, q);
					midp(v0, q, &c->pv[3 * n]);
					n++;
					panel_cellcen(pa, 0, j, &c->pv[3 * n]);
					n++;
				}
				c->npv = n;
				memcpy(c->papex, v0, 3 * sizeof(double));
				memcpy(c->x1, v0, 3 * sizeof(double));
				memcpy(c->x2, v0, 3 * sizeof(double));
				const double area = poly_area(c);
				c->len = sqrt(area);
				c->wid = sqrt(area);
				c->thick = pa->thick;
				c->sigma = pa->sigma;
				c->aP = 0.25 * (c->wid + pa->thick);
				p->chgnode[p->nchg] = nid[0];
				p->nchg++;
			}
			for (int i = 1; i <= na; i++) {
				for (int j = 0; j < nb; j++) {
					panel_chg(p, pa, nid[1 + ((i - 1) * nb) + j], i, j);
				}
			}
			free(nid);
		}
	}

	// ── 誘電体ブリック (dielectric) : 過剰容量の体積セル ─────────────
	// ノードは (na+1) x (nb+1) x (nt+1) の 3 次元格子 (導体と接する面は
	// nodetol マージで導体ノードと共有される)。枝は 3 方向で、それぞれ
	// 双対格子の断面をもつ体積セル。節点の双対矩形 (a-b 面内、ノードの
	// t 面上) を束縛電荷セルとして電位係数に参加させる (内部ノードの
	// 束縛電荷は一様場でほぼ 0 なので、面パネル近似は表面電荷が支配的な
	// 配置で正確)。
	for (int idl = 0; idl < p->ndiel; idl++) {
		const diel_t *dl = &p->diel[idl];
		const int na = dl->ndiva;
		const int nb = dl->ndivb;
		const int nt = dl->ndivt;
		const double la = norm3(dl->ea);
		const double lb = norm3(dl->eb);
		const double ha = la / na;
		const double hb = lb / nb;
		const double ht = dl->thick / nt;
		double ta[3], tb[3], nv[3];
		for (int c = 0; c < 3; c++) {
			ta[c] = dl->ea[c] / la;
			tb[c] = dl->eb[c] / lb;
		}
		nv[0] = (ta[1] * tb[2]) - (ta[2] * tb[1]);
		nv[1] = (ta[2] * tb[0]) - (ta[0] * tb[2]);
		nv[2] = (ta[0] * tb[1]) - (ta[1] * tb[0]);

		int *nid = (int *)malloc((size_t)(na + 1) * (nb + 1) * (nt + 1) * sizeof(int));
		if (nid == NULL) {
			printf("%s\n", "*** memory allocation error (dielectric)");
			return 1;
		}
#define DNID(ia, ib, it) nid[((((ia) * (nb + 1)) + (ib)) * (nt + 1)) + (it)]
		for (int ia = 0; ia <= na; ia++) {
		for (int ib = 0; ib <= nb; ib++) {
			for (int it = 0; it <= nt; it++) {
				double x[3];
				for (int c = 0; c < 3; c++) {
					x[c] = dl->org[c] + (ia * ha * ta[c]) + (ib * hb * tb[c])
					     + (it * ht * nv[c]);
				}
				DNID(ia, ib, it) = node_at(p, x, &autoid);
			}
		}
		}

		// a 方向の枝 (断面 = b, t の双対幅)
		for (int ib = 0; ib <= nb; ib++) {
		for (int it = 0; it <= nt; it++) {
			const double wb = dualw(ib, nb, hb);
			const double ob = dualo(ib, nb, hb);
			const double wt = dualw(it, nt, ht);
			const double ot = dualo(it, nt, ht);
			for (int k = 0; k < na; k++) {
				double x1[3], x2[3];
				for (int c = 0; c < 3; c++) {
					const double base = dl->org[c] + (((ib * hb) + ob) * tb[c])
					                  + (((it * ht) + ot) * nv[c]);
					x1[c] = base + (k * ha * ta[c]);
					x2[c] = base + ((k + 1) * ha * ta[c]);
				}
				diel_seg(p, x1, x2, DNID(k, ib, it), DNID(k + 1, ib, it),
					tb, wb, wt, dl->epsr);
			}
		}
		}
		// b 方向の枝 (断面 = a, t)
		for (int ia = 0; ia <= na; ia++) {
		for (int it = 0; it <= nt; it++) {
			const double wa = dualw(ia, na, ha);
			const double oa = dualo(ia, na, ha);
			const double wt = dualw(it, nt, ht);
			const double ot = dualo(it, nt, ht);
			for (int k = 0; k < nb; k++) {
				double x1[3], x2[3];
				for (int c = 0; c < 3; c++) {
					const double base = dl->org[c] + (((ia * ha) + oa) * ta[c])
					                  + (((it * ht) + ot) * nv[c]);
					x1[c] = base + (k * hb * tb[c]);
					x2[c] = base + ((k + 1) * hb * tb[c]);
				}
				diel_seg(p, x1, x2, DNID(ia, k, it), DNID(ia, k + 1, it),
					ta, wa, wt, dl->epsr);
			}
		}
		}
		// t 方向の枝 (断面 = a, b)
		for (int ia = 0; ia <= na; ia++) {
		for (int ib = 0; ib <= nb; ib++) {
			const double wa = dualw(ia, na, ha);
			const double oa = dualo(ia, na, ha);
			const double wb = dualw(ib, nb, hb);
			const double ob = dualo(ib, nb, hb);
			for (int k = 0; k < nt; k++) {
				double x1[3], x2[3];
				for (int c = 0; c < 3; c++) {
					const double base = dl->org[c] + (((ia * ha) + oa) * ta[c])
					                  + (((ib * hb) + ob) * tb[c]);
					x1[c] = base + (k * ht * nv[c]);
					x2[c] = base + ((k + 1) * ht * nv[c]);
				}
				diel_seg(p, x1, x2, DNID(ia, ib, k), DNID(ia, ib, k + 1),
					ta, wa, wb, dl->epsr);
			}
		}
		}

		// 束縛電荷セル : 各ノードの a-b 面内の双対矩形
		seg_t proto;
		memset(&proto, 0, sizeof(seg_t));
		proto.shape = SHAPE_PLATE;
		proto.thick = ht;
		for (int ia = 0; ia <= na; ia++) {
		for (int ib = 0; ib <= nb; ib++) {
			const double wa = dualw(ia, na, ha);
			const double oa = dualo(ia, na, ha);
			const double wb = dualw(ib, nb, hb);
			const double ob = dualo(ib, nb, hb);
			for (int it = 0; it <= nt; it++) {
				double x1[3], x2[3];
				for (int c = 0; c < 3; c++) {
					const double c0 = dl->org[c] + (((ia * ha) + oa) * ta[c])
					                + (((ib * hb) + ob) * tb[c]) + (it * ht * nv[c]);
					x1[c] = c0 - (0.5 * wa * ta[c]);
					x2[c] = c0 + (0.5 * wa * ta[c]);
				}
				add_chg(p, DNID(ia, ib, it), x1, x2, wb, tb, &proto);
			}
		}
		}
#undef DNID
		free(nid);
	}

	p->maxid = autoid;

	// 地板 (groundplane) : すべてのセルが地板より上 (z >= gpz) にあること。
	// 鏡像法は上半空間でのみ有効なので、下にはみ出す導体はエラーにする。
	if (p->gp) {
		for (int i = 0; i < p->nseg + p->nchg; i++) {
			const seg_t *s = (i < p->nseg) ? &p->seg[i] : &p->chg[i - p->nseg];
			double mz;
			if (s->npv > 0) {
				mz = s->pv[2];
				for (int k = 1; k < s->npv; k++) {
					if (s->pv[(3 * k) + 2] < mz) mz = s->pv[(3 * k) + 2];
				}
			}
			else {
				mz = (s->x1[2] < s->x2[2]) ? s->x1[2] : s->x2[2];
				if (s->wid > 0) {
					// リボンの横方向 (体積セルは厚み方向も) の広がり
					double lo = fabs(0.5 * s->wid * s->wv[2]);
					if (s->vol) {
						double t[3];
						for (int c = 0; c < 3; c++) {
							t[c] = (s->x2[c] - s->x1[c]) / s->len;
						}
						lo += fabs(0.5 * s->thick * ((t[0] * s->wv[1]) - (t[1] * s->wv[0])));
					}
					mz -= lo;
				}
			}
			if (mz < p->gpz - p->nodetol) {
				printf("*** conductor below ground plane (z = %.3e < %.3e)\n", mz, p->gpz);
				fprintf(fp_log, "*** conductor below ground plane (z = %.3e < %.3e)\n", mz, p->gpz);
				return 1;
			}
		}
		fprintf(fp_log, "ground plane : z = %g (image method)\n", p->gpz);
	}

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

	fprintf(fp_log, "geometry : %d wires + %d plates + %d panels + %d dielectric bricks -> %d cells, %d charge cells, %d nodes\n",
		p->nwire, p->nplate, p->npanel, p->ndiel, p->nseg, p->nchg, p->ngnode);

	return 0;
}
