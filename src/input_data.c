/*
input_data.c

.peec 入力ファイルの読み込み (key = value 形式、OpenFDTD sol/input_data.c の方式)
キー省略時の既定値は後方互換 (従来動作と一致) になるよう初期化する。
*/

#include "peec.h"

#define MAXTOKEN 100

// 動的配列の伸長 (C99 VLA は MSVC 非対応のため malloc/realloc + フラット配列を使う)
#define APPEND(ptr, num, cap, type) \
	do { \
		if ((num) >= (cap)) { \
			(cap) = (cap) ? (2 * (cap)) : 16; \
			(ptr) = (type *)realloc((ptr), (size_t)(cap) * sizeof(type)); \
		} \
	} while (0)

// 非負整数のノード id
static int nodeid(const char *str, int *err)
{
	char *endp = NULL;
	const long v = strtol(str, &endp, 10);
	if ((endp == str) || (*endp != '\0') || (v < 0) || (v > 100000000)) {
		*err = 1;
		return 0;
	}
	return (int)v;
}

static void updatemax(peec_t *p, int id)
{
	if (id > p->maxid) p->maxid = id;
}

/*
らせん (helix) を直線の wire 列に展開する。

  軸    : c を通り単位ベクトル ax の向き
  中心  : c を始点 s の回転面へ射影した c' = c + ((s - c)・ax) ax
  半径  : R = |s - c'|  (始点は軸上にあってはならない)
  基底  : e1 = (s - c')/R、e2 = ax x e1  (e1 -> e2 が +ax に対して右ねじ)
  点 k  : c' + R (cos th_k e1 + sin th_k e2) + (pitch turns k/n) ax、
          th_k = 2 pi turns k/n、k = 0 .. n、n = ceil(分割数 x 巻数)

始点・終点に id が付いていれば node = と同じ束縛座標として登録する
(座標は展開に使ったものそのものなので nodetol マージが必ず成立する)。
*/
static int helix_expand(peec_t *p, helix_t *h, double a, double sigma, int ndiv,
	int *cwire, int *cnode)
{
	double na = sqrt((h->ax[0] * h->ax[0]) + (h->ax[1] * h->ax[1]) + (h->ax[2] * h->ax[2]));
	if (na <= 0) return 1;
	for (int k = 0; k < 3; k++) h->ax[k] /= na;

	double d[3], e1[3], e2[3];
	double dp = 0;
	for (int k = 0; k < 3; k++) {
		d[k] = h->s[k] - h->c[k];
		dp += d[k] * h->ax[k];
	}
	for (int k = 0; k < 3; k++) {
		h->c[k] += dp * h->ax[k];          // 始点の回転面へ射影
		e1[k] = d[k] - (dp * h->ax[k]);
	}
	h->radius = sqrt((e1[0] * e1[0]) + (e1[1] * e1[1]) + (e1[2] * e1[2]));
	if (h->radius <= 0) return 1;
	for (int k = 0; k < 3; k++) e1[k] /= h->radius;
	e2[0] = (h->ax[1] * e1[2]) - (h->ax[2] * e1[1]);
	e2[1] = (h->ax[2] * e1[0]) - (h->ax[0] * e1[2]);
	e2[2] = (h->ax[0] * e1[1]) - (h->ax[1] * e1[0]);

	const double nseg = ceil((ndiv * h->turns) - 1e-9);
	h->nseg = (nseg < 1) ? 1 : (int)nseg;
	h->iwire = p->nwire;
	h->len = 0;

	double prev[3];
	for (int k = 0; k < 3; k++) prev[k] = h->s[k];
	for (int i = 1; i <= h->nseg; i++) {
		const double t = (double)i / h->nseg;
		const double th = 2 * PI * h->turns * t;
		const double z = h->pitch * h->turns * t;
		double pt[3];
		for (int k = 0; k < 3; k++) {
			pt[k] = h->c[k] + (h->radius * ((cos(th) * e1[k]) + (sin(th) * e2[k])))
			      + (z * h->ax[k]);
		}
		APPEND(p->wire, p->nwire, *cwire, wire_t);
		wire_t *w = &p->wire[p->nwire];
		memset(w, 0, sizeof(wire_t));
		for (int k = 0; k < 3; k++) {
			w->x1[k] = prev[k];
			w->x2[k] = pt[k];
		}
		w->shape = SHAPE_ROUND;
		w->radius = a;
		w->sigma = sigma;
		w->ndiv = 1;
		h->len += sqrt(((pt[0] - prev[0]) * (pt[0] - prev[0]))
		             + ((pt[1] - prev[1]) * (pt[1] - prev[1]))
		             + ((pt[2] - prev[2]) * (pt[2] - prev[2])));
		p->nwire++;
		for (int k = 0; k < 3; k++) prev[k] = pt[k];
	}
	for (int k = 0; k < 3; k++) h->e[k] = prev[k];

	// 端点のノード id (node = と同じ束縛座標として登録)
	for (int end = 0; end < 2; end++) {
		const int id = end ? h->nend : h->nstart;
		const double *x = end ? h->e : h->s;
		if (id < 0) continue;
		APPEND(p->ncid, p->nnodexyz, *cnode, int);
		p->ncxyz = (double *)realloc(p->ncxyz, (size_t)(*cnode) * 3 * sizeof(double));
		p->ncid[p->nnodexyz] = id;
		for (int k = 0; k < 3; k++) p->ncxyz[3 * p->nnodexyz + k] = x[k];
		updatemax(p, id);
		p->nnodexyz++;
	}
	return 0;
}

int input_data(FILE *fp, peec_t *p)
{
	int    nline = 0;
	char   strline[BUFSIZ], strkey[BUFSIZ], strsave[BUFSIZ];
	char   *token[MAXTOKEN];
	const char sep[] = " \t";
	const char errfmt2[] = "*** invalid %s data\n";
	int    cres = 0, ccap = 0, cind = 0, cmut = 0, csrc = 0, cport = 0, cwire = 0, cnode = 0;
	int    ctl = 0, cexcit = 0;
	int    cplate = 0, cpanel = 0, cdiel = 0, cnf = 0, chelix = 0;

	// initialize (既定値 : キー省略時は従来動作)
	memset(p, 0, sizeof(peec_t));
	p->nodetol = 1e-8;
	p->gmin = 0;
	p->refnode = -1;
	p->tranatt = 40;              // 過渡応答 : 帯域端でのガウス励振の減衰 [dB]
	p->ctol = 1e-6;               // 圧縮 (compression) の ACA 相対許容誤差

	// read
	while (fgets(strline, sizeof(strline), fp) != NULL) {
		// skip a empty line
		if (strlen(strline) <= 1) continue;

		// skip a comment line
		if (strline[0] == '#') continue;

		// delete "\n"
		if (strstr(strline, "\r\n") != NULL) {
			strline[strlen(strline) - 2] = '\0';
		}
		else if ((strstr(strline, "\r") != NULL) || (strstr(strline, "\n") != NULL)) {
			strline[strlen(strline) - 1] = '\0';
		}

		// "end" -> break
		if (!strncmp(strline, "end", 3)) break;

		// save "strline"
		strcpy(strsave, strline);

		// token ("strline" is destroyed)
		const int ntoken = tokenize(strline, sep, token, MAXTOKEN);

		// check number of data and "=" (exclude header)
		if ((nline > 0) && ((ntoken < 3) || strcmp(token[1], "="))) continue;

		// keyword
		strcpy(strkey, token[0]);

		int err = 0;

		if      (nline == 0) {
			if (strcmp(strkey, "OpenPEEC")) {
				printf("%s\n", "*** not OpenPEEC data");
				return 1;
			}
		}
		else if (!strcmp(strkey, "title")) {
			const char *s = strchr(strsave, '=') + 1;
			while (*s == ' ' || *s == '\t') s++;
			strcpy(p->title, s);
		}
		else if (!strcmp(strkey, "resistor")) {
			if (ntoken < 5) err = 1;
			else {
				APPEND(p->res, p->nres, cres, rc_t);
				rc_t *e = &p->res[p->nres];
				e->n1 = nodeid(token[2], &err);
				e->n2 = nodeid(token[3], &err);
				e->val = atof(token[4]);
				if (e->val <= 0) err = 1;
				updatemax(p, e->n1);
				updatemax(p, e->n2);
				if (!err) p->nres++;
			}
			if (err) {
				printf(errfmt2, "resistor");
				return 1;
			}
		}
		else if (!strcmp(strkey, "capacitor")) {
			if (ntoken < 5) err = 1;
			else {
				APPEND(p->cap, p->ncap, ccap, rc_t);
				rc_t *e = &p->cap[p->ncap];
				e->n1 = nodeid(token[2], &err);
				e->n2 = nodeid(token[3], &err);
				e->val = atof(token[4]);
				// 負の容量も受け付ける (C = P^-1 を枝に分けた PEEC の等価回路では同じ
				// 導体のセル間の枝が負になることがある : spice 出力の読み戻し)。0 は拒否
				if (e->val == 0) err = 1;
				updatemax(p, e->n1);
				updatemax(p, e->n2);
				if (!err) p->ncap++;
			}
			if (err) {
				printf(errfmt2, "capacitor");
				return 1;
			}
		}
		else if (!strcmp(strkey, "inductor")) {
			if (ntoken < 6) err = 1;
			else {
				APPEND(p->ind, p->nind, cind, ind_t);
				ind_t *e = &p->ind[p->nind];
				strncpy(e->name, token[2], NAMELEN - 1);
				e->name[NAMELEN - 1] = '\0';
				e->n1 = nodeid(token[3], &err);
				e->n2 = nodeid(token[4], &err);
				e->val = atof(token[5]);
				if (e->val <= 0) err = 1;
				updatemax(p, e->n1);
				updatemax(p, e->n2);
				if (!err) p->nind++;
			}
			if (err) {
				printf(errfmt2, "inductor");
				return 1;
			}
		}
		else if (!strcmp(strkey, "excitation")) {
			// excitation = ポート番号 (1 始まり) 振幅[A] [位相deg]
			if (ntoken < 4) err = 1;
			else {
				APPEND(p->exc, p->nexc, cexcit, exc_t);
				exc_t *e = &p->exc[p->nexc];
				e->port = nodeid(token[2], &err) - 1;
				e->amp = atof(token[3]);
				e->phase = (ntoken > 4) ? atof(token[4]) : 0;
				if (e->port < 0) err = 1;
				if (!err) p->nexc++;
			}
			if (err) {
				printf(errfmt2, "excitation");
				return 1;
			}
		}
		else if (!strcmp(strkey, "tline")) {
			// tline = n1 n2 n3 n4 Z0 遅延[s] [損失dB]
			// tline = n1 n2 n3 n4 rlgc 長さ[m] R[Ω/m] L[H/m] G[S/m] C[F/m] [Rs[Ω/(m√Hz)] [tanδ]]
			if (ntoken < 8) err = 1;
			else {
				APPEND(p->tl, p->ntl, ctl, tline_t);
				tline_t *e = &p->tl[p->ntl];
				memset(e, 0, sizeof(tline_t));
				e->n1 = nodeid(token[2], &err);
				e->n2 = nodeid(token[3], &err);
				e->n3 = nodeid(token[4], &err);
				e->n4 = nodeid(token[5], &err);
				if (!strcmp(token[6], "rlgc")) {
					if (ntoken < 12) err = 1;
					else {
						e->rlgc = 1;
						e->len = atof(token[7]);
						e->r = atof(token[8]);
						e->l = atof(token[9]);
						e->g = atof(token[10]);
						e->c = atof(token[11]);
						e->rs = (ntoken >= 13) ? atof(token[12]) : 0;
						e->tand = (ntoken >= 14) ? atof(token[13]) : 0;
						if ((e->len <= 0) || (e->l <= 0) || (e->c <= 0) || (e->r < 0) ||
						    (e->g < 0) || (e->rs < 0) || (e->tand < 0)) err = 1;
						// 参考値 (ログ用) : 損失を除いた Z0 と遅延
						if (!err) {
							e->z0 = sqrt(e->l / e->c);
							e->td = e->len * sqrt(e->l * e->c);
						}
						e->att = 1;
					}
				}
				else {
					e->z0 = atof(token[6]);
					e->td = atof(token[7]);
					const double db = (ntoken >= 9) ? atof(token[8]) : 0;
					if ((e->z0 <= 0) || (e->td <= 0) || (db < 0)) err = 1;
					e->att = pow(10, -db / 20);
				}
				if ((e->n1 == e->n2) || (e->n3 == e->n4)) err = 1;
				updatemax(p, e->n1);
				updatemax(p, e->n2);
				updatemax(p, e->n3);
				updatemax(p, e->n4);
				if (!err) p->ntl++;
			}
			if (err) {
				printf(errfmt2, "tline");
				return 1;
			}
		}
		else if (!strcmp(strkey, "mutual")) {
			if (ntoken < 5) err = 1;
			else {
				APPEND(p->mut, p->nmut, cmut, mut_t);
				mut_t *e = &p->mut[p->nmut];
				strncpy(e->name1, token[2], NAMELEN - 1);
				e->name1[NAMELEN - 1] = '\0';
				strncpy(e->name2, token[3], NAMELEN - 1);
				e->name2[NAMELEN - 1] = '\0';
				e->k = atof(token[4]);
				e->l1 = e->l2 = -1;
				// 0 < |k| <= 1 (負の k はドットの向きの反転 = 逆向きの結合。
				// PEEC の部分相互インダクタンスは区間の向きで負になる : spice 出力の読み戻し)
				if ((e->k == 0) || (fabs(e->k) > 1)) err = 1;
				if (!err) p->nmut++;
			}
			if (err) {
				printf(errfmt2, "mutual");
				return 1;
			}
		}
		else if (!strcmp(strkey, "vsource") || !strcmp(strkey, "isource")) {
			if (ntoken < 5) err = 1;
			else {
				APPEND(p->src, p->nsrc, csrc, src_t);
				src_t *e = &p->src[p->nsrc];
				e->isvsrc = !strcmp(strkey, "vsource");
				e->n1 = nodeid(token[2], &err);
				e->n2 = nodeid(token[3], &err);
				e->amp = atof(token[4]);
				e->phase = (ntoken > 5) ? atof(token[5]) : 0;
				e->ibr = -1;
				updatemax(p, e->n1);
				updatemax(p, e->n2);
				if (!err) p->nsrc++;
			}
			if (err) {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "node")) {
			if (ntoken < 6) err = 1;
			else {
				APPEND(p->ncid, p->nnodexyz, cnode, int);
				// 座標配列は id 配列と同じ個数だけ伸ばす
				p->ncxyz = (double *)realloc(p->ncxyz, (size_t)cnode * 3 * sizeof(double));
				p->ncid[p->nnodexyz] = nodeid(token[2], &err);
				for (int k = 0; k < 3; k++) {
					p->ncxyz[3 * p->nnodexyz + k] = atof(token[3 + k]);
				}
				updatemax(p, p->ncid[p->nnodexyz]);
				if (!err) p->nnodexyz++;
			}
			if (err) {
				printf(errfmt2, "node");
				return 1;
			}
		}
		else if (!strcmp(strkey, "wire") || !strcmp(strkey, "bar")) {
			// wire = x1 y1 z1 x2 y2 z2 半径     導電率 分割数 [断面分割数]
			// bar  = x1 y1 z1 x2 y2 z2 幅 厚さ  導電率 分割数 [分割数w [分割数t]]
			const int isbar = !strcmp(strkey, "bar");
			if (ntoken < (isbar ? 12 : 11)) err = 1;
			else {
				APPEND(p->wire, p->nwire, cwire, wire_t);
				wire_t *e = &p->wire[p->nwire];
				memset(e, 0, sizeof(wire_t));
				for (int k = 0; k < 3; k++) {
					e->x1[k] = atof(token[2 + k]);
					e->x2[k] = atof(token[5 + k]);
				}
				if (isbar) {
					e->shape = SHAPE_BAR;
					e->width = atof(token[8]);
					e->thick = atof(token[9]);
					e->sigma = atof(token[10]);
					e->ndiv = atoi(token[11]);
					e->ndivw = (ntoken >= 13) ? atoi(token[12]) : 1;
					e->ndivt = (ntoken >= 14) ? atoi(token[13]) : 1;
					if ((e->width <= 0) || (e->thick <= 0) || (e->ndivw < 1) || (e->ndivt < 1)) err = 1;
				}
				else {
					e->shape = SHAPE_ROUND;
					e->radius = atof(token[8]);
					e->sigma = atof(token[9]);
					e->ndiv = atoi(token[10]);
					// 断面分割 (n x n 格子のうち円内のセル。1 = 分割しない)
					e->ndivw = (ntoken >= 12) ? atoi(token[11]) : 1;
					if ((e->radius <= 0) || (e->ndivw < 1)) err = 1;
				}
				const double dx = e->x2[0] - e->x1[0];
				const double dy = e->x2[1] - e->x1[1];
				const double dz = e->x2[2] - e->x1[2];
				const double len = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if ((e->sigma < 0) || (e->ndiv < 1) || (len <= 0)) err = 1;
				if (!err) p->nwire++;
			}
			if (err) {
				printf(errfmt2, strkey);
				return 1;
			}
		}
		else if (!strcmp(strkey, "helix")) {
			// helix = cx cy cz ax ay az sx sy sz ピッチ 巻数 半径 導電率 分割数 [n1 n2]
			//   c : 軸上の点、ax : 軸方向、s : 始点 (導体上、コイル半径は s と軸の距離)
			//   ピッチ : 1 周あたりの軸方向の進み (負 = 左巻き、0 = 平面ループ)
			//   分割数 : 1 周あたりの区間数 (巻数が小数なら ceil(分割数 x 巻数) 区間)
			//   n1 n2 : 始点・終点に付けるノード id (省略可)
			if (ntoken < 16) err = 1;
			else {
				APPEND(p->helix, p->nhelix, chelix, helix_t);
				helix_t *h = &p->helix[p->nhelix];
				memset(h, 0, sizeof(helix_t));
				for (int k = 0; k < 3; k++) {
					h->c[k] = atof(token[2 + k]);
					h->ax[k] = atof(token[5 + k]);
					h->s[k] = atof(token[8 + k]);
				}
				h->pitch = atof(token[11]);
				h->turns = atof(token[12]);
				const double a = atof(token[13]);
				const double sigma = atof(token[14]);
				const int ndiv = atoi(token[15]);
				h->nstart = (ntoken >= 17) ? nodeid(token[16], &err) : -1;
				h->nend = (ntoken >= 18) ? nodeid(token[17], &err) : -1;
				if ((h->turns <= 0) || (a <= 0) || (sigma < 0) || (ndiv < 1)) err = 1;
				if (!err) err = helix_expand(p, h, a, sigma, ndiv, &cwire, &cnode);
				if (!err) p->nhelix++;
			}
			if (err) {
				printf(errfmt2, "helix");
				return 1;
			}
		}
		else if (!strcmp(strkey, "plate")) {
			// plate = ox oy oz  ax ay az  bx by bz  厚さ 導電率 ndiv_a ndiv_b [ndiv_t]
			// ndiv_t (省略時 1) >= 2 で厚み方向にも分割する (体積セル)
			if (ntoken < 15) err = 1;
			else {
				APPEND(p->plate, p->nplate, cplate, plate_t);
				plate_t *e = &p->plate[p->nplate];
				memset(e, 0, sizeof(plate_t));
				for (int k = 0; k < 3; k++) {
					e->org[k] = atof(token[2 + k]);
					e->ea[k] = atof(token[5 + k]);
					e->eb[k] = atof(token[8 + k]);
				}
				e->thick = atof(token[11]);
				e->sigma = atof(token[12]);
				e->ndiva = atoi(token[13]);
				e->ndivb = atoi(token[14]);
				e->ndivt = (ntoken >= 16) ? atoi(token[15]) : 1;
				const double la = sqrt((e->ea[0] * e->ea[0]) + (e->ea[1] * e->ea[1]) + (e->ea[2] * e->ea[2]));
				const double lb = sqrt((e->eb[0] * e->eb[0]) + (e->eb[1] * e->eb[1]) + (e->eb[2] * e->eb[2]));
				const double ab = (e->ea[0] * e->eb[0]) + (e->ea[1] * e->eb[1]) + (e->ea[2] * e->eb[2]);
				if ((la <= 0) || (lb <= 0) || (e->thick <= 0) || (e->sigma < 0)
				 || (e->ndiva < 1) || (e->ndivb < 1) || (e->ndivt < 1)) err = 1;
				// 2 辺は直交していること (矩形セルを前提にしている)
				if (!err && (fabs(ab) > 1e-9 * la * lb)) {
					printf("%s\n", "*** plate : the two edge vectors must be perpendicular");
					return 1;
				}
				if (!err) p->nplate++;
			}
			if (err) {
				printf(errfmt2, "plate");
				return 1;
			}
		}
		else if (!strcmp(strkey, "dielectric")) {
			// dielectric = ox oy oz  ax ay az  bx by bz  厚さ epsr
			//              ndiv_a ndiv_b ndiv_t [tand [eps_inf 分散]]
			// org を底面として法線 (ea x eb) 方向へ thick 押し出した直方体。
			// tand (省略時 0) は誘電正接 : epsr* = epsr (1 - j tand)。
			//
			// eps_inf の後ろに分散を書くと多極分散になる (どの極も因果的で
			// Kramers-Kronig を満たす)。2 通りの書き方を残り語数で判別する :
			//
			//   f_relax                     (1 語)  単極 Debye : 強度は epsr - eps_inf
			//   極の並び                    (それ以外) 多極 : 強度を明示
			//
			// 極の並びは次の 2 種類を任意の順・組合せで書ける :
			//
			//   f_k de_k                    Debye 極   de_k/(1 + j f/f_k)
			//   lorentz f_k de_k g_k        Lorentz 極 de_k f_k^2/(f_k^2 - f^2 + j f g_k)
			//
			//   epsr*(f) = eps_inf + sum_k (極の寄与)
			//
			// Lorentz 極は共振周波数 f_k、減衰幅 g_k [Hz] (> 0。Q = f_k/g_k で、
			// 高 Q では損失ピークの半値全幅にほぼ等しい) の共振型で、
			// f_k の少し上では epsr' が 1 を割り負にもなる (Debye では
			// 起きない)。どちらの極も f -> 0 で de_k に戻るので、直流極限
			// eps_inf + sum de_k が epsr と一致している必要がある (書き間違いを
			// 静かに通さないため一致を要求する)。損失は極が生むので定数 tand
			// との併用は二重計上 : tand = 0 を書くこと。
			if (ntoken < 16) err = 1;
			else {
				APPEND(p->diel, p->ndiel, cdiel, diel_t);
				diel_t *e = &p->diel[p->ndiel];
				memset(e, 0, sizeof(diel_t));
				for (int k = 0; k < 3; k++) {
					e->org[k] = atof(token[2 + k]);
					e->ea[k] = atof(token[5 + k]);
					e->eb[k] = atof(token[8 + k]);
				}
				e->thick = atof(token[11]);
				e->epsr = atof(token[12]);
				e->ndiva = atoi(token[13]);
				e->ndivb = atoi(token[14]);
				e->ndivt = atoi(token[15]);
				e->tand = (ntoken >= 17) ? atof(token[16]) : 0;
				if (e->tand < 0) err = 1;
				if (ntoken >= 19) {
					e->epsinf = atof(token[17]);
					const int nrest = ntoken - 18;   // eps_inf の後ろの語数
					if (nrest == 1) {
						// 単極 (従来の書き方) : 強度は epsr - eps_inf
						e->npole = 1;
						e->pf[0] = atof(token[18]);
						e->pde[0] = e->epsr - e->epsinf;
					}
					else {
						// 多極 : Debye 極 (f_k de_k) と Lorentz 極
						// (lorentz f_k de_k g_k) の並び。語が余ったら誤り
						int it = 18;
						while (!err && (it < ntoken)) {
							if (e->npole >= DIEL_MAXPOLE) {
								printf("%s\n", "*** dielectric : too many poles");
								return 1;
							}
							// 種類 : lorentz f0 de g / drude fp g / sigma s / (既定) Debye f de
							int kind = POLE_DEBYE;
							if (!strcmp(token[it], "lorentz")) kind = POLE_LORENTZ;
							else if (!strcmp(token[it], "drude")) kind = POLE_DRUDE;
							else if (!strcmp(token[it], "sigma")) kind = POLE_SIGMA;
							if (kind != POLE_DEBYE) it++;
							const int nw = (kind == POLE_LORENTZ) ? 3 : (kind == POLE_SIGMA) ? 1 : 2;
							if ((it + nw) > ntoken) {
								err = 1;
								break;
							}
							const int k = e->npole;
							e->pkind[k] = kind;
							if (kind == POLE_SIGMA) {
								e->pde[k] = atof(token[it]);
							}
							else if (kind == POLE_DRUDE) {
								e->pf[k] = atof(token[it]);
								e->pg[k] = atof(token[it + 1]);
							}
							else {
								e->pf[k] = atof(token[it]);
								e->pde[k] = atof(token[it + 1]);
								if (kind == POLE_LORENTZ) e->pg[k] = atof(token[it + 2]);
							}
							e->npole++;
							it += nw;
						}
					}

					// eps_inf は 1 以上 epsr 以下 (強度 de_k > 0 の和だけ下がる)。
					// 定数 tand と分散の併用は損失の二重計上なので拒否する
					if ((e->epsinf < 1) || (e->epsinf > e->epsr) || (e->tand != 0)) err = 1;
					// 各項の正値性。Drude / 導電率は直流で発散するので直流極限の和に入れない
					double sde = 0;
					for (int k = 0; k < e->npole; k++) {
						const int kind = e->pkind[k];
						if (kind == POLE_SIGMA) {
							if (e->pde[k] <= 0) err = 1;
							continue;
						}
						if (e->pf[k] <= 0) err = 1;
						if (((kind == POLE_LORENTZ) || (kind == POLE_DRUDE)) && (e->pg[k] <= 0)) err = 1;
						if (kind == POLE_DRUDE) continue;
						if (e->pde[k] <= 0) err = 1;
						sde += e->pde[k];
					}
					// 直流極限の整合 : eps_inf + sum de_k = epsr
					// (単極では定義上厳密に成り立つ。多極は書き間違いの番人)
					if (!err && (fabs(sde - (e->epsr - e->epsinf)) > 1e-9 * e->epsr)) {
						printf("%s\n", "*** dielectric : sum of pole strengths must equal epsr - eps_inf");
						return 1;
					}
				}
				const double la = sqrt((e->ea[0] * e->ea[0]) + (e->ea[1] * e->ea[1]) + (e->ea[2] * e->ea[2]));
				const double lb = sqrt((e->eb[0] * e->eb[0]) + (e->eb[1] * e->eb[1]) + (e->eb[2] * e->eb[2]));
				const double ab = (e->ea[0] * e->eb[0]) + (e->ea[1] * e->eb[1]) + (e->ea[2] * e->eb[2]);
				if ((la <= 0) || (lb <= 0) || (e->thick <= 0) || (e->epsr < 1)
				 || (e->ndiva < 1) || (e->ndivb < 1) || (e->ndivt < 1)) err = 1;
				if (!err && (fabs(ab) > 1e-9 * la * lb)) {
					printf("%s\n", "*** dielectric : the two edge vectors must be perpendicular");
					return 1;
				}
				// epsr = 1 は真空と同じ : セルを作らない (完全に無効果、空気と bit 一致)。
				// ただし導電性 (Drude / 導電率) があれば epsr = 1 でも電流が流れる
				int cond = 0;
				for (int k = 0; k < e->npole; k++) {
					if ((e->pkind[k] == POLE_DRUDE) || (e->pkind[k] == POLE_SIGMA)) cond = 1;
				}
				if (!err && ((e->epsr > 1) || cond)) p->ndiel++;
			}
			if (err) {
				printf(errfmt2, "dielectric");
				return 1;
			}
		}
		else if (!strcmp(strkey, "quad")) {
			// quad = x1 y1 z1  x2 y2 z2  x3 y3 z3  x4 y4 z4  厚さ 導電率 分割a 分割b [分割t]
			// 凸な一般四辺形 (頂点は一周順)。双一次写像の構造格子で分割する。
			// 分割 a は辺 1-2 (と 4-3)、分割 b は辺 1-4 (と 2-3) の方向。
			if (ntoken < 18) err = 1;
			else {
				APPEND(p->panel, p->npanel, cpanel, panel_t);
				panel_t *e = &p->panel[p->npanel];
				memset(e, 0, sizeof(panel_t));
				e->kind = PANEL_QUAD;
				for (int k = 0; k < 12; k++) {
					e->v[k] = atof(token[2 + k]);
				}
				e->thick = atof(token[14]);
				e->sigma = atof(token[15]);
				e->ndiva = atoi(token[16]);
				e->ndivb = atoi(token[17]);
				e->ndivt = (ntoken >= 19) ? atoi(token[18]) : 1;
				if ((e->thick <= 0) || (e->sigma < 0)
				 || (e->ndiva < 1) || (e->ndivb < 1) || (e->ndivt < 1)) err = 1;
				// 平面性と凸性 (双一次格子と双対セルの前提)
				if (!err) {
					double nrm[3] = {0, 0, 0};
					double size = 0;
					for (int k = 0; k < 4; k++) {
						const double *a = &e->v[3 * k];
						const double *b = &e->v[3 * ((k + 1) % 4)];
						nrm[0] += (a[1] * b[2]) - (a[2] * b[1]);
						nrm[1] += (a[2] * b[0]) - (a[0] * b[2]);
						nrm[2] += (a[0] * b[1]) - (a[1] * b[0]);
					}
					const double nl = sqrt((nrm[0] * nrm[0]) + (nrm[1] * nrm[1]) + (nrm[2] * nrm[2]));
					if (nl <= 0) err = 1;
					for (int k = 0; !err && (k < 4); k++) {
						const double *a = &e->v[3 * k];
						const double *b = &e->v[3 * ((k + 1) % 4)];
						const double *c = &e->v[3 * ((k + 2) % 4)];
						double u[3], w[3], x[3];
						for (int i = 0; i < 3; i++) {
							u[i] = b[i] - a[i];
							w[i] = c[i] - b[i];
						}
						const double ul = sqrt((u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]));
						if (ul > size) size = ul;
						double d0 = 0;
						for (int i = 0; i < 3; i++) {
							d0 += (b[i] - e->v[i]) * nrm[i] / nl;
						}
						if (fabs(d0) > 1e-9 * (size + fabs(d0))) {
							printf("%s\n", "*** quad : vertices must be coplanar");
							return 1;
						}
						x[0] = (u[1] * w[2]) - (u[2] * w[1]);
						x[1] = (u[2] * w[0]) - (u[0] * w[2]);
						x[2] = (u[0] * w[1]) - (u[1] * w[0]);
						if (((x[0] * nrm[0]) + (x[1] * nrm[1]) + (x[2] * nrm[2])) <= 0) {
							printf("%s\n", "*** quad : must be convex (vertices in cyclic order)");
							return 1;
						}
					}
				}
				if (!err) p->npanel++;
			}
			if (err) {
				printf(errfmt2, "quad");
				return 1;
			}
		}
		else if (!strcmp(strkey, "disk")) {
			// disk = cx cy cz  nx ny nz  半径 厚さ 導電率 nring nsec [分割t]
			// 円板。中心 + nring リング x nsec セクタの極格子で分割する。
			if (ntoken < 13) err = 1;
			else {
				APPEND(p->panel, p->npanel, cpanel, panel_t);
				panel_t *e = &p->panel[p->npanel];
				memset(e, 0, sizeof(panel_t));
				e->kind = PANEL_DISK;
				for (int k = 0; k < 3; k++) {
					e->org[k] = atof(token[2 + k]);
					e->nrm[k] = atof(token[5 + k]);
				}
				e->radius = atof(token[8]);
				e->thick = atof(token[9]);
				e->sigma = atof(token[10]);
				e->ndiva = atoi(token[11]);       // nring
				e->ndivb = atoi(token[12]);       // nsec
				e->ndivt = (ntoken >= 14) ? atoi(token[13]) : 1;
				if (e->ndivt < 1) err = 1;
				const double nl = sqrt((e->nrm[0] * e->nrm[0])
					+ (e->nrm[1] * e->nrm[1]) + (e->nrm[2] * e->nrm[2]));
				// 中心の双対環は 2 x nsec 頂点なので nsec <= POLY_MAX/2
				if ((nl <= 0) || (e->radius <= 0) || (e->thick <= 0) || (e->sigma < 0)
				 || (e->ndiva < 1) || (e->ndivb < 3) || (e->ndivb > POLY_MAX / 2)) err = 1;
				if (!err) p->npanel++;
			}
			if (err) {
				printf(errfmt2, "disk");
				return 1;
			}
		}
		else if (!strcmp(strkey, "port")) {
			if (ntoken < 5) err = 1;
			else {
				APPEND(p->port, p->nport, cport, port_t);
				port_t *e = &p->port[p->nport];
				e->n1 = nodeid(token[2], &err);
				e->n2 = nodeid(token[3], &err);
				e->z0 = atof(token[4]);
				if (e->n1 == e->n2) err = 1;
				updatemax(p, e->n1);
				updatemax(p, e->n2);
				if (!err) p->nport++;
			}
			if (err) {
				printf(errfmt2, "port");
				return 1;
			}
		}
		else if (!strcmp(strkey, "mixedmode")) {
			// mixedmode = p1 n1 [p2 n2 ...] : 単端ポート番号 (1 始まり) の対
			// (正側, 負側) を差動対にして混合モード S を出す (mixed.csv)。
			// ポート番号の範囲・重複・対内の基準抵抗の一致はポートが
			// 出揃ってから (ファイル末尾で) 検査する
			const int nw = ntoken - 2;
			if ((nw < 2) || ((nw % 2) != 0) || (nw > 2 * MM_MAXPAIR)) err = 1;
			else {
				p->nmm = nw / 2;
				for (int k = 0; k < nw; k++) p->mmport[k] = atoi(token[2 + k]) - 1;
			}
			if (err) {
				printf(errfmt2, "mixedmode");
				return 1;
			}
		}
		else if (!strcmp(strkey, "frequency")) {
			if (ntoken < 5) err = 1;
			else {
				p->f0 = atof(token[2]);
				p->f1 = atof(token[3]);
				const int ndiv = atoi(token[4]);
				p->nfreq = ndiv + 1;
				// 第 4 引数 "log" で等比 (対数) 掃引 (省略時は従来どおり線形)
				if ((ntoken > 5) && !strcmp(token[5], "log")) p->flog = 1;
				if ((p->f0 <= 0) || (p->f1 < p->f0) || (ndiv < 0)) err = 1;
			}
			if (err) {
				printf(errfmt2, "frequency");
				return 1;
			}
		}
		else if (!strcmp(strkey, "nodetol")) {
			p->nodetol = atof(token[2]);
			if (p->nodetol <= 0) {
				printf(errfmt2, "nodetol");
				return 1;
			}
		}
		else if (!strcmp(strkey, "skineffect")) {
			p->skin = atoi(token[2]);
		}
		else if (!strcmp(strkey, "capacitance")) {
			p->capacitance = atoi(token[2]);
		}
		else if (!strcmp(strkey, "distribution")) {
			p->dist = atoi(token[2]);
		}
		else if (!strcmp(strkey, "retardation")) {
			p->retardation = atoi(token[2]);
		}
		else if (!strcmp(strkey, "acceleration")) {
			p->accel = atoi(token[2]);
		}
		else if (!strcmp(strkey, "transient")) {
			// transient = 1 [帯域端の減衰dB] : 掃引の逆フーリエ変換で時間波形
			p->tran = atoi(token[2]);
			if (ntoken >= 4) p->tranatt = atof(token[3]);
			if (p->tranatt <= 0) {
				printf(errfmt2, "transient");
				return 1;
			}
		}
		else if (!strcmp(strkey, "planewave")) {
			// planewave = theta phi 偏波 振幅 [位相deg]
			// (theta, phi) は到来方向 (伝搬方向はその逆)、偏波 1 = theta / 2 = phi。
			// 複数行は同時に入射する波 (重ね合わせ、PW_MAX 行まで)。1 行目が
			// leff (pw.csv) と RCS (farpw.csv) の振幅の規格化に使われる
			if ((ntoken < 6) || (p->npw >= PW_MAX)) err = 1;
			else {
				const int k = p->npw;
				p->pwi_th[k] = atof(token[2]);
				p->pwi_ph[k] = atof(token[3]);
				p->pwi_pol[k] = atoi(token[4]);
				p->pwi_amp[k] = atof(token[5]);
				p->pwi_phase[k] = (ntoken >= 7) ? atof(token[6]) : 0;
				if ((p->pwi_pol[k] != 1) && (p->pwi_pol[k] != 2)) err = 1;
				if (p->pwi_amp[k] <= 0) err = 1;
				if (!err) {
					if (k == 0) {
						p->pwth = p->pwi_th[0];
						p->pwph = p->pwi_ph[0];
						p->pwpol = p->pwi_pol[0];
						p->pwamp = p->pwi_amp[0];
						p->pwphase = p->pwi_phase[0];
					}
					p->npw++;
					p->pw = 1;
				}
			}
			if (err) {
				printf(errfmt2, "planewave");
				return 1;
			}
		}
		else if (!strcmp(strkey, "farfield")) {
			if (ntoken < 4) err = 1;
			else {
				p->ffnth = atoi(token[2]);
				p->ffnph = atoi(token[3]);
				if ((p->ffnth < 2) || (p->ffnph < 2)) err = 1;
			}
			if (err) {
				printf(errfmt2, "farfield");
				return 1;
			}
		}
		else if (!strcmp(strkey, "nearfield")) {
			// nearfield = x1 y1 z1 x2 y2 z2 n1 n2 n3
			//   角 x1 -> 角 x2 の直方体を各軸 n 分割した格子点 (n = 0 は固定)
			if (ntoken < 11) err = 1;
			else {
				APPEND(p->nf, p->nnf, cnf, nf_t);
				nf_t *e = &p->nf[p->nnf];
				memset(e, 0, sizeof(nf_t));
				for (int k = 0; k < 3; k++) {
					e->p1[k] = atof(token[2 + k]);
					e->p2[k] = atof(token[5 + k]);
					e->n[k] = atoi(token[8 + k]);
					if (e->n[k] < 0) err = 1;
				}
				if (!err) p->nnf++;
			}
			if (err) {
				printf(errfmt2, "nearfield");
				return 1;
			}
		}
		else if (!strcmp(strkey, "grading")) {
			p->grading = atoi(token[2]);
		}
		else if (!strcmp(strkey, "spice")) {
			// spice = 1 [kmin] : 準静的な等価回路を peec.cir に書く (|k| <= kmin の結合は省く)
			p->spice = atoi(token[2]);
			p->spkmin = (ntoken >= 4) ? atof(token[3]) : 0;
			if (((p->spice != 0) && (p->spice != 1)) || (p->spkmin < 0) || (p->spkmin >= 1)) {
				printf(errfmt2, "spice");
				return 1;
			}
		}
		else if (!strcmp(strkey, "basis")) {
			// basis = 0 (区間ごとに一定の電流、既定) | 1 (区分線形の電流、細線のみ)
			p->basis = atoi(token[2]);
			if ((p->basis != 0) && (p->basis != 1)) {
				printf(errfmt2, "basis");
				return 1;
			}
		}
		else if (!strcmp(strkey, "compression")) {
			// compression = 1 [tol] : 部分インダクタンスを H 行列で圧縮
			p->compress = atoi(token[2]);
			if (ntoken >= 4) p->ctol = atof(token[3]);
			if (p->ctol <= 0) {
				printf(errfmt2, "compression");
				return 1;
			}
		}
		else if (!strcmp(strkey, "groundplane")) {
			// groundplane = z0            (z = z0 の面 : 従来の書き方)
			// groundplane = x|y|z c       (軸に垂直な面。互いに直交する 3 枚まで)
			int ax = 2;
			double c;
			if (!strcmp(token[2], "x") || !strcmp(token[2], "y") || !strcmp(token[2], "z")) {
				ax = token[2][0] - 'x';
				if (ntoken < 4) err = 1;
				c = (ntoken >= 4) ? atof(token[3]) : 0;
			}
			else {
				c = atof(token[2]);
			}
			// 平行な 2 面 (同じ軸) は鏡像が無限級数になるので扱わない
			for (int k = 0; k < p->ngp; k++) {
				if (p->gpax[k] == ax) err = 1;
			}
			if (p->ngp >= 3) err = 1;
			if (err) {
				printf(errfmt2, "groundplane");
				return 1;
			}
			p->gpax[p->ngp] = ax;
			p->gpc[p->ngp] = c;
			p->ngp++;
			p->gp = 1;
		}
		else if (!strcmp(strkey, "gmin")) {
			p->gmin = atof(token[2]);
			if (p->gmin < 0) {
				printf(errfmt2, "gmin");
				return 1;
			}
		}
		// 未知キーは無視する (前方互換)

		nline++;
	}

	// mutual の名前を inductor index に解決
	for (int i = 0; i < p->nmut; i++) {
		mut_t *e = &p->mut[i];
		for (int j = 0; j < p->nind; j++) {
			if (!strcmp(e->name1, p->ind[j].name)) e->l1 = j;
			if (!strcmp(e->name2, p->ind[j].name)) e->l2 = j;
		}
		if ((e->l1 < 0) || (e->l2 < 0) || (e->l1 == e->l2)) {
			printf("*** mutual : unknown or same inductor : %s %s\n", e->name1, e->name2);
			return 1;
		}
	}

	// 必須キーの確認
	if (p->nport <= 0) {
		printf("%s\n", "*** no port data");
		return 1;
	}
	if (p->nfreq <= 0) {
		printf("%s\n", "*** no frequency data");
		return 1;
	}
	// mixedmode : 各ポートは 1 つの対にしか入れない。差動 / 同相の基準
	// (2 z0 / z0/2) は対の 2 ポートの z0 が等しいときだけ定義できる
	for (int k = 0; k < 2 * p->nmm; k++) {
		const int ip = p->mmport[k];
		int bad = (ip < 0) || (ip >= p->nport);
		for (int l = 0; !bad && (l < k); l++) {
			if (p->mmport[l] == ip) bad = 1;
		}
		if (!bad && ((k % 2) == 1) && (p->port[ip].z0 != p->port[p->mmport[k - 1]].z0)) {
			printf("*** mixedmode : ports %d and %d must have the same reference impedance\n",
				p->mmport[k - 1] + 1, ip + 1);
			return 1;
		}
		if (bad) {
			printf("*** mixedmode : invalid or duplicated port number %d\n", ip + 1);
			return 1;
		}
	}
	// excitation : 範囲外・重複のポートは拒否 (重ね合わせの重みが曖昧になる)
	for (int k = 0; k < p->nexc; k++) {
		int bad = (p->exc[k].port >= p->nport);
		for (int l = 0; !bad && (l < k); l++) {
			if (p->exc[l].port == p->exc[k].port) bad = 1;
		}
		if (bad) {
			printf("*** excitation : invalid or duplicated port number %d\n", p->exc[k].port + 1);
			return 1;
		}
	}
	if ((p->nres + p->ncap + p->nind + p->ntl + p->nwire + p->nplate + p->npanel) <= 0) {
		printf("%s\n", "*** no element data (resistor/capacitor/inductor/tline/wire/bar/plate/polygon/disk)");
		return 1;
	}
	// 誘電体は電位係数 (束縛電荷) が無いと意味を持たない
	if ((p->ndiel > 0) && !p->capacitance) {
		printf("%s\n", "*** dielectric requires capacitance = 1");
		return 1;
	}
	// acceleration は密 LU の再利用が前提なので圧縮とは併用できない
	if (p->compress && p->accel) {
		printf("%s\n", "*** compression = 1 and acceleration = 1 are mutually exclusive");
		return 1;
	}
	if (p->compress && ((p->nwire + p->nplate + p->npanel) <= 0)) {
		printf("%s\n", "*** compression = 1 requires conductor geometry (wire/bar/plate/quad/disk)");
		return 1;
	}
	// spice = 1 は準静的な R / L / C で書ける構成だけ : 密な Lp / C が要る (圧縮不可)、
	// basis = 1 の傾斜関数・誘電損と分散・損失のある線路は周波数依存で書けない
	if (p->spice) {
		const char *why = NULL;
		if (p->compress) why = "compression = 0";
		else if (p->basis) why = "basis = 0";
		for (int i = 0; (why == NULL) && (i < p->ndiel); i++) {
			if ((p->diel[i].tand != 0) || (p->diel[i].npole > 0)) why = "lossless non-dispersive dielectric";
		}
		for (int i = 0; (why == NULL) && (i < p->ntl); i++) {
			if (p->tl[i].rlgc || (p->tl[i].att != 1)) why = "lossless tline (Z0 / delay)";
		}
		if (why != NULL) {
			printf("*** spice = 1 requires %s\n", why);
			return 1;
		}
	}
	// basis = 1 は細線 (分割しない wire / bar) の電荷で駆動する傾斜関数なので、
	// 容量性 PEEC が必須。面導体・体積セル・誘電体・圧縮は扱わない
	if (p->basis) {
		const char *why = NULL;
		if (!p->capacitance) why = "capacitance = 1";
		else if (p->nwire <= 0) why = "wire or bar conductors";
		else if ((p->nplate + p->npanel) > 0) why = "thin wires only (no plate/quad/disk)";
		else if (p->ndiel > 0) why = "no dielectric";
		else if (p->compress) why = "compression = 0";
		for (int i = 0; (why == NULL) && (i < p->nwire); i++) {
			if ((p->wire[i].ndivw > 1) || (p->wire[i].ndivt > 1)) why = "no cross-section subdivision";
		}
		if (why != NULL) {
			printf("*** basis = 1 requires %s\n", why);
			return 1;
		}
	}

	return 0;
}


// peec_t が抱える動的配列をまとめて解放する
//
// 確保箇所はファイルをまたぐ (input_data / wire_build / lp_fill / pot_fill /
// mna_numbering / solve) が、いずれも peec_t のメンバとして最後まで保持される
// ので解放も 1 箇所にまとめる。input_data() 冒頭の memset で全ポインタが NULL
// から始まるため、途中で失敗した状態から呼んでも安全 (free(NULL) は無害)。
//
// バッチ実行ではプロセス終了時にまとめて回収されるが、密行列 (lp / cmat) が
// 支配的で規模に比例して増えるため、LeakSanitizer を CI に掛けられるように
// 明示的に解放する。
void peec_free(peec_t *p)
{
	if (p == NULL) return;

	// ネットリスト (input_data.c)
	free(p->plate);
	free(p->panel);
	free(p->diel);
	free(p->res);
	free(p->cap);
	free(p->ind);
	free(p->mut);
	free(p->src);
	free(p->tl);
	free(p->exc);
	free(p->port);
	free(p->wire);
	free(p->helix);
	free(p->ncid);
	free(p->ncxyz);
	free(p->nf);

	// 形状・分割 (wire.c)
	free(p->seg);
	free(p->chg);
	free(p->chgnode);
	free(p->gid);
	free(p->gxyz);

	// 部分要素 (partial.c / potential.c / hmatrix.c)
	hmat_free(p->hlp);
	free(p->lp);
	free(p->cellid);
	free(p->cellof);
	free(p->carea);
	free(p->csoff);
	free(p->csidx);
	free(p->cpt);
	free(p->cmat);
	free(p->bk);
	free(p->bx);

	// MNA と結果 (mna.c / solve.c)
	free(p->nodemap);
	free(p->zin);
	free(p->zmat);
	free(p->smat);
	free(p->segi);
	free(p->cellq);
	free(p->voc);
	free(p->segipw);
	free(p->cellqpw);

	memset(p, 0, sizeof(peec_t));
}
