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

int input_data(FILE *fp, peec_t *p)
{
	int    nline = 0;
	char   strline[BUFSIZ], strkey[BUFSIZ], strsave[BUFSIZ];
	char   *token[MAXTOKEN];
	const char sep[] = " \t";
	const char errfmt2[] = "*** invalid %s data\n";
	int    cres = 0, ccap = 0, cind = 0, cmut = 0, csrc = 0, cport = 0, cwire = 0, cnode = 0;
	int    cplate = 0, cpanel = 0;

	// initialize (既定値 : キー省略時は従来動作)
	memset(p, 0, sizeof(peec_t));
	p->nodetol = 1e-8;
	p->gmin = 0;
	p->refnode = -1;

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
				if (e->val <= 0) err = 1;
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
				if ((e->k <= 0) || (e->k > 1)) err = 1;
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
			// wire = x1 y1 z1 x2 y2 z2 半径     導電率 分割数
			// bar  = x1 y1 z1 x2 y2 z2 幅 厚さ  導電率 分割数
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
					if ((e->width <= 0) || (e->thick <= 0)) err = 1;
				}
				else {
					e->shape = SHAPE_ROUND;
					e->radius = atof(token[8]);
					e->sigma = atof(token[9]);
					e->ndiv = atoi(token[10]);
					if (e->radius <= 0) err = 1;
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
		else if (!strcmp(strkey, "quad")) {
			// quad = x1 y1 z1  x2 y2 z2  x3 y3 z3  x4 y4 z4  厚さ 導電率 分割a 分割b
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
				if ((e->thick <= 0) || (e->sigma < 0)
				 || (e->ndiva < 1) || (e->ndivb < 1)) err = 1;
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
			// disk = cx cy cz  nx ny nz  半径 厚さ 導電率 nring nsec
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
		else if (!strcmp(strkey, "frequency")) {
			if (ntoken < 5) err = 1;
			else {
				p->f0 = atof(token[2]);
				p->f1 = atof(token[3]);
				const int ndiv = atoi(token[4]);
				p->nfreq = ndiv + 1;
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
	if ((p->nres + p->ncap + p->nind + p->nwire + p->nplate + p->npanel) <= 0) {
		printf("%s\n", "*** no element data (resistor/capacitor/inductor/wire/bar/plate/polygon/disk)");
		return 1;
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
	free(p->res);
	free(p->cap);
	free(p->ind);
	free(p->mut);
	free(p->src);
	free(p->port);
	free(p->wire);
	free(p->ncid);
	free(p->ncxyz);

	// 形状・分割 (wire.c)
	free(p->seg);
	free(p->chg);
	free(p->chgnode);
	free(p->gid);
	free(p->gxyz);

	// 部分要素 (partial.c / potential.c)
	free(p->lp);
	free(p->cellid);
	free(p->cellof);
	free(p->carea);
	free(p->cmat);

	// MNA と結果 (mna.c / solve.c)
	free(p->nodemap);
	free(p->zin);
	free(p->zmat);
	free(p->smat);
	free(p->segi);
	free(p->cellq);

	memset(p, 0, sizeof(peec_t));
}
