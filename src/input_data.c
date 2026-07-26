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
		else if (!strcmp(strkey, "wire")) {
			if (ntoken < 11) err = 1;
			else {
				APPEND(p->wire, p->nwire, cwire, wire_t);
				wire_t *e = &p->wire[p->nwire];
				for (int k = 0; k < 3; k++) {
					e->x1[k] = atof(token[2 + k]);
					e->x2[k] = atof(token[5 + k]);
				}
				e->radius = atof(token[8]);
				e->sigma = atof(token[9]);
				e->ndiv = atoi(token[10]);
				const double dx = e->x2[0] - e->x1[0];
				const double dy = e->x2[1] - e->x1[1];
				const double dz = e->x2[2] - e->x1[2];
				const double len = sqrt((dx * dx) + (dy * dy) + (dz * dz));
				if ((e->radius <= 0) || (e->sigma < 0) || (e->ndiv < 1) || (len <= 0)) err = 1;
				if (!err) p->nwire++;
			}
			if (err) {
				printf(errfmt2, "wire");
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
	if ((p->nres + p->ncap + p->nind + p->nwire) <= 0) {
		printf("%s\n", "*** no element data (resistor/capacitor/inductor/wire)");
		return 1;
	}

	return 0;
}
