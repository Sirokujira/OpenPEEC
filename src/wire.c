/*
wire.c

ワイヤの区間分割と座標ノードマージ
- node = で束縛された座標を種として、ワイヤ端点/分割点を nodetol 以内で照合する
- 一致しない点には maxid+1 から決定的に自動採番する (ファイル記載順)
*/

#include "peec.h"

static double dist3(const double *a, const double *b)
{
	const double dx = a[0] - b[0];
	const double dy = a[1] - b[1];
	const double dz = a[2] - b[2];
	return sqrt((dx * dx) + (dy * dy) + (dz * dz));
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

int wire_build(peec_t *p, FILE *fp_log)
{
	if (p->nwire <= 0) return 0;

	// ノード表の最大数 = node 束縛数 + 全分割点数
	int npoint = p->nnodexyz;
	int nseg = 0;
	for (int i = 0; i < p->nwire; i++) {
		npoint += p->wire[i].ndiv + 1;
		nseg += p->wire[i].ndiv;
	}
	p->gid = (int *)malloc((size_t)npoint * sizeof(int));
	p->gxyz = (double *)malloc((size_t)npoint * 3 * sizeof(double));
	p->seg = (seg_t *)malloc((size_t)nseg * sizeof(seg_t));
	p->ngnode = 0;
	p->nseg = 0;

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

	// ワイヤ分割
	int autoid = p->maxid;
	for (int i = 0; i < p->nwire; i++) {
		const wire_t *w = &p->wire[i];
		const int ndiv = w->ndiv;
		int idprev = -1;
		for (int k = 0; k <= ndiv; k++) {
			double x[3];
			for (int c = 0; c < 3; c++) {
				x[c] = w->x1[c] + (w->x2[c] - w->x1[c]) * k / ndiv;
			}
			int id;
			const int found = find_gnode(p, x);
			if (found >= 0) {
				id = p->gid[found];
			}
			else {
				id = ++autoid;
				add_gnode(p, id, x);
			}
			if (k > 0) {
				seg_t *s = &p->seg[p->nseg];
				s->n1 = idprev;
				s->n2 = id;
				for (int c = 0; c < 3; c++) {
					s->x1[c] = w->x1[c] + (w->x2[c] - w->x1[c]) * (k - 1) / ndiv;
					s->x2[c] = x[c];
				}
				s->len = dist3(s->x1, s->x2);
				s->radius = w->radius;
				s->sigma = w->sigma;
				s->res = (w->sigma > 0) ? s->len / (w->sigma * PI * w->radius * w->radius) : 0;
				p->nseg++;
			}
			idprev = id;
		}
	}
	p->maxid = autoid;

	// 近接ノード警告 (マージされなかったが nodetol の 10 倍未満 : 意図しない分離の可能性)
	for (int i = 0; i < p->ngnode; i++) {
		for (int j = i + 1; j < p->ngnode; j++) {
			const double d = dist3(&p->gxyz[3 * i], &p->gxyz[3 * j]);
			if (d < 10 * p->nodetol) {
				fprintf(fp_log, "*** warning : node %d and node %d are separated by only %.3e m\n",
					p->gid[i], p->gid[j], d);
			}
		}
	}

	fprintf(fp_log, "wire : %d wires -> %d segments, %d geometry nodes\n",
		p->nwire, p->nseg, p->ngnode);

	return 0;
}
