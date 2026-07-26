/*
peec.h

OpenPEEC : 準静的 PEEC (部分要素等価回路) 回路ソルバー
定数・構造体・全プロトタイプ
*/

#ifndef _PEEC_H_
#define _PEEC_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "complex.h"

#define PROGRAM "OpenPEEC"
#define VERSION_MAJOR 1
#define VERSION_MINOR 0
#define VERSION_BUILD 0

#define FN_LOG "peec.log"
#define FN_CSV "zin.csv"

// 数学・物理定数 (自前マクロ : <math.h> の M_PI には依存しない)
#define PI   (4.0 * atan(1.0))
#define C0   (2.99792458e8)
#define MU0  (4.0 * PI * 1e-7)
#define EPS0 (1.0 / (C0 * C0 * MU0))
#define EPS  (1e-12)

#define NAMELEN 32

// 集中定数素子 (ノード id は非負整数、0 = GND)
typedef struct {int n1, n2; double val;} rc_t;                                // resistor / capacitor
typedef struct {char name[NAMELEN]; int n1, n2; double val;} ind_t;           // inductor
typedef struct {char name1[NAMELEN], name2[NAMELEN]; int l1, l2; double k;} mut_t;  // M = k*sqrt(L1*L2)
typedef struct {int n1, n2, isvsrc, ibr; double amp, phase;} src_t;           // phase [deg]
typedef struct {int n1, n2; double z0;} port_t;
typedef struct {double x1[3], x2[3], radius, sigma; int ndiv;} wire_t;

// ワイヤ分割後の 1 区間 = PEEC の 1 枝 (向き n1 -> n2)
typedef struct {
	int    n1, n2;
	double x1[3], x2[3];
	double len, radius, res;      // res : DC 抵抗 = len / (sigma * pi * radius^2)
} seg_t;

typedef struct {
	char   title[BUFSIZ];

	// ネットリスト
	int    nres, ncap, nind, nmut, nsrc, nport, nwire, nnodexyz;
	rc_t   *res, *cap;
	ind_t  *ind;
	mut_t  *mut;
	src_t  *src;
	port_t *port;
	wire_t *wire;
	int    *ncid;                 // node = : ノード id
	double *ncxyz;                // node = : 座標 (フラット [3*i])

	// ワイヤ形状から導出
	int    nseg;
	seg_t  *seg;
	double *lp;                   // 部分インダクタンス密行列 (フラット [i*nseg+j])

	// 形状ノード表 (座標マージ用)
	int    ngnode;
	int    *gid;
	double *gxyz;                 // フラット [3*i]

	// MNA 番号付け
	int    maxid;                 // 使用された最大ノード id (自動採番含む)
	int    nnode;                 // 基準ノードを除く行列ノード数
	int    refnode;               // 基準ノード id
	int    *nodemap;              // [maxid+1] : id -> 行列 index / -1 = 基準 / -2 = 未使用
	int    offL, offV, offS;      // 枝電流の先頭 index (L 素子 / 電圧源 / ワイヤ区間)
	int    nunknown;

	// 掃引・オプション
	double f0, f1;
	int    nfreq;
	double nodetol, gmin;

	// 結果
	d_complex_t *zin;             // [nport * nfreq]
} peec_t;

// 掃引周波数 (ifreq = 0 ... nfreq-1)
static inline double freq_at(const peec_t *p, int ifreq)
{
	return (p->nfreq <= 1) ? p->f0
	     : p->f0 + (p->f1 - p->f0) * ifreq / (p->nfreq - 1);
}

// utils.c
int  tokenize(char *str, const char *tokensep, char *token[], int maxtoken);

// input_data.c
int  input_data(FILE *fp, peec_t *p);

// wire.c
int  wire_build(peec_t *p, FILE *fp_log);

// partial.c
double lp_self(double l, double a);
double lp_pair(const seg_t *s1, const seg_t *s2);
void lp_fill(peec_t *p, FILE *fp_log);

// mna.c
int  mna_numbering(peec_t *p, FILE *fp_log);
void mna_assemble(const peec_t *p, double f, d_complex_t *a);
void mna_rhs_port(const peec_t *p, int iport, d_complex_t *b);
void mna_rhs_sources(const peec_t *p, d_complex_t *b);

// lu.c
int  lu_decomp(int n, d_complex_t *a, int *piv);
void lu_solve(int n, const d_complex_t *a, const int *piv, d_complex_t *b);

// solve.c
int  solve(peec_t *p, FILE *fp_log);

// output.c
void output_zin(const peec_t *p, FILE *fp_log);
int  output_csv(const peec_t *p, const char *fn);

#endif		// _PEEC_H_
