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
#define FN_DIST "dist.csv"
#define FN_FAR "far.csv"
#define FN_PW "pw.csv"

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
// 導体断面
//   SHAPE_ROUND = 丸線 (半径 a)
//   SHAPE_BAR   = 角線 (幅 w x 厚さ t)  — 細線 (等価半径) として扱う
//   SHAPE_PLATE = 面導体のセル (幅 wid の帯 = リボン) — 面積分で扱う
//   SHAPE_POLY  = 三角形メッシュ由来の平面多角形セル — 多角形面積分で扱う
#define SHAPE_ROUND 0
#define SHAPE_BAR   1
#define SHAPE_PLATE 2
#define SHAPE_POLY  3

// 多角形セルの最大頂点数 (disk 中心の双対セルは 2 x nsec なので nsec <= 32)
#define POLY_MAX 64

typedef struct {
	double x1[3], x2[3];
	int    shape;
	double radius, width, thick, sigma;
	int    ndiv;
} wire_t;

// 構造格子でメッシュ化する平面パネル (quad / disk)。
// plate (直交矩形格子) の一般化 : セルは平面多角形になり、多角形の
// 幾何二重積分 (polygon.c) で結合する。
//   PANEL_QUAD : 凸な一般四辺形 (双一次写像の ndiva x ndivb 格子)
//   PANEL_DISK : 円板 (中心 + nring リング x nsec セクタの極格子)
#define PANEL_QUAD 0
#define PANEL_DISK 1
typedef struct {
	int    kind;
	double v[12];                 // quad の 4 頂点 (一周順)
	double org[3], nrm[3];        // disk の中心と法線
	double radius;
	int    ndiva, ndivb;          // quad の分割数 / disk は (nring, nsec)
	double thick, sigma;
} panel_t;

// 平面矩形導体 : o + s*ea + t*eb (s, t は 0..1)
// ndivt >= 2 で厚み方向にも分割し、セルを体積バー (Hoer-Love) として扱う。
// ndivt = 1 (既定) は従来どおりリボン (面積分) 1 層。
typedef struct {
	double org[3], ea[3], eb[3];
	double thick, sigma;
	int    ndiva, ndivb, ndivt;
} plate_t;

// 誘電体ブリック (Ruehli の過剰容量による誘電体 PEEC)。
// org を底面の角として o + s*ea + t*eb + u*thick*n (n = ea x eb 方向、
// s,t,u は 0..1) の直方体を占める (plate と違い法線方向へ片側に押し出す)。
// 3 方向の枝 (体積セル) に過剰容量 C_e = eps0 (epsr-1) A/len を直列に置き、
// 節点の束縛電荷が電位係数 P に参加する。導体と接する面のノードは
// nodetol マージで導体ノードと共有される。
typedef struct {
	double org[3], ea[3], eb[3];
	double thick, epsr;
	int    ndiva, ndivb, ndivt;
} diel_t;

// 分割後の 1 セル。導体電流 (枝) にも電荷セルにも使う。
// wid = 0 なら細線フィラメント、wid > 0 なら幅 wid の帯 (リボン)。
// リボンは x1->x2 を軸、wv を面内の横方向単位ベクトルとする矩形。
typedef struct {
	int    n1, n2;
	double x1[3], x2[3];
	double len;
	int    shape;
	double radius, width, thick;
	double wid, wv[3];            // リボンの幅と横方向単位ベクトル
	int    vol;                   // 1 = 体積セル (断面 wid x thick の矩形バー)
	// 多角形セル (npv > 0 のとき有効) : 平面多角形の頂点リング pv と
	// 扇分割の基点 papex (セルはこの点から星形)。wid = 面積/len。
	int    npv;
	double papex[3];
	double pv[3 * POLY_MAX];
	double aL, aP;                // 細線の等価半径 : インダクタンス用 (GMD) / 容量用
	double area, perim;           // 断面積・周長 (DC 抵抗・表皮効果)
	double sigma;
	double res;                   // DC 抵抗 = len / (sigma * area)
	int    diel;                  // 1 = 誘電体セル (枝の直列インピーダンス = 1/(jw cexc))
	double cexc;                  // 過剰容量 C_e = eps0 (epsr-1) area / len [F]
} seg_t;

typedef struct {
	char   title[BUFSIZ];

	// ネットリスト
	int    nres, ncap, nind, nmut, nsrc, nport, nwire, nplate, npanel, ndiel, nnodexyz;
	plate_t *plate;
	panel_t *panel;
	diel_t *diel;
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
	d_complex_t *lp;              // 部分インダクタンス密行列 (フラット [i*nseg+j])
	int    lpwarn;                // |k| > 1 警告を出したか

	// 容量性 PEEC (capacitance = 1 のときのみ)
	// 電荷セル (幾何段で作る : ワイヤは半区間、面導体は節点の双対矩形)
	int    nchg;                  // 電荷セルの個数
	seg_t  *chg;                  // [nchg] 電荷セル
	int    *chgnode;              // [nchg] 電荷セルが属するノード id
	int    ncell;                 // 容量セル (= 幾何ノード) 数
	int    *cellid;               // [ncell] セルのノード id
	int    *cellof;               // [nchg] 電荷セル -> セル index
	double *carea;                // [ncell] セルの長さ (細線) / 面積 (面導体)
	d_complex_t *cmat;            // [ncell*ncell] 節点容量行列 C = P^-1
	int    clogged;               // 総容量をログ出力したか

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
	int    flog;                  // frequency ... log : 対数 (等比) 掃引
	double nodetol, gmin;
	int    skin;                  // skineffect = 1 : 表皮効果 + 内部インダクタンス
	int    capacitance;           // capacitance = 1 : 容量性 PEEC (電位係数)
	int    retardation;           // retardation = 1 : 遅延 (フルウェーブ PEEC)
	int    accel;                 // acceleration = 1 : 掃引で LU を再利用する GMRES
	int    dist;                  // distribution = 1 : 電流・電荷分布を出力
	int    gp;                    // groundplane = z : 無限 PEC 地板 (鏡像法)
	double gpz;                   // 地板の z 座標 (gp = 1 のとき有効)
	int    ffnth, ffnph;          // farfield = 分割数theta 分割数phi (0 = 無効)
	int    pw;                    // planewave : 平面波入射 (外部界励振)
	double pwth, pwph;            // 平面波の到来方向 (theta, phi) [deg]
	int    pwpol;                 // 偏波 : 1 = theta 偏波、2 = phi 偏波
	double pwamp, pwphase;        // 振幅 [V/m] と位相 [deg]

	// 結果
	d_complex_t *zin;             // [nport * nfreq] 各ポートの入力インピーダンス
	// 多ポートの Z 行列と S 行列 (どちらも [ifreq][i][j] = [(ifreq*nport + i)*nport + j])
	// Z はポート j に 1A を注入し (他ポートは開放) ポート i の電圧を読んで得る。
	// S は電力波の定義 (Kurokawa) でポートごとの実数基準抵抗 z0 から変換する。
	d_complex_t *zmat;            // [nfreq * nport * nport]
	d_complex_t *smat;            // [nfreq * nport * nport]
	// 分布 (distribution = 1 のときのみ)。port #1 を 1A で励振したときの値。
	d_complex_t *segi;            // [nfreq * nseg]  区間電流 [A]
	d_complex_t *cellq;           // [nfreq * ncell] セル電荷 [C] (capacitance = 1 のときのみ)
	// 平面波入射 (planewave) : 各ポートの端子電圧 [V] と誘起電流 [A]
	// (ポート間に素子が無ければ端子電圧 = 開放端電圧 Voc)
	d_complex_t *voc;             // [nport * nfreq]
	d_complex_t *segipw;          // [nfreq * nseg] 誘起区間電流 (distribution = 1)
} peec_t;

// Z / S 行列の添字 (周波数 ifreq、行 i、列 j)
#define ZIDX(p, ifreq, i, j) ((size_t)((ifreq) * (p)->nport + (i)) * (p)->nport + (j))

// 掃引周波数 (ifreq = 0 ... nfreq-1)。flog = 1 なら等比 (対数) 掃引。
static inline double freq_at(const peec_t *p, int ifreq)
{
	if (p->nfreq <= 1) return p->f0;
	if (p->flog) {
		return p->f0 * pow(p->f1 / p->f0, (double)ifreq / (p->nfreq - 1));
	}
	return p->f0 + (p->f1 - p->f0) * ifreq / (p->nfreq - 1);
}

// utils.c
int  tokenize(char *str, const char *tokensep, char *token[], int maxtoken);

// input_data.c
int  input_data(FILE *fp, peec_t *p);
void peec_free(peec_t *p);   // peec_t の動的配列を一括解放する

// wire.c
int  wire_build(peec_t *p, FILE *fp_log);

// surface.c
double rect_potential(const seg_t *s, const double *pt);
double ribbon_static(const seg_t *s1, const seg_t *s2, int nsub);
d_complex_t ribbon_corr(const seg_t *s1, const seg_t *s2, double kw, int nsub);

// volume.c
int    bar_use_hl(const seg_t *s1, const seg_t *s2);
double bar_pair(const seg_t *s1, const seg_t *s2);

// polygon.c
double poly_area(const seg_t *s);
double poly_potential(const seg_t *s, const double *pt);
double poly_static(const seg_t *s1, const seg_t *s2, int nsub);
d_complex_t poly_corr(const seg_t *s1, const seg_t *s2, double kw);

// partial.c
void seg_mirror(const seg_t *s, double gpz, seg_t *out);
double neumann_self(double l, double a);
double neumann_pair(const seg_t *s1, const seg_t *s2, double a1, double a2);
d_complex_t neumann_self_k(const seg_t *s, double a, double kw);
d_complex_t neumann_pair_k(const seg_t *s1, const seg_t *s2, double a1, double a2, double kw);
double lp_self(double l, double a);
void lp_fill(peec_t *p, double f, FILE *fp_log);

// potential.c
int  pot_fill(peec_t *p, double f, FILE *fp_log);

// skin.c
d_complex_t zint_round(double len, double a, double sigma, double freq);
d_complex_t zint_bar(double len, double area, double perim, double sigma, double freq);
d_complex_t zint_seg(const seg_t *s, double freq);

// mna.c
int  mna_numbering(peec_t *p, FILE *fp_log);
void mna_assemble(const peec_t *p, double f, d_complex_t *a);
void mna_rhs_port(const peec_t *p, int iport, d_complex_t *b);
void mna_rhs_sources(const peec_t *p, d_complex_t *b);
void mna_rhs_planewave(const peec_t *p, double f, d_complex_t *b);

// lu.c
int  lu_decomp(int n, d_complex_t *a, int *piv);
void lu_solve(int n, const d_complex_t *a, const int *piv, d_complex_t *b);

// iterative.c
int  gmres_solve(int n, const d_complex_t *a, const d_complex_t *alu, const int *piv,
	d_complex_t *b, double tol);

// solve.c
int  solve(peec_t *p, FILE *fp_log);

// output.c
void output_zin(const peec_t *p, FILE *fp_log);
int  output_csv(const peec_t *p, const char *fn);
void output_spara(const peec_t *p, FILE *fp_log);
int  output_touchstone(const peec_t *p);
int  output_dist(const peec_t *p, const char *fn);
int  output_pw(const peec_t *p, const char *fn, FILE *fp_log);

// farfield.c
int  output_far(const peec_t *p, const char *fn, FILE *fp_log);

#endif		// _PEEC_H_
