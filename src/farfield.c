/*
farfield.c

遠方界の後処理 (farfield = 分割数theta 分割数phi)

port #1 を 1 A で励振したときの区間電流 I_m (solve.c が捕捉する p->segi)
から放射ベクトル

  N(r^) = sum_m I_m len_m t^_m sinc(k len_m (r^・t^_m)/2) exp(+j k r^・r_cm)

を求める (区間内は一定電流なので、線形位相の区間積分は厳密に sinc になる)。
遠方界は横成分だけが残り

  r E = -j (omega mu0 / 4 pi) [N - (N・r^) r^]   … exp(-jkr)/r を除いた値 [V]

で、far.csv には r E の theta / phi 成分と指向性 D・利得 G を出力する。

放射電力は球面 (地板ありは上半球) のグリッド数値積分

  Prad = ∫∫ (|Et|^2 + |Ep|^2) / (2 eta0) sin(theta) dtheta dphi

(theta は台形則、phi は周期矩形則)。指向性 D = 4 pi U / Prad、
利得 G = 4 pi U / Pin (Pin = Re(Zin_port1) / 2 : 1 A 励振の入力電力)、
放射効率 eff = Prad / Pin。G は導体損を含み、D は含まない。

地板 (groundplane) があるときは鏡像電流を加算して上半球のみを評価する。
鏡像は幾何鏡像セル (方向 t' = (tx, ty, -tz)) に -I を流したもの
(= 水平電流反転・垂直電流保存。partial.c の鏡像規約と同一)。

retardation = 0 の準静的電流でもパターンは計算できるが、電流に放射の
反作用が含まれないため効率は物理的な意味を持たない (ログに注意を出す)。
*/

#include "peec.h"

// sin(x)/x (x -> 0 で 1)
static double sinc(double x)
{
	return (fabs(x) < 1e-8) ? 1.0 : (sin(x) / x);
}

int output_far(const peec_t *p, const char *fn, FILE *fp_log)
{
	if ((p->ffnth <= 0) || (p->ffnph <= 0)) return 0;
	if ((p->segi == NULL) || (p->nseg <= 0) || (p->nport <= 0)) return 0;

	const int nth = p->ffnth;
	const int nph = p->ffnph;
	const double thmax = p->gp ? (PI / 2) : PI;
	const double eta0 = MU0 * C0;

	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fn);
		return 1;
	}

	if (!p->retardation) {
		fprintf(fp_log, "*** warning : farfield with retardation = 0 (quasi-static currents do not radiate; efficiency is not meaningful)\n");
	}

	const int ndir = (nth + 1) * nph;
	d_complex_t *et = (d_complex_t *)malloc((size_t)ndir * sizeof(d_complex_t));
	d_complex_t *ep = (d_complex_t *)malloc((size_t)ndir * sizeof(d_complex_t));
	double *uu = (double *)malloc((size_t)ndir * sizeof(double));
	if ((et == NULL) || (ep == NULL) || (uu == NULL)) {
		printf("%s\n", "*** memory allocation error (farfield)");
		free(et); free(ep); free(uu);
		fclose(fp);
		return 1;
	}

	fprintf(fp, "frequency[Hz],theta[deg],phi[deg],rEtheta_real[V],rEtheta_imag[V],rEphi_real[V],rEphi_imag[V],D[dBi],G[dBi]\n");

	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		const double f = freq_at(p, ifreq);
		const double kw = 2 * PI * f / C0;
		const double coef = 2 * PI * f * MU0 / (4 * PI);   // omega mu0 / 4 pi

		// 方向ごとの遠方界と放射強度 U
		for (int it = 0; it <= nth; it++) {
			const double th = thmax * it / nth;
			const double st = sin(th);
			const double ct = cos(th);
			for (int ip = 0; ip < nph; ip++) {
				const double phv = 2 * PI * ip / nph;
				const double cp = cos(phv);
				const double sp = sin(phv);
				const double rh[3] = {st * cp, st * sp, ct};
				const double thh[3] = {ct * cp, ct * sp, -st};
				const double phh[3] = {-sp, cp, 0};

				d_complex_t nv[3] = {{0, 0}, {0, 0}, {0, 0}};
				const int npass = p->gp ? 2 : 1;
				for (int pass = 0; pass < npass; pass++) {
					for (int m = 0; m < p->nseg; m++) {
						const seg_t *s = &p->seg[m];
						double tv[3], rc[3];
						for (int c = 0; c < 3; c++) {
							tv[c] = (s->x2[c] - s->x1[c]) / s->len;
							rc[c] = 0.5 * (s->x1[c] + s->x2[c]);
						}
						double sgn = 1;
						if (pass == 1) {
							// 鏡像 : 幾何鏡像 + 電流 -I (水平反転・垂直保存)
							tv[2] = -tv[2];
							rc[2] = (2 * p->gpz) - rc[2];
							sgn = -1;
						}
						const double rt = (rh[0] * tv[0]) + (rh[1] * tv[1]) + (rh[2] * tv[2]);
						const double rr = (rh[0] * rc[0]) + (rh[1] * rc[1]) + (rh[2] * rc[2]);
						const double amp = sgn * s->len * sinc(0.5 * kw * s->len * rt);
						const d_complex_t ph = d_complex(cos(kw * rr), sin(kw * rr));
						const d_complex_t w = d_rmul(amp,
							d_mul(p->segi[DIDX(p, ifreq, 0, m)], ph));
						for (int c = 0; c < 3; c++) {
							nv[c] = d_add(nv[c], d_rmul(tv[c], w));
						}
					}
				}

				// 横成分 : Etheta = -j coef (N・theta^)、Ephi = -j coef (N・phi^)
				d_complex_t nt = d_complex(0, 0), np = d_complex(0, 0);
				for (int c = 0; c < 3; c++) {
					nt = d_add(nt, d_rmul(thh[c], nv[c]));
					np = d_add(np, d_rmul(phh[c], nv[c]));
				}
				const int id = (it * nph) + ip;
				et[id] = d_complex(coef * nt.i, -coef * nt.r);
				ep[id] = d_complex(coef * np.i, -coef * np.r);
				uu[id] = (d_norm(et[id]) + d_norm(ep[id])) / (2 * eta0);
			}
		}

		// 放射電力 (theta : 台形則、phi : 周期矩形則)
		double prad = 0;
		for (int it = 0; it <= nth; it++) {
			const double wt = ((it == 0) || (it == nth)) ? 0.5 : 1;
			const double st = sin(thmax * it / nth);
			double row = 0;
			for (int ip = 0; ip < nph; ip++) {
				row += uu[(it * nph) + ip];
			}
			prad += wt * st * row;
		}
		prad *= (thmax / nth) * (2 * PI / nph);

		// 入力電力 (port #1 に 1 A) と最大指向性
		const double pin = 0.5 * p->zin[ifreq].r;
		int imax = 0;
		for (int id = 1; id < ndir; id++) {
			if (uu[id] > uu[imax]) imax = id;
		}

		for (int it = 0; it <= nth; it++) {
			for (int ip = 0; ip < nph; ip++) {
				const int id = (it * nph) + ip;
				const double d = (prad > 0) ? (4 * PI * uu[id] / prad) : 0;
				const double g = (pin > 0) ? (4 * PI * uu[id] / pin) : 0;
				fprintf(fp, "%.9e,%.4f,%.4f,%.9e,%.9e,%.9e,%.9e,%.4f,%.4f\n",
					f, 180 / PI * thmax * it / nth, 360.0 * ip / nph,
					et[id].r, et[id].i, ep[id].r, ep[id].i,
					(d > 1e-30) ? (10 * log10(d)) : -600,
					(g > 1e-30) ? (10 * log10(g)) : -600);
			}
		}

		const double dmax = (prad > 0) ? (4 * PI * uu[imax] / prad) : 0;
		fprintf(fp_log, "farfield : f = %.5e Hz, Prad = %.5e W, Pin = %.5e W, eff = %.5f, max D = %.3f dBi at (theta, phi) = (%.1f, %.1f) deg\n",
			f, prad, pin, (pin > 0) ? (prad / pin) : 0,
			(dmax > 1e-30) ? (10 * log10(dmax)) : -600,
			180 / PI * thmax * (imax / nph) / nth, 360.0 * (imax % nph) / nph);
	}
	fflush(fp_log);

	free(et);
	free(ep);
	free(uu);
	fclose(fp);

	return 0;
}
