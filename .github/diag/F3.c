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

/*
放射ベクトル N(r^) = sum_m I_m len_m t^_m sinc(k len_m (r^・t^_m)/2) exp(+j k r^・r_cm)
(地板ありは鏡像セルの寄与も足す)。cur は区間電流 [0 .. nseg) の先頭
(port #1 励振なら segi、平面波なら segipw の周波数ごとの区画)。
*/
static void ff_nvec(const peec_t *p, const d_complex_t *cur, double kw,
	const double rh[3], d_complex_t nv[3])
{
	for (int c = 0; c < 3; c++) nv[c] = d_complex(0, 0);
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
			const d_complex_t w = d_rmul(amp, d_mul(cur[m], ph));
			for (int c = 0; c < 3; c++) {
				nv[c] = d_add(nv[c], d_rmul(tv[c], w));
			}
		}
	}
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

				d_complex_t nv[3];
				ff_nvec(p, &p->segi[DIDX(p, ifreq, 0, 0)], kw, rh, nv);

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


/*
平面波の散乱遠方界とレーダー断面積 (farpw.csv) — planewave + farfield のとき

誘起区間電流 (solve.c が捕捉する p->segipw) を放射源として、port #1 励振と
同じ式 (ff_nvec) で散乱界 r E_s を求める。入射振幅 |E0| = pwamp に対する
双方向 (bistatic) レーダー断面積は

  sigma(r^) = 4 pi |r E_s(r^)|^2 / |E0|^2                 [m^2]

で、far.csv と同じ格子・同じ求積で全散乱断面積

  sigma_s = (1/|E0|^2) ∫∫ |r E_s|^2 dOmega

を出す (地板ありは上半球)。地板が無いときは前方散乱振幅から光学定理で
消光断面積を出す。exp(+j omega t) 規約では (Bohren-Huffman の e^{-i omega t}
版の複素共役)

  sigma_e = -(4 pi / k) Im[ e^・r E_s(d) / E0 ],   d = -r^_inc (伝搬方向)

で、エネルギー保存により sigma_e = sigma_s + sigma_a (sigma_a は導体損と
負荷で吸収された電力 / 入射電力密度 |E0|^2/(2 eta0))。peec.log には
sigma_a = sigma_e - sigma_s を出す。

前方散乱振幅は離散系では回路が入射波から受け取る電力
P = (1/2) Re sum_m emf_m conj(I_m) と恒等的に結ばれる (起電力 emf_m と
放射ベクトルが同じ sinc・同じ位相規約を共有するため)。したがって
sigma_e = sigma_s + sigma_a が成り立つのは誘起電流に放射の反作用が入って
いるとき (retardation = 1) だけで、retardation = 0 では sigma_e は吸収分しか
含まない (ログに注意を出す)。地板があるときは入射が 2 波 (直接 + 反射) に
なって前方方向が 1 つに決まらないので sigma_e は出さない。
*/
int output_farpw(const peec_t *p, const char *fn, FILE *fp_log)
{
	if (!p->pw || (p->ffnth <= 0) || (p->ffnph <= 0)) return 0;
	if ((p->segipw == NULL) || (p->nseg <= 0)) return 0;

	const int nth = p->ffnth;
	const int nph = p->ffnph;
	const double thmax = p->gp ? (PI / 2) : PI;
	const double e2 = p->pwamp * p->pwamp;

	FILE *fp = fopen(fn, "w");
	if (fp == NULL) {
		printf("*** file %s open error.\n", fn);
		return 1;
	}

	const int ndir = (nth + 1) * nph;
	double *ss = (double *)malloc((size_t)ndir * sizeof(double));
	d_complex_t *et = (d_complex_t *)malloc((size_t)ndir * sizeof(d_complex_t));
	d_complex_t *ep = (d_complex_t *)malloc((size_t)ndir * sizeof(d_complex_t));
	if ((ss == NULL) || (et == NULL) || (ep == NULL)) {
		printf("%s\n", "*** memory allocation error (farfield, planewave)");
		free(ss); free(et); free(ep);
		fclose(fp);
		return 1;
	}

	if (!p->retardation) {
		fprintf(fp_log, "*** warning : farfield + planewave with retardation = 0 (no radiation reaction; sigma_ext then holds only the absorbed part)\n");
	}

	fprintf(fp, "frequency[Hz],theta[deg],phi[deg],rEtheta_real[V],rEtheta_imag[V],rEphi_real[V],rEphi_imag[V],RCS[m^2],RCS[dBsm]\n");

	for (int ifreq = 0; ifreq < p->nfreq; ifreq++) {
		const double f = freq_at(p, ifreq);
		const double kw = 2 * PI * f / C0;
		const double coef = 2 * PI * f * MU0 / (4 * PI);   // omega mu0 / 4 pi
		const d_complex_t *cur = &p->segipw[(size_t)ifreq * p->nseg];

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

				d_complex_t nv[3];
				ff_nvec(p, cur, kw, rh, nv);
				d_complex_t nt = d_complex(0, 0), np = d_complex(0, 0);
				for (int c = 0; c < 3; c++) {
					nt = d_add(nt, d_rmul(thh[c], nv[c]));
					np = d_add(np, d_rmul(phh[c], nv[c]));
				}
				const int id = (it * nph) + ip;
				et[id] = d_complex(coef * nt.i, -coef * nt.r);
				ep[id] = d_complex(coef * np.i, -coef * np.r);
				ss[id] = 4 * PI * (d_norm(et[id]) + d_norm(ep[id])) / e2;
			}
		}

		// 全散乱断面積 (far.csv の放射電力と同じ求積)
		double sigs = 0;
		for (int it = 0; it <= nth; it++) {
			const double wt = ((it == 0) || (it == nth)) ? 0.5 : 1;
			const double st = sin(thmax * it / nth);
			double row = 0;
			for (int ip = 0; ip < nph; ip++) {
				row += ss[(it * nph) + ip];
			}
			sigs += wt * st * row;
		}
		sigs *= (thmax / nth) * (2 * PI / nph) / (4 * PI);

		for (int it = 0; it <= nth; it++) {
			for (int ip = 0; ip < nph; ip++) {
				const int id = (it * nph) + ip;
				fprintf(fp, "%.9e,%.4f,%.4f,%.9e,%.9e,%.9e,%.9e,%.9e,%.4f\n",
					f, 180 / PI * thmax * it / nth, 360.0 * ip / nph,
					et[id].r, et[id].i, ep[id].r, ep[id].i, ss[id],
					(ss[id] > 1e-300) ? (10 * log10(ss[id])) : -3000);
			}
		}

		if (p->gp) {
			fprintf(fp_log, "farfield (planewave) : f = %.5e Hz, sigma_scat = %.9e m^2 (upper half space; sigma_ext not defined with ground plane)\n",
				f, sigs);
			continue;
		}

		// 光学定理 : 伝搬方向 d = -r^_inc の前方散乱振幅
		pw_wave_t w[2];
		pw_waves(p, kw, w);
		const double dv[3] = {-w[0].rh[0], -w[0].rh[1], -w[0].rh[2]};
		d_complex_t nv[3];
		ff_nvec(p, cur, kw, dv, nv);
		d_complex_t en = d_complex(0, 0);              // e^・N(d)
		for (int c = 0; c < 3; c++) {
			en = d_add(en, d_rmul(w[0].ev[c], nv[c]));
		}
		const d_complex_t efw = d_complex(coef * en.i, -coef * en.r);   // e^・r E_s(d)
		const double sige = -(4 * PI / kw) * d_div(efw, w[0].e0).i;
		fprintf(fp_log, "farfield (planewave) : f = %.5e Hz, sigma_scat = %.9e m^2, sigma_ext = %.9e m^2 (optical theorem), sigma_abs = %.9e m^2\n",
			f, sigs, sige, sige - sigs);
	}
	fflush(fp_log);

	free(ss);
	free(et);
	free(ep);
	fclose(fp);

	return 0;
}
