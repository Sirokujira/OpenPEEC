#!/usr/bin/env python3
"""
basis_study.py — 細線の電流基底の次数の検討 (ソルバー本体とは独立な試作、純 Python)

OpenPEEC の細線は「区間ごとに一定の電流 (パルス) + 節点セルごとに一定の電荷
(パルス)」の互い違い格子 (staggered) で、節点電位はセル平均 (平均電位法) である。
ここでは、その電流にだけ 1 次の補正を足した高次基底を Galerkin で組み、
収束が良くなるかを直線ダイポールで調べる。ソルバー本体には採用していない
(結論は下記)。README の「既知の制限」から参照している。

高次基底 (basis = linear) の定式化
  電流 J = Σ_k I_k p_k + Σ_n c_n χ_n,  c_n = jω q_n
    p_k : 区間 k のパルス (従来の PEEC)
    χ_n : 節点 n の電荷セル (両隣の半区間) 上の傾斜関数。半区間 h で
          (L_h - s)/L_n (s は節点からの距離、向きは節点から外向き、L_n = Σ L_h)
  KCL (Σ I_out + jω q = 注入) と合わせると電流は連続な区分線形
  (節点と区間中点で折れる) になり、電荷はセル内一様のまま連続の式を満たす。
  試験関数も {p_k, χ_n} (Galerkin)。χ_n で試験すると
    <χ_n, -∇φ> = φ(節点) - φ̄_n (セル平均)
  となるので、節点電位を「節点での点電位」とすれば未知数 (φ, I, q) と式の数が
  揃い、回路素子は点電位に繋がる。従来の PEEC の式に足されるのは
    KVL_k  : - jω Σ_n K_kn c_n
    電荷行 n : φ_n = Σ_m P_nm q_m + jω (Σ_l K_ln I_l + Σ_m X_nm c_m)
  K_kn = <p_k, A[χ_n]>、X_nm = <χ_n, A[χ_m]> (縮約カーネル + 遅延) で、
  ω -> 0 では従来の PEEC に一致する。

結果 (半波長ダイポール、l = 1 m、a = 1 mm、ギャップ 1 mm、f = 143.76 MHz)
  片腕の分割数    4         8         16        32       (Rin [ohm])
  パルス (従来)   80.243    74.216    72.716    72.368    (64 : 72.338)
  線形補正        70.058    71.434    71.903    72.114
  パルスは誤差 +7.9 -> +1.9 -> +0.37 -> +0.02 で 2 次以上、線形補正は
  -2.3 -> -0.92 -> -0.45 -> -0.24 で 1 次 (粗い格子だけ有利)。Xin はデルタギャップ
  給電のギャップ容量が格子を細かくすると対数的に発散する (既知の性質) ので
  比較に使わない。
  給電の無い散乱 (1 m の線に broadside の平面波、誘起双極子モーメント
  M = ∫I dz) でも、N = 4/8/16/32 のパルスの列から外挿した極限
  M∞ ≈ (5.966 - 0.515j) x 1e-3 A m に対して、N = 32 の誤差はパルス 1.5e-5 /
  線形補正 4.6e-5 (相対 2.5e-3 / 7.7e-3)。収束もパルス (各段 1/6) の方が速く、
  線形補正は虚部が各段 1/3 程度。

結論 : 電流だけを高次にすると、電荷が区分一定のままなので互い違い格子の
  対称性 (パルス電流 + パルス電荷の誤差の打ち消し) が崩れて 1 次に落ちる。
  次数を上げるには電荷も高次 (区分線形) にする必要があり、その場合は節点電位を
  セル平均で定義する PEEC の回路解釈 (集中定数・ポートが繋がる節点) が成り立た
  ない。よってソルバーは従来の基底のままとし、精度は分割数で確保する
  (半波長ダイポールで片腕 16 分割 → Rin 誤差 0.5%)。

使い方 (numpy 不要。32 分割で数分かかる)
  python3 tools/basis_study.py 4 8 16          # ダイポールの Zin
  python3 tools/basis_study.py --scat 4 8 16   # 散乱 (誘起双極子モーメント)
パルス基底の値は OpenPEEC 本体 (dipole_ff.peec の分割数を変えたもの) と
1e-5 で一致する (遅延補正の求積の差だけ)。
"""

import cmath
import math
import sys

C0 = 2.99792458e8
MU0 = 4e-7 * math.pi
EPS0 = 1 / (MU0 * C0 * C0)
A = 1e-3          # 線半径 (縮約カーネルの正則化)
GAP = 0.0005      # 給電ギャップの半幅
SIGMA = 5.8e7


def gauss_legendre(n):
    """[0, 1] 上の Gauss-Legendre 点と重み"""
    xs, ws = [], []
    for i in range(1, n + 1):
        z = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p1, p2 = 1.0, 0.0
            for j in range(1, n + 1):
                p3 = p2
                p2 = p1
                p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
            pp = n * (z * p1 - p2) / (z * z - 1)
            z1 = z
            z = z1 - p1 / pp
            if abs(z - z1) < 1e-15:
                break
        xs.append(0.5 * (1 - z))
        ws.append(1.0 / ((1 - z * z) * pp * pp))
    return xs, ws


GX, GW = gauss_legendre(8)   # 遅延補正 (なめらかな被積分関数)
TX, TW = gauss_legendre(6)   # 傾斜関数をパルスの重ね合わせで表す t 積分


def prim2(u):
    """1/sqrt(u^2 + a^2) の 2 階原始関数"""
    return u * math.asinh(u / A) - math.sqrt(u * u + A * A)


def pulse_pair(a1, b1, a2, b2):
    """∫_{a1}^{b1} ∫_{a2}^{b2} ds ds' / sqrt((s-s')^2 + a^2) (同一直線上、閉形式)"""
    return prim2(b1 - a2) - prim2(a1 - a2) - prim2(b1 - b2) + prim2(a1 - b2)


def retard(seg1, w1, seg2, w2, k):
    """∫∫ w1 w2 (exp(-jkR) - 1)/R (縮約距離、なめらかなので Gauss で十分)"""
    if k == 0:
        return 0
    a1, b1 = seg1
    a2, b2 = seg2
    tot = 0
    for xi, wi in zip(GX, GW):
        s = a1 + (b1 - a1) * xi
        for xj, wj in zip(GX, GW):
            t = a2 + (b2 - a2) * xj
            r = math.sqrt((s - t) ** 2 + A * A)
            tot += wi * wj * w1(s) * w2(t) * (cmath.exp(-1j * k * r) - 1) / r
    return tot * (b1 - a1) * (b2 - a2)


def one(_s):
    return 1.0


def lu_solve(m, b):
    n = len(b)
    m = [row[:] for row in m]
    b = b[:]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(m[r][c]))
        m[c], m[p] = m[p], m[c]
        b[c], b[p] = b[p], b[c]
        for r in range(c + 1, n):
            fct = m[r][c] / m[c][c]
            if fct != 0:
                mr, mc = m[r], m[c]
                for j in range(c, n):
                    mr[j] -= fct * mc[j]
                b[r] -= fct * b[c]
    x = [0j] * n
    for r in range(n - 1, -1, -1):
        x[r] = (b[r] - sum(m[r][j] * x[j] for j in range(r + 1, n))) / m[r][r]
    return x


def solve(nper, f, linear, scat=False):
    """nper : 片腕 (scat なら全体) の分割数。戻り値 : Zin、または散乱の M = ∫I dz"""
    k = 2 * math.pi * f / C0
    w = 2 * math.pi * f
    if scat:
        nodes = [-0.5 + i / nper for i in range(nper + 1)]
        segs = [(i, i + 1, nodes[i], nodes[i + 1]) for i in range(nper)]
    else:
        z1 = [-0.5 + (0.5 - GAP) * i / nper for i in range(nper + 1)]
        z2 = [GAP + (0.5 - GAP) * i / nper for i in range(nper + 1)]
        nodes = z1 + z2
        segs = [(i, i + 1, z1[i], z1[i + 1]) for i in range(nper)]
        segs += [(nper + 1 + i, nper + 2 + i, z2[i], z2[i + 1]) for i in range(nper)]
    nn, ns = len(nodes), len(segs)
    # 電荷セル : 節点ごとの半区間 (始点, 終点, 節点の z, 外向きの符号)
    cells = [[] for _ in range(nn)]
    for (n1, n2, za, zb) in segs:
        zm = 0.5 * (za + zb)
        cells[n1].append((za, zm, za, +1))
        cells[n2].append((zm, zb, zb, -1))
    lc = [sum(h[1] - h[0] for h in c) for c in cells]

    def ramp(h):
        lh = h[1] - h[0]
        return lambda s: lh - abs(s - h[2])

    def sub(h, t):
        # 傾斜 L_h - |s - 節点| = ∫_0^{L_h} [節点から長さ t のパルス] dt
        return (h[2], h[2] + t) if h[3] > 0 else (h[2] - t, h[2])

    def pulse_ramp(s1, h):
        lh = h[1] - h[0]
        return lh * sum(wt * pulse_pair(s1[0], s1[1], *sub(h, lh * x)) for x, wt in zip(TX, TW))

    def ramp_ramp(h, g):
        lh, lg = h[1] - h[0], g[1] - g[0]
        tot = 0
        for x, wx in zip(TX, TW):
            u = sub(h, lh * x)
            for y, wy in zip(TX, TW):
                tot += wx * wy * pulse_pair(u[0], u[1], *sub(g, lg * y))
        return tot * lh * lg

    cm = MU0 / (4 * math.pi)
    cp = 1 / (4 * math.pi * EPS0)
    lp = [[cm * (pulse_pair(si[2], si[3], sj[2], sj[3])
                 + retard((si[2], si[3]), one, (sj[2], sj[3]), one, k))
           for sj in segs] for si in segs]
    pm = [[cp * sum(pulse_pair(h[0], h[1], g[0], g[1]) + retard((h[0], h[1]), one, (g[0], g[1]), one, k)
                    for h in cells[n] for g in cells[m]) / (lc[n] * lc[m])
           for m in range(nn)] for n in range(nn)]
    kk = [[0] * nn for _ in range(ns)]
    xx = [[0] * nn for _ in range(nn)]
    if linear:
        for i, si in enumerate(segs):
            for n in range(nn):
                kk[i][n] = cm * sum(h[3] * (pulse_ramp((si[2], si[3]), h)
                                            + retard((si[2], si[3]), one, (h[0], h[1]), ramp(h), k))
                                    for h in cells[n]) / lc[n]
        for n in range(nn):
            for m in range(nn):
                xx[n][m] = cm * sum(h[3] * g[3] * (ramp_ramp(h, g)
                                                   + retard((h[0], h[1]), ramp(h), (g[0], g[1]), ramp(g), k))
                                    for h in cells[n] for g in cells[m]) / (lc[n] * lc[m])
    # 未知数 [φ (nn), I (ns), q (nn)]
    nt = 2 * nn + ns
    oi, oq = nn, nn + ns
    mat = [[0j] * nt for _ in range(nt)]
    for i, (n1, n2, za, zb) in enumerate(segs):
        r = oi + i
        mat[n1][r] += 1
        mat[n2][r] -= 1
        mat[r][n1] += 1
        mat[r][n2] -= 1
        mat[r][r] -= (zb - za) / (SIGMA * math.pi * A * A)
        for j in range(ns):
            mat[r][oi + j] -= 1j * w * lp[i][j]
        for m in range(nn):
            mat[r][oq + m] += w * w * kk[i][m]
    for n in range(nn):
        mat[n][oq + n] += 1j * w
        r = oq + n
        mat[r][n] -= 1
        for m in range(nn):
            mat[r][oq + m] += pm[n][m] - w * w * xx[n][m]
        for j in range(ns):
            mat[r][oi + j] += 1j * w * kk[j][n]
    rhs = [0j] * nt
    if scat:
        # broadside の E_z = 1 V/m : KVL に -<p_k, E>、電荷行に +<χ_n, E>
        mom_c = [sum(h[3] * (h[1] - h[0]) ** 2 / 2 for h in cells[n]) / lc[n] for n in range(nn)]
        for i, (n1, n2, za, zb) in enumerate(segs):
            rhs[oi + i] = -(zb - za)
        if linear:
            for n in range(nn):
                rhs[oq + n] = mom_c[n]
        x = lu_solve(mat, rhs)
        mom = sum(x[oi + i] * (s[3] - s[2]) for i, s in enumerate(segs))
        if linear:
            mom += sum(1j * w * x[oq + n] * mom_c[n] for n in range(nn))
        return mom
    rhs[nper + 1] = 1     # 1 A を腕 2 の根元へ注入し、腕 1 の根元から戻す
    rhs[nper] = -1
    x = lu_solve(mat, rhs)
    return x[nper + 1] - x[nper]


def main(argv):
    scat = "--scat" in argv
    ns = [int(v) for v in argv if v != "--scat"] or [4, 8, 16]
    f = 1.4376e8
    for n in ns:
        z0 = solve(n, f, False, scat)
        z1 = solve(n, f, True, scat)
        label = "M" if scat else "Zin"
        print("n=%3d  pulse %s = %.7g %+.7gj   linear %s = %.7g %+.7gj"
              % (n, label, z0.real, z0.imag, label, z1.real, z1.imag), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
