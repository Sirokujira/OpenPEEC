#!/usr/bin/env python3
"""OpenPEEC の CSV 出力を 1 個の HDF5 にまとめる (画面表示・後処理用)。

OpenPEEC 本体は外部ライブラリに依存しない方針なので (AGENTS.md の規則 8)、
HDF5 化はこの変換スクリプトが受け持つ。**ソルバーのビルドには不要**で、
このスクリプトを使うときだけ numpy と h5py が要る:

    pip install numpy h5py
    peec input.peec                       # CSV 一式が出る
    python tools/peec2h5.py . -o peec.h5  # まとめて HDF5 へ

入力 (あるものだけ読む):
    zin.csv    各ポートの Zin(f)
    tran.csv   過渡波形 (transient)          … 時系列
    pw.csv     平面波入射の応答 (planewave)
    far.csv    遠方界パターン (farfield)
    dist.csv   電流・電荷分布 (distribution)
    peec.sNp   S パラメータ (Touchstone 1.1)
    peec.log   タイトルなどのメタ情報

出力の構造 (すべて単位を attrs に持たせる):
    /                     title, source
    /frequency            [nf]                     Hz
    /zin                  [nport, nf] complex      ohm
    /sparam               [nf, nport, nport] complex
    /transient/time       [nt]                     s
    /transient/excitation [nt]                     規格化パルス (ピーク 1)
    /transient/s          [nport, nport, nt]       Sij の時間応答
    /transient/v          [nport, nt]              誘起端子電圧 V
    /planewave/voc        [nport, nf] complex      V
    /planewave/leff       [nport, nf]              m
    /planewave/pav        [nport, nf]              W
    /farfield/theta       [ntheta]                 deg
    /farfield/phi         [nphi]                   deg
    /farfield/e_theta     [nf, ntheta, nphi] complex   rE  V
    /farfield/e_phi       [nf, ntheta, nphi] complex   rE  V
    /farfield/d           [nf, ntheta, nphi]       dBi
    /farfield/g           [nf, ntheta, nphi]       dBi
    /distribution/<type>/{position,value,frequency}
        type = I (ポート励振の区間電流) / Ipw (平面波による誘起電流) /
               Q (容量セルの電荷)
"""

import argparse
import csv
import os
import sys

try:
    import numpy as np
except ImportError:                                     # pragma: no cover
    sys.exit("numpy が要る: pip install numpy h5py")
try:
    import h5py
except ImportError:                                     # pragma: no cover
    sys.exit("h5py が要る: pip install numpy h5py")


def read_csv(path):
    """ヘッダ付き CSV を dict のリストで返す。無ければ None。"""
    if not os.path.isfile(path):
        return None
    with open(path, newline="", encoding="utf-8") as fp:
        return list(csv.DictReader(fp))


def uniq(values):
    """出現順を保った一意化 (格子軸の復元用)。"""
    seen, out = set(), []
    for v in values:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def add(group, name, data, **attrs):
    dset = group.create_dataset(name, data=data, compression="gzip",
                                shuffle=True) if np.size(data) > 16 else \
        group.create_dataset(name, data=data)
    for key, val in attrs.items():
        dset.attrs[key] = val
    return dset


def conv_zin(h5, rows):
    """zin.csv -> /frequency, /zin"""
    ports = uniq(int(r["port"]) for r in rows)
    freqs = uniq(float(r["frequency[Hz]"]) for r in rows)
    z = np.zeros((len(ports), len(freqs)), dtype=np.complex128)
    for r in rows:
        i = ports.index(int(r["port"]))
        k = freqs.index(float(r["frequency[Hz]"]))
        z[i, k] = complex(float(r["Rin[ohm]"]), float(r["Xin[ohm]"]))
    add(h5, "frequency", np.array(freqs), units="Hz")
    add(h5, "zin", z, units="ohm", ports=np.array(ports, dtype="i4"))
    return len(ports), len(freqs)


def conv_touchstone(h5, path, nport, nfreq):
    """peec.sNp -> /sparam。2 ポートだけ列順が転置なのが Touchstone の仕様。"""
    if (nport is None) or (not os.path.isfile(path)):
        return False
    vals = []
    with open(path, encoding="utf-8") as fp:
        for line in fp:
            line = line.strip()
            if (not line) or line[0] in "!#":
                continue
            vals.extend(float(t) for t in line.split())

    stride = 1 + (2 * nport * nport)
    if (len(vals) % stride) != 0:
        return False
    nf = len(vals) // stride
    s = np.zeros((nf, nport, nport), dtype=np.complex128)
    for k in range(nf):
        row = vals[(k * stride) + 1:(k + 1) * stride]
        pairs = [complex(row[2 * m], row[(2 * m) + 1]) for m in range(nport * nport)]
        if nport == 2:
            # S11 S21 S12 S22 の順で格納されている
            s[k] = np.array([[pairs[0], pairs[2]], [pairs[1], pairs[3]]])
        else:
            s[k] = np.array(pairs).reshape(nport, nport)
    add(h5, "sparam", s, units="ratio",
        note="power-wave definition; see peec.log for per-port z0")
    return True


def conv_tran(h5, rows):
    """tran.csv -> /transient (時系列)"""
    grp = h5.create_group("transient")
    times = uniq(float(r["time[s]"]) for r in rows)
    tindex = {t: n for n, t in enumerate(times)}
    add(grp, "time", np.array(times), units="s")

    exc = np.zeros(len(times))
    sij, vpi = {}, {}
    for r in rows:
        n = tindex[float(r["time[s]"])]
        i, j, val = int(r["i"]), int(r["j"]), float(r["value"])
        if r["type"] == "X":
            exc[n] = val
        elif r["type"] == "S":
            sij.setdefault((i, j), np.zeros(len(times)))[n] = val
        elif r["type"] == "V":
            vpi.setdefault(i, np.zeros(len(times)))[n] = val

    add(grp, "excitation", exc, units="1",
        note="gaussian pulse, peak 1; responses are for this excitation")
    if sij:
        np_ = max(max(k) for k in sij)
        arr = np.zeros((np_, np_, len(times)))
        for (i, j), w in sij.items():
            arr[i - 1, j - 1] = w
        add(grp, "s", arr, units="1", note="time response of Sij")
    if vpi:
        arr = np.zeros((max(vpi), len(times)))
        for i, w in vpi.items():
            arr[i - 1] = w
        add(grp, "v", arr, units="V", note="plane-wave induced terminal voltage")
    return len(times)


def conv_pw(h5, rows):
    """pw.csv -> /planewave"""
    grp = h5.create_group("planewave")
    ports = uniq(int(r["port"]) for r in rows)
    freqs = uniq(float(r["frequency[Hz]"]) for r in rows)
    shape = (len(ports), len(freqs))
    voc = np.zeros(shape, dtype=np.complex128)
    leff = np.zeros(shape)
    pav = np.zeros(shape)
    for r in rows:
        i = ports.index(int(r["port"]))
        k = freqs.index(float(r["frequency[Hz]"]))
        voc[i, k] = complex(float(r["Voc_real[V]"]), float(r["Voc_imag[V]"]))
        leff[i, k] = float(r["leff[m]"])
        pav[i, k] = float(r["Pav[W]"])
    add(grp, "voc", voc, units="V",
        note="terminal voltage; equals open-circuit voltage when the port is unloaded")
    add(grp, "leff", leff, units="m")
    add(grp, "pav", pav, units="W")


def conv_far(h5, rows):
    """far.csv -> /farfield"""
    grp = h5.create_group("farfield")
    freqs = uniq(float(r["frequency[Hz]"]) for r in rows)
    thetas = uniq(float(r["theta[deg]"]) for r in rows)
    phis = uniq(float(r["phi[deg]"]) for r in rows)
    shape = (len(freqs), len(thetas), len(phis))
    et = np.zeros(shape, dtype=np.complex128)
    ep = np.zeros(shape, dtype=np.complex128)
    dd = np.zeros(shape)
    gg = np.zeros(shape)
    for r in rows:
        k = freqs.index(float(r["frequency[Hz]"]))
        it = thetas.index(float(r["theta[deg]"]))
        ip = phis.index(float(r["phi[deg]"]))
        et[k, it, ip] = complex(float(r["rEtheta_real[V]"]), float(r["rEtheta_imag[V]"]))
        ep[k, it, ip] = complex(float(r["rEphi_real[V]"]), float(r["rEphi_imag[V]"]))
        dd[k, it, ip] = float(r["D[dBi]"])
        gg[k, it, ip] = float(r["G[dBi]"])
    add(grp, "theta", np.array(thetas), units="deg")
    add(grp, "phi", np.array(phis), units="deg")
    add(grp, "e_theta", et, units="V", note="r*E with exp(-jkr)/r removed")
    add(grp, "e_phi", ep, units="V", note="r*E with exp(-jkr)/r removed")
    add(grp, "d", dd, units="dBi")
    add(grp, "g", gg, units="dBi")


def conv_dist(h5, rows):
    """dist.csv -> /distribution/<type>"""
    grp = h5.create_group("distribution")
    units = {"I": "A", "Ipw": "A", "Q": "C"}
    kinds = uniq(r["type"] for r in rows)
    for kind in kinds:
        sel = [r for r in rows if r["type"] == kind]
        freqs = uniq(float(r["frequency[Hz]"]) for r in sel)
        idx = uniq(int(r["index"]) for r in sel)
        val = np.zeros((len(freqs), len(idx)), dtype=np.complex128)
        pos = np.zeros((len(idx), 3))
        for r in sel:
            k = freqs.index(float(r["frequency[Hz]"]))
            m = idx.index(int(r["index"]))
            val[k, m] = complex(float(r["real"]), float(r["imag"]))
            pos[m] = (float(r["x[m]"]), float(r["y[m]"]), float(r["z[m]"]))
        sub = grp.create_group(kind)
        add(sub, "frequency", np.array(freqs), units="Hz")
        add(sub, "position", pos, units="m")
        add(sub, "value", val, units=units.get(kind, ""))


def read_title(path):
    if not os.path.isfile(path):
        return ""
    with open(path, encoding="utf-8", errors="replace") as fp:
        for line in fp:
            if line.startswith("title = "):
                return line[len("title = "):].strip()
    return ""


def main():
    ap = argparse.ArgumentParser(
        description="OpenPEEC の CSV 出力を HDF5 にまとめる")
    ap.add_argument("directory", nargs="?", default=".",
                    help="peec を実行したディレクトリ (既定 : カレント)")
    ap.add_argument("-o", "--output", default="peec.h5", help="出力 HDF5 (既定 : peec.h5)")
    args = ap.parse_args()

    dd = args.directory
    path = lambda name: os.path.join(dd, name)

    zin = read_csv(path("zin.csv"))
    if zin is None:
        sys.exit(f"{path('zin.csv')} が無い (peec を実行したディレクトリを指定する)")

    written = []
    with h5py.File(args.output, "w") as h5:
        h5.attrs["title"] = read_title(path("peec.log"))
        h5.attrs["source"] = "OpenPEEC"

        nport, nfreq = conv_zin(h5, zin)
        written.append("zin")

        if conv_touchstone(h5, path(f"peec.s{nport}p"), nport, nfreq):
            written.append("sparam")

        rows = read_csv(path("tran.csv"))
        if rows:
            nt = conv_tran(h5, rows)
            written.append(f"transient ({nt} samples)")

        rows = read_csv(path("pw.csv"))
        if rows:
            conv_pw(h5, rows)
            written.append("planewave")

        rows = read_csv(path("far.csv"))
        if rows:
            conv_far(h5, rows)
            written.append("farfield")

        rows = read_csv(path("dist.csv"))
        if rows:
            conv_dist(h5, rows)
            written.append("distribution")

    print(f"{args.output} : {', '.join(written)}")


if __name__ == "__main__":
    main()
