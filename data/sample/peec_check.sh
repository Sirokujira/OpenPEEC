#!/bin/sh
# peec_check.sh — OpenPEEC 検証 (CI 用)
#
# data/sample/ の各ケースを実行し、zin.csv / peec.log を解析解と比較する。
# 期待値の導出は各 .peec ファイルのコメント参照。
#
# 使い方 : peec_check.sh <peec 実行ファイル(絶対パス)> [作業ディレクトリ]

set -e

PEEC="$1"
WORK="${2:-.}"
SRC="$(cd "$(dirname "$0")" && pwd)"
TWOPI=6.283185307179586

if [ -z "$PEEC" ]; then
	echo "Usage: peec_check.sh <peec> [workdir]" >&2
	exit 2
fi

mkdir -p "$WORK"
status=0
csv="$WORK/zin.csv"
log="$WORK/peec.log"

# run <input.peec> : 実行して normal end を確認
run() {
	(cd "$WORK" && "$PEEC" -n 2 "$1" > /dev/null)
	grep -q "normal end" "$log"
}

# chk <label> <actual> <expected> <tol> : 相対誤差判定
chk() {
	awk -v a="$2" -v e="$3" -v tol="$4" -v lb="$1" 'BEGIN {
		d = (a - e) / e; ad = (d < 0) ? -d : d;
		printf "%-24s actual=%.6g expected=%.6g -> %s %+.4f%%\n", lb, a, e, (ad <= tol) ? "OK" : "NG", d * 100;
		exit (ad <= tol) ? 0 : 1
	}' || status=1
}

# 1 行目のデータ行から L = Xin / (2 pi f) と Rin を取り出す
getL() { awk -F, -v tp="$TWOPI" 'NR==2{printf "%.9e", $4/(tp*$2)}' "$csv"; }
getR() { awk -F, 'NR==2{print $3}' "$csv"; }
# "PEEC: capacitive cells = N, total capacitance = C F" から C を取り出す
getC() { grep "total capacitance" "$log" | tail -1 | awk '{print $(NF-1)}'; }

# peec.log の "S<i><j>" 表の 1 行目から実部 / 虚部を取り出す
getS() { awk -v k="$1" -v c="$2" '$0 == k {found = 1; next} found && (++n == 2) {print $(c + 1); exit}' "$log"; }

echo "--- MNA (lumped elements)"
# (a) 直列 RLC : Rin(f0)=50, Xin(f0)~0, Xin(2f0)=47.43416
cp "$SRC/rlc_series.peec" "$WORK/"
run rlc_series.peec
chk "RLC Rin(f0)" "$(awk -F, 'NR==2{print $3}' "$csv")" 50.0 0.005
awk -F, 'NR==2 {
	x = $4; ax = (x < 0) ? -x : x;
	printf "%-24s actual=%.6g -> %s (|Xin| <= 0.05)\n", "RLC Xin(f0)", x, (ax <= 0.05) ? "OK" : "NG";
	exit (ax <= 0.05) ? 0 : 1
}' "$csv" || status=1
chk "RLC Xin(2f0)" "$(awk -F, 'NR==3{print $4}' "$csv")" 47.43416 0.005

echo "--- PEEC partial inductance"
# (b) 単線ワイヤ : Lp=1.320380e-6 H, Rin(DC)=5.48810e-3 ohm
cp "$SRC/wire_single.peec" "$WORK/"
run wire_single.peec
chk "wire L (ndiv=1)" "$(getL)" 1.320380e-6 0.001
chk "wire Rin (DC)" "$(getR)" 5.48810e-3 0.005
# 細線カーネルを自己項・相互項で統一しているので L は分割数に依存しない
sed 's/5.8e7 1$/5.8e7 8/' "$SRC/wire_single.peec" > "$WORK/wire_ndiv8.peec"
run wire_ndiv8.peec
chk "wire L (ndiv=8)" "$(getL)" 1.320380e-6 0.001

# (c) 正方ループ : L=7.24689e-7 H
cp "$SRC/loop_square.peec" "$WORK/"
run loop_square.peec
chk "loop L" "$(getL)" 7.24689e-7 0.02

echo "--- skin effect"
# (d) 低周波 : Rin=5.494087e-3 ohm, L_tot=1.370349e-6 H (内部 L -> mu0 l/(8 pi))
cp "$SRC/wire_skin.peec" "$WORK/"
run wire_skin.peec
chk "skin LF Rin" "$(getR)" 5.494087e-3 0.005
chk "skin LF L(ext+int)" "$(getL)" 1.370349e-6 0.005

# (e) 高周波 : Rin -> l/(2 pi a delta sigma) = 0.0415230 ohm (曲率補正 +0.33%)
cp "$SRC/wire_skin_hf.peec" "$WORK/"
run wire_skin_hf.peec
chk "skin HF Rin" "$(getR)" 0.0415230 0.01

echo "--- capacitive PEEC (potential coefficients)"
# (f) 単線の対無限遠容量 : 8.42673e-12 F
cp "$SRC/wire_cap.peec" "$WORK/"
run wire_cap.peec
chk "wire Ctotal" "$(getC)" 8.42673e-12 0.001

# (g) 平行 2 線の対無限遠容量 (相互電位係数の検証) : 1.0202867e-11 F
cp "$SRC/twowire_cap.peec" "$WORK/"
run twowire_cap.peec
chk "two-wire Ctotal" "$(getC)" 1.0202867e-11 0.001

# (h) 電気的に小さいループでは容量を入れても Zin が変わらないこと
awk '{print} /^title = /{print "capacitance = 1"}' "$SRC/loop_square.peec" > "$WORK/loop_cap.peec"
run loop_cap.peec
chk "loop L (with C)" "$(getL)" 7.24689e-7 0.02

# (i) フル PEEC (L+P+R) : ダイポールが容量性から誘導性へ移る直列共振を持つこと
cp "$SRC/dipole_full.peec" "$WORK/"
run dipole_full.peec
awk -F, 'NR>1 {
	if ($4 < 0) neg = 1;
	if (neg && ($4 > 0)) pos = 1;
} END {
	printf "%-24s %s (series resonance present)\n", "dipole Xin sign change", pos ? "OK" : "NG";
	exit pos ? 0 : 1
}' "$csv" || status=1

echo "--- retardation (full-wave PEEC)"
# (j) 半波長ダイポール : Xin の零交差で l/lambda = 0.478、Rin = 73.1 ohm (教科書値)
cp "$SRC/dipole_halfwave.peec" "$WORK/"
run dipole_halfwave.peec
res=$(awk -F, 'NR>1 {
	f = $2; r = $3; x = $4;
	if (pf && (px < 0) && (x >= 0) && !done) {
		u = -px / (x - px);              # Xin の零交差を線形補間
		printf "%.9e %.9e", (pf + u*(f-pf)) / 2.99792458e8, pr + u*(r-pr);
		done = 1;
	}
	pf = f; pr = r; px = x;
}' "$csv")
if [ -z "$res" ]; then
	echo "*** no resonance (Xin zero crossing) in dipole_halfwave sweep" >&2
	status=1
else
	chk "dipole l/lambda" "${res% *}" 0.478 0.05
	chk "dipole Rin at res." "${res#* }" 73.1 0.10
fi

# (k) 微小ダイポール : R_rad = 20 pi^2 (l/lambda)^2 = 0.123541 ohm
cp "$SRC/dipole_short.peec" "$WORK/"
run dipole_short.peec
chk "short dipole R_rad" "$(getR)" 0.123541 0.05

# 低周波では遅延ありでも準静的な結果に一致すること
awk '{print} /^title = /{print "retardation = 1"}' "$SRC/loop_square.peec" > "$WORK/loop_ret.peec"
run loop_ret.peec
chk "loop L (retarded)" "$(getL)" 7.24689e-7 0.02

echo "--- acceleration (sweep-reuse GMRES)"
# 掃引 LU 再利用 GMRES (acceleration = 1) が密 LU と同じ結果になること。
# 収束判定 1e-10 なので Zin の差は ~1e-9 以下になるはず (判定は 1e-8)。
accel_cmp() {
	paste -d, "$WORK/zin_ref.csv" "$csv" | awk -F, -v lb="$1" 'NR > 1 {
		dr = $3 - $10; di = $4 - $11;
		e = sqrt((dr * dr) + (di * di)) / (sqrt(($3 * $3) + ($4 * $4)) + 1);
		if (e > mx) mx = e; n++
	} END {
		ok = (n > 0) && (mx <= 1e-8);
		printf "%-24s rows=%d max|dZ|/(|Z|+1)=%.3e -> %s\n", lb, n, mx, ok ? "OK" : "NG";
		exit ok ? 0 : 1
	}' || status=1
}
# (v2) フル PEEC (L+P+R) の掃引
cp "$SRC/dipole_full.peec" "$WORK/"
run dipole_full.peec
cp "$csv" "$WORK/zin_ref.csv"
awk '{print} /^title = /{print "acceleration = 1"}' "$SRC/dipole_full.peec" > "$WORK/dipole_full_ac.peec"
run dipole_full_ac.peec
accel_cmp "accel full PEEC sweep"
# (w2) 遅延あり (行列が毎周波数変わる) でも一致すること
cp "$SRC/dipole_halfwave.peec" "$WORK/"
run dipole_halfwave.peec
cp "$csv" "$WORK/zin_ref.csv"
awk '{print} /^title = /{print "acceleration = 1"}' "$SRC/dipole_halfwave.peec" > "$WORK/dipole_hw_ac.peec"
run dipole_hw_ac.peec
accel_cmp "accel retarded sweep"
# (x2) スレッド数不変性 : GMRES の内積・Gram-Schmidt は直列、
#      行列ベクトル積は行ごとに独立なのでビット単位で一致する
(cd "$WORK" && "$PEEC" -n 1 dipole_full_ac.peec > /dev/null)
cp "$csv" "$WORK/zin_ac_n1.csv"
(cd "$WORK" && "$PEEC" -n 4 dipole_full_ac.peec > /dev/null)
if cmp -s "$WORK/zin_ac_n1.csv" "$csv"; then
	printf "%-24s -> OK (-n 1 と -n 4 が完全一致)\n" "accel thread invariance"
else
	printf "%-24s -> NG (-n 1 と -n 4 が不一致)\n" "accel thread invariance" >&2
	status=1
fi

echo "--- rectangular bar conductors"
# (l) 角線の自己インダクタンス (Grover) と DC 抵抗
cp "$SRC/bar_single.peec" "$WORK/"
run bar_single.peec
chk "bar L (Grover)" "$(getL)" 1.141093e-6 0.005
chk "bar Rin (DC)" "$(getR)" 1.724138e-3 0.005

# (m) 薄板の対無限遠容量 (薄板 <-> 半径 w/4 の丸線)
cp "$SRC/bar_strip_cap.peec" "$WORK/"
run bar_strip_cap.peec
chk "strip Ctotal" "$(getC)" 9.782209e-12 0.005

# (n) 角線の表皮効果 (高周波極限 R -> l/(sigma delta P))
cp "$SRC/bar_skin.peec" "$WORK/"
run bar_skin.peec
chk "bar skin HF Rin" "$(getR)" 0.118589 0.01

echo "--- surface (plate) conductors"
# (o) 正方形平板の対無限遠容量 : 8 分割と 16 分割から Richardson 補外
#     C_inf = 2 C_16 - C_8 を文献値 0.3667892 * 4 pi eps0 a = 4.08100e-11 F と比較
cp "$SRC/plate_cap.peec" "$WORK/"
run plate_cap.peec
c8=$(getC)
sed 's/5.8e7 8 8/5.8e7 16 16/' "$SRC/plate_cap.peec" > "$WORK/plate_cap16.peec"
run plate_cap16.peec
c16=$(getC)
chk "plate C (extrapolated)" "$(awk -v a="$c8" -v b="$c16" 'BEGIN{printf "%.9e", 2*b-a}')" 4.08100e-11 0.01

# (p) 帯のインダクタンス : 面積分 (plate) と GMD 近似 (bar) の相互検証
cp "$SRC/plate_strip.peec" "$WORK/"
run plate_strip.peec
chk "strip L (surface)" "$(getL)" 1.159664e-6 0.01

echo "--- volume conductor (plate ndivt)"
# (v) 厚板の帯 : 体積セル (Hoer-Love) の自己 L が Grover の角線式と一致し、
#     DC 抵抗が解析値かつ厚み層数に厳密に依存しないこと
cp "$SRC/plate_thick.peec" "$WORK/"
run plate_thick.peec
chk "thick plate L (Grover)" "$(getL)" 1.1237356e-6 0.005
chk "thick plate Rin (DC)" "$(getR)" 8.6206897e-4 0.001
rnt4=$(getR)
for nt in 1 2; do
	sed "s/16 1 4\$/16 1 $nt/" "$SRC/plate_thick.peec" > "$WORK/plate_thick_nt$nt.peec"
	run plate_thick_nt$nt.peec
	chk "thick plate R (ndivt=$nt)" "$(getR)" "$rnt4" 1e-6
done

# (w) 厚み方向の表皮効果 : 対称サンドイッチの R_ac/R_dc が層数を増やすと
#     1D スラブ解 Re[(1+j)X coth((1+j)X)] (X = t/2delta = 2) = 1.897806 に収束
cp "$SRC/plate_skin_layers.peec" "$WORK/"
run plate_skin_layers.peec
r12=$(awk -F, 'NR==2{rdc=$3} NR==3{printf "%.9e", $3/rdc}' "$csv")
chk "slab Rdc (sandwich)" "$(getR)" 1.724138e-4 0.001
for nt in 6 3; do
	sed "s/5.8e7 4 4 12\$/5.8e7 4 4 $nt/" "$SRC/plate_skin_layers.peec" > "$WORK/plate_skin_nt$nt.peec"
	run plate_skin_nt$nt.peec
	rat=$(awk -F, 'NR==2{rdc=$3} NR==3{printf "%.9e", $3/rdc}' "$csv")
	if [ "$nt" = 6 ]; then r6=$rat; else r3=$rat; fi
done
# Richardson 補外 (O(層厚^2)) を解析値と比較。許容 3.5% の根拠は .peec を参照
chk "slab skin ratio (extr.)" "$(awk -v a="$r12" -v b="$r6" 'BEGIN{printf "%.9e", a + (a-b)/3}')" 1.897806 0.035
awk -v r3="$r3" -v r6="$r6" -v r12="$r12" 'BEGIN {
	d1 = r6 - r3; if (d1 < 0) d1 = -d1;
	d2 = r12 - r6; if (d2 < 0) d2 = -d2;
	ok = (d2 <= 0.35 * d1);
	printf "%-24s |r12-r6|=%.4f |r6-r3|=%.4f -> %s (<= 0.35)\n", "slab skin convergence", d2, d1, ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' || status=1

echo "--- panel conductors (quad / disk)"
# (x) 一般四辺形 : 正方形シートの DC 抵抗 (1 スクエアの解析値) と、
#     L の plate (矩形リボン閉形式 = 別経路の実装) との相互検証
cp "$SRC/quad_square.peec" "$WORK/"
run quad_square.peec
chk "quad sheet Rin (DC)" "$(getR)" 1.724138e-4 0.001
lq=$(getL)
cp "$SRC/quad_square_ref.peec" "$WORK/"
run quad_square_ref.peec
chk "quad sheet L vs plate" "$lq" "$(getL)" 0.005

# (y) 台形シート : R = l ln(W2/W1)/(sigma t (W2-W1)) (1 次元 + くさび補正)
cp "$SRC/quad_taper.peec" "$WORK/"
run quad_taper.peec
chk "taper Rin (DC)" "$(getR)" 4.780325e-4 0.015

# (z) 円板の対無限遠容量 : 4/8 リングから Richardson 補外 -> 8 eps0 a
cp "$SRC/disk_cap.peec" "$WORK/"
run disk_cap.peec
cd4=$(getC)
sed 's/5.8e7 4 32/5.8e7 8 32/' "$SRC/disk_cap.peec" > "$WORK/disk_cap8.peec"
run disk_cap8.peec
cd8=$(getC)
chk "disk C (extrapolated)" "$(awk -v a="$cd4" -v b="$cd8" 'BEGIN{printf "%.9e", 2*b-a}')" 7.08335e-12 0.01

# (aa) 円環の広がり抵抗 : ln(2)/(2 pi sigma t)
cp "$SRC/disk_annulus.peec" "$WORK/"
run disk_annulus.peec
chk "annulus Rin (DC)" "$(getR)" 1.902031e-5 0.03

echo "--- multi-port S parameters"
# (q) 抵抗性 T 型 2 ポート : Z=[[75,50],[50,75]], Z0=50 -> S11=1/21, S21=8/21 (実数)
cp "$SRC/tnetwork_spara.peec" "$WORK/"
run tnetwork_spara.peec
chk "T-net S11 (real)" "$(getS S11 1)" 0.047619048 0.001
chk "T-net S21 (real)" "$(getS S21 1)" 0.380952381 0.001
# 相反性 S12 = S21 と対称性 S11 = S22 (この回路では厳密に成り立つ)
chk "T-net S12 = S21"  "$(getS S12 1)" "$(getS S21 1)" 0.001
chk "T-net S22 = S11"  "$(getS S22 1)" "$(getS S11 1)" 0.001

# (r) 直列 L 入り T 型 2 ポート (wL = 50 ohm) : S が複素数になる非対称・相反回路
cp "$SRC/tnetwork_l_spara.peec" "$WORK/"
run tnetwork_l_spara.peec
chk "T-net(L) Re S11" "$(getS S11 1)"  0.101123595 0.001
chk "T-net(L) Im S11" "$(getS S11 2)"  0.561797753 0.001
chk "T-net(L) Re S21" "$(getS S21 1)"  0.359550562 0.001
chk "T-net(L) Im S21" "$(getS S21 2)" -0.224719101 0.001
chk "T-net(L) Re S22" "$(getS S22 1)"  0.056179775 0.001
chk "T-net(L) Im S22" "$(getS S22 2)"  0.089887640 0.001
# 相反性 (非対称回路でも S12 = S21)
chk "T-net(L) S12=S21 Re" "$(getS S12 1)" "$(getS S21 1)" 0.001
chk "T-net(L) S12=S21 Im" "$(getS S12 2)" "$(getS S21 2)" 0.001

# Touchstone 出力の体裁 (2 ポートは 1 行 = 周波数 + 4 成分 x 実虚 = 9 列)
if [ ! -s "$WORK/peec.s2p" ]; then
	echo "*** peec.s2p not generated" >&2
	status=1
else
	awk '/^[^!#]/ {if (NF != 9) {printf "%-24s -> NG (columns=%d)\n", "Touchstone s2p format", NF; exit 1}}
	     END {printf "%-24s -> OK\n", "Touchstone s2p format"}' "$WORK/peec.s2p" || status=1
fi

echo "--- current / charge distribution"
# (s) 単線 8 分割、容量なし : 電荷が溜まらないので全区間に同じ 1A が流れる
#     (キルヒホッフの電流則。分割数・周波数・材料に依らず厳密)
cp "$SRC/wire_dist.peec" "$WORK/"
run wire_dist.peec
awk -F, '$1 == "I" {
	n++; d = $9 - 1.0; ad = (d < 0) ? -d : d; if (ad > mx) mx = ad;
	ai = ($8 < 0) ? -$8 : $8; if (ai > mi) mi = ai;
}
END {
	ok = (n == 8) && (mx <= 1e-9) && (mi <= 1e-9);
	printf "%-24s segments=%d max||I|-1|=%.3e max|Im|=%.3e -> %s\n",
		"wire I continuity", n, mx, mi, ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' "$WORK/dist.csv" || status=1

# (t) 同じ構成に capacitance = 1 : ポートは 1A を入れて 1A を出すので
#     構造の総電荷は 0 のまま (max|q| に対する相対値で判定)
sed 's/^distribution = 1/capacitance = 1\ndistribution = 1/' "$SRC/wire_dist.peec" > "$WORK/wire_dist_cap.peec"
run wire_dist_cap.peec
awk -F, '$1 == "Q" {
	n++; sr += $7; si += $8; if ($9 > mx) mx = $9;
}
END {
	rel = (mx > 0) ? sqrt(sr * sr + si * si) / mx : 1;
	ok = (n > 0) && (rel <= 1e-6);
	printf "%-24s cells=%d |sum q|/max|q|=%.3e -> %s\n",
		"wire charge neutrality", n, rel, ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' "$WORK/dist.csv" || status=1

echo "--- thread invariance (OpenMP)"
# (u) スレッド数を変えても出力が一致すること。
#     並列化しているのは部分インダクタンス充填・電位係数充填・LU 分解の
#     残余行列更新で、いずれも要素ごとに独立 (順序依存の加算が無い) なので
#     結果はビット単位で一致する。
#     plate_cap は未知数 225 で LU の並列経路 (n-k > 64) も通る。
cp "$SRC/plate_cap.peec" "$WORK/"
(cd "$WORK" && "$PEEC" -n 1 plate_cap.peec > /dev/null)
cp "$csv" "$WORK/zin_n1.csv"
(cd "$WORK" && "$PEEC" -n 4 plate_cap.peec > /dev/null)
if cmp -s "$WORK/zin_n1.csv" "$csv"; then
	printf "%-24s -> OK (-n 1 と -n 4 が完全一致)\n" "thread invariance"
else
	printf "%-24s -> NG (-n 1 と -n 4 が不一致)\n" "thread invariance" >&2
	status=1
fi

if [ "$status" -ne 0 ]; then
	echo "*** PEEC validation FAILED" >&2
else
	echo "PEEC validation passed"
fi
exit $status
