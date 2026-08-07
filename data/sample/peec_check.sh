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

# (a2) 対数掃引 : f0/10 .. 10 f0 を 2 分割 -> 幾何中点がちょうど f0 になる。
#      X(u f0) = sqrt(L/C) (u - 1/u) : X(0.1 f0) = -313.06549, X(10 f0) = +313.06549
sed 's|^frequency = .*|frequency = 5.0329212e5 5.0329212e7 2 log|' \
	"$SRC/rlc_series.peec" > "$WORK/rlc_log.peec"
run rlc_log.peec
chk "RLC log Xin(f0/10)" "$(awk -F, 'NR==2{print $4}' "$csv")" -313.06549 0.005
chk "RLC log Rin(mid=f0)" "$(awk -F, 'NR==3{print $3}' "$csv")" 50.0 0.005
awk -F, 'NR==3 {
	x = $4; ax = (x < 0) ? -x : x;
	printf "%-24s actual=%.6g -> %s (|Xin| <= 0.05)\n", "RLC log Xin(mid=f0)", x, (ax <= 0.05) ? "OK" : "NG";
	exit (ax <= 0.05) ? 0 : 1
}' "$csv" || status=1
chk "RLC log Xin(10f0)" "$(awk -F, 'NR==4{print $4}' "$csv")" 313.06549 0.005

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

# (an) 遅延 x 面セル : ストリップダイポール。遅延ありの検証はこれまで全て
#      細線だったので、面セル (リボン) 経路にも番人を置く。
#      教科書値 (l = 0.478 lambda、Rin = 73.1 ohm) と、同じ帯を bar
#      (等価半径の細線 = 独立経路) でモデル化したものとの相互検証。
# resonance <label> : Xin の零交差を線形補間して "l/lambda Rin" を返す
resonance() {
	awk -F, 'NR>1 {
		f = $2; r = $3; x = $4;
		if (pf && (px < 0) && (x >= 0) && !done) {
			u = -px / (x - px);
			printf "%.9e %.9e", (pf + u*(f-pf)) / 2.99792458e8, pr + u*(r-pr);
			done = 1;
		}
		pf = f; pr = r; px = x
	}' "$csv"
}
cp "$SRC/strip_dipole_ret.peec" "$WORK/"
run strip_dipole_ret.peec
sp=$(resonance)
cp "$SRC/strip_dipole_ret_bar.peec" "$WORK/"
run strip_dipole_ret_bar.peec
sb=$(resonance)
if [ -z "$sp" ] || [ -z "$sb" ]; then
	echo "*** no resonance in strip dipole sweep (plate='$sp' bar='$sb')" >&2
	status=1
else
	chk "strip dip l/lambda" "${sp% *}" 0.478 0.05
	chk "strip dip Rin at res." "${sp#* }" 73.1 0.10
	chk "strip dip f vs bar" "${sp% *}" "${sb% *}" 0.03
	chk "strip dip Rin vs bar" "${sp#* }" "${sb#* }" 0.05
fi

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
#     許容 0.1% : 同じ格子を多角形経路とリボン経路で解いた結果なので実測は
#     1e-5 レベルで一致する。ここを締めておくと、両経路の外側求積を距離で
#     粗くする最適化 (ribbon_nq / POLY_R7) が壊れたときに検出できる。
cp "$SRC/quad_square.peec" "$WORK/"
run quad_square.peec
chk "quad sheet Rin (DC)" "$(getR)" 1.724138e-4 0.001
lq=$(getL)
cp "$SRC/quad_square_ref.peec" "$WORK/"
run quad_square_ref.peec
chk "quad sheet L vs plate" "$lq" "$(getL)" 0.001

# 表皮効果ありの高周波抵抗も plate と一致すること。パネルセルの断面は plate と
# 同じく格子の双対幅から幾何的に取る (不変条件 7) ので、合成式が同じ値を返す。
# 適用漏れがあると厚み方向の表皮効果 (delta < thick) が丸ごと落ちて R が
# 過小になる (1 GHz・厚さ 0.1 mm の銅シートで実測 -24% だった)。
for c in quad_square quad_square_ref; do
	sed -e 's/^title = /skineffect = 1\ntitle = /' -e 's/^frequency = .*/frequency = 1e9 1e9 0/' \
		"$SRC/$c.peec" > "$WORK/sk_$c.peec"
	run "sk_$c.peec"
	if [ "$c" = quad_square ]; then rq=$(getR); else rp=$(getR); fi
done
chk "quad sheet R vs plate (HF)" "$rq" "$rp" 0.001
# 合成式の DC 極限は幾何断面の抵抗に厳密に一致する (低周波で不連続が無いこと)
sed 's/^title = /skineffect = 1\ntitle = /' "$SRC/quad_square.peec" > "$WORK/sk_dc.peec"
run sk_dc.peec
chk "quad sheet R (skin at DC)" "$(getR)" 1.724138e-4 0.001

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

echo "--- ground plane (image method)"
# (ab) 地板上の水平単線 : L = L_self - M(2h) (Grover)、R は鏡像の影響を受けない
cp "$SRC/wire_gp.peec" "$WORK/"
run wire_gp.peec
chk "wire-gp L (images)" "$(getL)" 5.953664e-7 0.001
chk "wire-gp Rin (DC)" "$(getR)" 5.48810e-3 0.005
# 対地板容量 (鏡像電荷 -q、ndiv = 1 で平均電位法と厳密に一致)
awk '{sub(/5.8e7 8$/, "5.8e7 1"); print} /^title = /{print "capacitance = 1"}' \
	"$SRC/wire_gp.peec" > "$WORK/wire_gp_cap.peec"
run wire_gp_cap.peec
chk "wire-gp C (to plane)" "$(getC)" 1.868849e-11 0.001

# (ac) λ/4 モノポール : イメージ理論により Zin = ダイポール/2 が
#      離散化レベルで厳密に成り立つ (許容は丸め誤差 + LU 順序差のマージン)
cp "$SRC/dipole_halfwave.peec" "$WORK/"
run dipole_halfwave.peec
cp "$csv" "$WORK/zin_dip.csv"
cp "$SRC/monopole_gp.peec" "$WORK/"
run monopole_gp.peec
paste -d, "$WORK/zin_dip.csv" "$csv" | awk -F, 'NR > 1 {
	dr = $10 - ($3 / 2); di = $11 - ($4 / 2);
	e = sqrt((dr * dr) + (di * di)) / ((sqrt(($3 * $3) + ($4 * $4)) / 2) + 1);
	if (e > mx) mx = e; n++
} END {
	ok = (n > 0) && (mx <= 1e-6);
	printf "%-24s rows=%d max|Zm-Zd/2|rel=%.3e -> %s\n", "monopole = dipole/2", n, mx, ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' || status=1
# 文献値アンカー : 共振で Rin = 36.55 ohm、共振周波数はダイポールと同一
res=$(awk -F, 'NR>1 {
	f = $2; r = $3; x = $4;
	if (pf && (px < 0) && (x >= 0) && !done) {
		u = -px / (x - px);
		printf "%.9e %.9e", (pf + u*(f-pf)) / 2.99792458e8, pr + u*(r-pr);
		done = 1;
	}
	pf = f; pr = r; px = x;
}' "$csv")
if [ -z "$res" ]; then
	echo "*** no resonance (Xin zero crossing) in monopole_gp sweep" >&2
	status=1
else
	chk "monopole res. f/c" "${res% *}" 0.478 0.05
	chk "monopole Rin at res." "${res#* }" 36.55 0.10
fi

echo "--- far field"
# far.csv (単一周波数) から最大 D [linear] を取り出す
getD() { awk -F, 'NR>1 {d = exp($8/10*log(10)); if (d > mx) mx = d} END {printf "%.6f", mx}' "$WORK/far.csv"; }
# peec.log の farfield 行から放射効率を取り出す
getEff() { grep "^farfield :" "$log" | tail -1 | sed 's/.*eff = \([0-9.eE+-]*\).*/\1/'; }

# (ad) 半波長ダイポール : D = 1.628 (正弦電流 kh = 1.5065 の解析値、
#      D = 2 max|F|^2/∫|F|^2 sin th dth, F = [cos(kh cos th)-cos(kh)]/sin th)
#      効率 ~ 1 : 遠方界の規格化 + 球面積分がポインティングの定理と整合すること
cp "$SRC/dipole_ff.peec" "$WORK/"
run dipole_ff.peec
dd=$(getD)
chk "dipole ff D" "$dd" 1.628 0.015
chk "dipole ff efficiency" "$(getEff)" 1.0 0.005

# (ae) 微小ダイポール : 軸方向電流なら分布に依らずパターンは sin^2 th で D = 1.5 (厳密)
awk '{print} /^title = /{print "farfield = 36 24"}' "$SRC/dipole_short.peec" > "$WORK/dipole_short_ff.peec"
run dipole_short_ff.peec
chk "short dipole D" "$(getD)" 1.5 0.005

# (af) モノポール (地板) : 上半球に同じ界・電力は半分なので D = 2 x ダイポール。
#      鏡像電流の遠方界と上半球積分の判定 (ダイポール側 (ad) と同一周波数で比較)
awk '{sub(/^frequency = .*/, "frequency = 1.4376e8 1.4376e8 0"); print}
     /^title = /{print "farfield = 18 24"}' "$SRC/monopole_gp.peec" > "$WORK/monopole_ff.peec"
run monopole_ff.peec
chk "monopole D = 2 x dipole" "$(awk -v m="$(getD)" -v d="$dd" 'BEGIN{printf "%.6f", m/d}')" 2.0 0.01
chk "monopole ff efficiency" "$(getEff)" 1.0 0.005

echo "--- plane wave incidence (external field excitation)"
# pw.csv の 1 行目 (port#1、単一周波数) から |Voc| と実効長を取り出す
getVoc() { awk -F, 'NR==2{print $5}' "$WORK/pw.csv"; }
getLeff() { awk -F, 'NR==2{print $6}' "$WORK/pw.csv"; }

# (ah) 半波長ダイポール : |l_eff| = lambda/pi = 0.663793 m (正弦電流分布の教科書値)
cp "$SRC/dipole_pw.peec" "$WORK/"
run dipole_pw.peec
chk "dipole pw leff" "$(getLeff)" 0.66379304 0.01

# (ai) 微小ループ : |Voc| = k E0 A = 8.383380e-3 V (ファラデーの法則 = 磁界結合。
#      ダイポールの電界結合とは別経路の判定)
cp "$SRC/loop_pw.peec" "$WORK/"
run loop_pw.peec
chk "loop pw |Voc| (kE0A)" "$(getVoc)" 8.383380e-3 0.005
# 磁界結合は E と 90 deg ずれるので Voc はほぼ純虚数 (残差は O((ks)^2))
awk -F, 'NR==2 {
	ar = ($3 < 0) ? -$3 : $3; rel = ar / $5;
	printf "%-24s |Re Voc|/|Voc|=%.4f -> %s (<= 0.05)\n", "loop pw quadrature", rel, (rel <= 0.05) ? "OK" : "NG";
	exit (rel <= 0.05) ? 0 : 1
}' "$WORK/pw.csv" || status=1

# 誘起電流 : ループを 50 ohm で閉じると |I| = |Voc| / |R + R_wire + j wL|。
# 期待値はファラデーの法則と検証済みの L (loop_square) の組合せで、
# 誘起電流の実装とは独立。直列ループなので 4 区間の電流は等しい。
sed 's/^planewave = .*/planewave = 90 0 2 1.0\ndistribution = 1\nresistor = 1 5 50/' \
	"$SRC/loop_pw.peec" > "$WORK/loop_pw_load.peec"
run loop_pw_load.peec
awk -F, '$1 == "Ipw" {
	n++; if ($10 > mx) mx = $10; if ((mn == 0) || ($10 < mn)) mn = $10;
} END {
	spread = (mx > 0) ? ((mx - mn) / mx) : 1;
	ok = (n == 4) && (spread <= 1e-9);
	printf "%-24s segments=%d spread=%.3e -> %s\n", "loop pw I continuity", n, spread, ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' "$WORK/dist.csv" || status=1
chk "loop pw induced |I|" "$(awk -F, '$1 == "Ipw" {print $10; exit}' "$WORK/dist.csv")" 1.239605e-4 0.005

# (aj) 相反定理 : 送信パス (farfield の放射ベクトル) と受信パス (planewave の
#      起電力) は独立な実装だが、離散化した系でも次の恒等式が厳密に成り立つ :
#        l_eff = rE * 2 lambda / (-j eta0),   Voc = -E0 (e^・l_eff)
#      (符号は mna_rhs_port の +1A 注入と Voc = v(n1)-v(n2) の向きの差)
#      到来方向・偏波を変えて突き合わせることで、位相規約・sinc 積分・
#      偏波ベクトル・角度依存が同時に検証される。
# recip <入力.peec のパス> <far.csv> <f[Hz]> <theta> <phi> <pol>
recip() {
	sed "s/^planewave = .*/planewave = $4 $5 $6 1.0/" "$1" > "$WORK/recip_t.peec"
	run recip_t.peec
	awk -F, -v th="$4" -v ph="$5" -v pol="$6" -v fr="$3" -v lb="$1" '
	BEGIN {c = 2.99792458e8; lam = c / fr; eta = 4 * atan2(0, -1) * 1e-7 * c; fac = 2 * lam / eta}
	NR == FNR && FNR > 1 && (($2 + 0) == th) && (($3 + 0) == ph) {
		if (pol == 1) {a = $4; b = $5} else {a = $6; b = $7}
		got = 1; next
	}
	FNR == 2 && got {
		# E0 = 1 : Re = fac*Im(rE)、Im = -fac*Re(rE)
		er = fac * b; ei = -fac * a;
		dr = $3 - er; di = $4 - ei;
		d = sqrt((dr * dr) + (di * di)) / (sqrt((er * er) + (ei * ei)) + 1e-12);
		ok = (d <= 1e-6);
		printf "%-24s rel.diff=%.3e -> %s\n", sprintf("recip %d/%d/p%d", th, ph, pol), d, ok ? "OK" : "NG";
		exit ok ? 0 : 1
	}
	END {if (!got) {printf "*** far.csv row (%d, %d) not found for %s\n", th, ph, lb; exit 1}}
	' "$2" "$WORK/pw.csv" || status=1
}

# 非対称な折れ曲がりダイポール : パターンが theta / phi 両成分を持つ
cp "$SRC/bent_pw.peec" "$WORK/"
run bent_pw.peec
cp "$WORK/far.csv" "$WORK/far_bent.csv"
for spec in "30 45 1" "60 120 1" "90 30 2" "120 255 2" "45 90 2"; do
	# shellcheck disable=SC2086
	recip "$SRC/bent_pw.peec" "$WORK/far_bent.csv" 1.5e8 $spec
done

# (ak) 地板ありの相反 : 送信側の鏡像 (farfield) と受信側の PEC 反射波
#      (planewave) が整合していることの判定。上半球のみ (theta <= 90)
sed -e 's/^frequency = .*/frequency = 1.4376e8 1.4376e8 0/' \
    -e 's/^title = .*/title = monopole receiving\nfarfield = 18 24\nplanewave = 60 0 1 1.0/' \
    "$SRC/monopole_gp.peec" > "$WORK/mono_pw.peec"
run mono_pw.peec
cp "$WORK/far.csv" "$WORK/far_mono.csv"
recip "$WORK/mono_pw.peec" "$WORK/far_mono.csv" 1.4376e8 60 0 1

echo "--- transient (inverse FFT of the sweep)"
# (al) 周波数に依らない反射係数 : y(t) = S11 x(t) が全時刻で厳密に成り立つ。
#      励振と応答を同じ合成式で作るので帯域打ち切りも相殺する。
cp "$SRC/tdr_resistor.peec" "$WORK/"
for R in 50 150 1e-6 1e9; do
	sed "s/^resistor = .*/resistor = 1 0 $R/" "$SRC/tdr_resistor.peec" > "$WORK/tdr_r.peec"
	run tdr_r.peec
	awk -F, -v R="$R" 'NR > 1 && $1 == "X" {x[$4] = $5}
	NR > 1 && $1 == "S" {s[$4] = $5}
	END {
		s11 = (R - 50) / (R + 50);
		for (t in x) {d = s[t] - (s11 * x[t]); ad = (d < 0) ? -d : d; if (ad > mx) mx = ad}
		ok = (mx <= 1e-9);
		printf "%-24s max|y-S11*x|=%.3e -> %s\n", sprintf("tdr R=%s (S11=%+.1f)", R, s11), mx, ok ? "OK" : "NG";
		exit ok ? 0 : 1
	}' "$WORK/tran.csv" || status=1
done

# 励振パルスが解析的なガウス (ピーク 1、t0 = 4 sigma) になっていること
run tdr_resistor.peec
awk -F, 'BEGIN {sig = 0.483008 / 6.4e9; t0 = 4 * sig}
$1 == "X" {
	g = exp(-($4 - t0) * ($4 - t0) / (2 * sig * sig));
	d = $5 - g; ad = (d < 0) ? -d : d; if (ad > mx) mx = ad
} END {
	ok = (mx <= 5e-3);
	printf "%-24s max|x-gaussian|=%.3e -> %s (<= 5e-3)\n", "transient excitation", mx, ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' "$WORK/tran.csv" || status=1

# 時間軸 : dt = 1/(2 fmax)、サンプル数 2N
awk -F, '$1 == "X" {n++; if (n == 2) dt = $4} END {
	ok = (n == 128) && ((dt - 7.8125e-11) < 1e-16) && ((7.8125e-11 - dt) < 1e-16);
	printf "%-24s samples=%d dt=%.6e -> %s\n", "transient time axis", n, dt, ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' "$WORK/tran.csv" || status=1

# 対数掃引では時間軸が定義できないので拒否されること
sed 's/^frequency = .*/frequency = 1e8 6.4e9 63 log/' "$SRC/tdr_resistor.peec" > "$WORK/tdr_log.peec"
if (cd "$WORK" && "$PEEC" -n 1 tdr_log.peec > tdr_log.out 2>&1); then
	printf "%-24s -> NG (log sweep must be rejected)\n" "transient sweep guard" >&2
	status=1
else
	printf "%-24s -> OK (log sweep rejected)\n" "transient sweep guard"
fi

# (am) 微小ループ : v(t) = -(E0 A/c) dx/dt (Voc ∝ jω = 微分)。
#      合成が exp(+jwt) 規約で一貫していることの番人 (純抵抗の判定は
#      S11 が実数なので共役の取り違えを検出できない)
cp "$SRC/loop_tran.peec" "$WORK/"
run loop_tran.peec
awk -F, 'BEGIN {
	PI = atan2(0, -1); fmax = 1e7; att = 80;
	sig = sqrt(att * log(10) / (20 * 2 * PI * PI)) / fmax;   # transient.c と同じ式
	t0 = 4 * sig; K = 0.04 / 2.99792458e8                    # K = E0 A / c
}
$1 == "V" {
	# v(t) = -(E0 A/c) dx/dt,  dx/dt = -((t-t0)/sig^2) x(t)
	a = K * (($4 - t0) / (sig * sig)) * exp(-($4 - t0) * ($4 - t0) / (2 * sig * sig));
	d = $5 - a; ad = (d < 0) ? -d : d; if (ad > mx) mx = ad;
	aa = (a < 0) ? -a : a; if (aa > pk) pk = aa
} END {
	rel = mx / pk;
	ok = (rel <= 0.01);
	printf "%-24s max|v-(-(E0A/c)dx/dt)|/peak=%.4f -> %s (<= 0.01)\n", "loop transient (d/dt)", rel, ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' "$WORK/tran.csv" || status=1

echo "--- dielectric (excess capacitance)"
# (ag) 平行平板 + εr = 4 ブリック : ΔC = (εr-1) eps0 A/d = 1.06250e-10 F。
#      平板は等電位面なので過剰容量ネットワークの合計は格子に依らず厳密
#      (縁効果は差で相殺、実測残差 ~1e-6)
getCzin() { awk -F, 'NR==2{printf "%.9e", -1/(6.283185307179586e6*$4)}' "$csv"; }
cp "$SRC/diel_pp.peec" "$WORK/"
run diel_pp.peec
cdiel=$(getCzin)
cp "$csv" "$WORK/zin_diel_lossless.csv"
sed '/^dielectric/d' "$SRC/diel_pp.peec" > "$WORK/diel_air.peec"
run diel_air.peec
cair=$(getCzin)
cp "$csv" "$WORK/zin_diel_air.csv"
chk "dielectric dC (pp)" "$(awk -v a="$cdiel" -v b="$cair" 'BEGIN{printf "%.9e", a-b}')" 1.0625025e-10 0.005
# εr = 1 のブリックはセルを作らない = ブリック無しと bit 単位で一致
sed 's/ 4 8 8 2$/ 1 8 8 2/' "$SRC/diel_pp.peec" > "$WORK/diel_eps1.peec"
run diel_eps1.peec
if cmp -s "$WORK/zin_diel_air.csv" "$csv"; then
	printf "%-24s -> OK (epsr = 1 はブリック無しと完全一致)\n" "dielectric epsr=1 noop"
else
	printf "%-24s -> NG (epsr = 1 がブリック無しと不一致)\n" "dielectric epsr=1 noop" >&2
	status=1
fi

# 誘電正接 : epsr* = epsr (1 - j tand) の損失分は枝のコンダクタンスになる。
# 損失性平行平板の解析値 G = w eps0 epsr tand A/d
#   = 2pi*1e6 * 8.8541878e-12 * 4 * 0.02 * 0.0016/4e-4 = 1.780240e-5 S
# 真空部もフリンジも無損失なので Re(Yin) はこの値そのものになる
# (導体損の寄与は Re(Y) 換算 9e-10 S で 4 桁下)。コードとは独立な期待値。
sed 's/ 4 8 8 2$/ 4 8 8 2 0.02/' "$SRC/diel_pp.peec" > "$WORK/diel_tand.peec"
run diel_tand.peec
chk "dielectric G (tand)" \
	"$(awk -F, 'NR==2{d=($3*$3)+($4*$4); printf "%.9e", $3/d}' "$csv")" 1.780240e-5 0.005
# tand = 0 を明示しても従来 (省略時) と bit 単位で一致すること
sed 's/ 4 8 8 2$/ 4 8 8 2 0/' "$SRC/diel_pp.peec" > "$WORK/diel_tand0.peec"
run diel_tand0.peec
if cmp -s "$WORK/zin_diel_lossless.csv" "$csv"; then
	printf "%-24s -> OK (tand = 0 は省略時と完全一致)\n" "dielectric tand=0 noop"
else
	printf "%-24s -> NG (tand = 0 が省略時と不一致)\n" "dielectric tand=0 noop" >&2
	status=1
fi

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
awk -F, '($1 == "I") && ($2 == 1) {
	n++; d = $10 - 1.0; ad = (d < 0) ? -d : d; if (ad > mx) mx = ad;
	ai = ($9 < 0) ? -$9 : $9; if (ai > mi) mi = ai;
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
awk -F, '($1 == "Q") && ($2 == 1) {
	n++; sr += $8; si += $9; if ($10 > mx) mx = $10;
}
END {
	rel = (mx > 0) ? sqrt(sr * sr + si * si) / mx : 1;
	ok = (n > 0) && (rel <= 1e-6);
	printf "%-24s cells=%d |sum q|/max|q|=%.3e -> %s\n",
		"wire charge neutrality", n, rel, ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' "$WORK/dist.csv" || status=1

# (ao) 多ポート : dist.csv はポートごとに 1 組出る (Z 行列の第 j 列)。
#      独立な 2 線に 1 ポートずつ置くと、ポート j 励振では線 j に厳密に 1A、
#      他方は開放なので厳密に 0A になる (KCL のみで決まり、分割数・周波数・
#      部分要素の値に依存しない)。port 列の取り違えがあれば即座に落ちる。
cp "$SRC/twoport_dist.peec" "$WORK/"
run twoport_dist.peec
awk -F, '$1 == "I" {
	n++; want = (((($2 + 0) == 1) && (($3 + 0) < 4)) ||
	             ((($2 + 0) == 2) && (($3 + 0) >= 4))) ? 1 : 0;
	d = $10 - want; ad = (d < 0) ? -d : d; if (ad > mx) mx = ad;
}
END {
	ok = (n == 16) && (mx <= 1e-12);
	printf "%-24s rows=%d max||I|-expected|=%.3e -> %s\n",
		"two-port I per port", n, mx, ok ? "OK" : "NG";
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

echo "--- HDF5 converter (tools/peec2h5.py)"
# ソルバー本体は依存ゼロなので、変換スクリプトは numpy + h5py がある環境
# でだけ検証する (無ければスキップ : CI の 3 OS には入れない)。
# 専用ディレクトリで実行して他ケースの残骸を拾わないようにし、変換結果が
# 元の CSV と整合すること (S11 が Zin と一致すること) まで確かめる。
H5DIR="$WORK/h5"
mkdir -p "$H5DIR"
cp "$SRC/tdr_resistor.peec" "$H5DIR/"
(cd "$H5DIR" && "$PEEC" -n 1 tdr_resistor.peec > /dev/null)
grep -q "normal end" "$H5DIR/peec.log"
if python3 -c "import numpy, h5py" 2>/dev/null; then
	if python3 "$SRC/../../tools/peec2h5.py" "$H5DIR" -o "$H5DIR/peec.h5" > /dev/null 2>&1 &&
	   python3 - "$H5DIR/peec.h5" <<'PYEOF'
import sys, h5py, numpy as np
with h5py.File(sys.argv[1]) as f:
    for name in ("frequency", "zin", "transient/time", "transient/excitation"):
        assert name in f, "missing dataset: " + name
    z = f["zin"][0]
    s = f["sparam"][:, 0, 0]
    err = np.abs(s - ((z - 50) / (z + 50))).max()
    assert err < 1e-6, "S11 vs Zin mismatch: %g" % err
    t = f["transient/time"][:]
    assert len(t) > 1 and np.all(np.diff(t) > 0), "bad time axis"
PYEOF
	then
		printf "%-24s -> OK (h5 written and consistent with the CSV)\n" "peec2h5 converter"
	else
		printf "%-24s -> NG (conversion or consistency check failed)\n" "peec2h5 converter" >&2
		status=1
	fi
else
	printf "%-24s -> SKIP (numpy/h5py not installed)\n" "peec2h5 converter"
fi

if [ "$status" -ne 0 ]; then
	echo "*** PEEC validation FAILED" >&2
else
	echo "PEEC validation passed"
fi
exit $status
