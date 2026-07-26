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

if [ "$status" -ne 0 ]; then
	echo "*** PEEC validation FAILED" >&2
else
	echo "PEEC validation passed"
fi
exit $status
