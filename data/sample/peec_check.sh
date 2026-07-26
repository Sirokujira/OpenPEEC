#!/bin/sh
# peec_check.sh — OpenPEEC 検証 (CI 用)
#
# data/sample/ の 3 ケースを実行し、zin.csv を解析解と比較する。
# 期待値の導出は各 .peec ファイルのコメント参照。
#
# 使い方 : peec_check.sh <peec 実行ファイル(絶対パス)> [作業ディレクトリ]

set -e

PEEC="$1"
WORK="${2:-.}"
SRC="$(cd "$(dirname "$0")" && pwd)"

if [ -z "$PEEC" ]; then
	echo "Usage: peec_check.sh <peec> [workdir]" >&2
	exit 2
fi

mkdir -p "$WORK"
status=0

# run <input.peec> : 実行して normal end を確認
run() {
	(cd "$WORK" && "$PEEC" -n 2 "$1" > /dev/null)
	grep -q "normal end" "$WORK/peec.log"
}

# chk <label> <actual> <expected> <tol> : 相対誤差判定
chk() {
	awk -v a="$2" -v e="$3" -v tol="$4" -v lb="$1" 'BEGIN {
		d = (a - e) / e; ad = (d < 0) ? -d : d;
		printf "%-22s actual=%.6g expected=%.6g -> %s %+.3f%%\n", lb, a, e, (ad <= tol) ? "OK" : "NG", d * 100;
		exit (ad <= tol) ? 0 : 1
	}' || status=1
}

csv="$WORK/zin.csv"

# ── (a) 直列 RLC : Rin(f0)=50, Xin(f0)~0, Xin(2f0)=47.43416 ─────────
cp "$SRC/rlc_series.peec" "$WORK/"
run rlc_series.peec
rin1=$(awk -F, 'NR==2{print $3}' "$csv")
xin1=$(awk -F, 'NR==2{print $4}' "$csv")
xin2=$(awk -F, 'NR==3{print $4}' "$csv")
chk "RLC Rin(f0)" "$rin1" 50.0 0.005
awk -v x="$xin1" 'BEGIN {
	ax = (x < 0) ? -x : x;
	printf "%-22s actual=%.6g -> %s (|Xin| <= 0.05)\n", "RLC Xin(f0)", x, (ax <= 0.05) ? "OK" : "NG";
	exit (ax <= 0.05) ? 0 : 1
}' || status=1
chk "RLC Xin(2f0)" "$xin2" 47.43416 0.005

# ── (b) 単線ワイヤ : L=1.32038e-6 H, Rin=5.48810e-3 ohm ─────────────
cp "$SRC/wire_single.peec" "$WORK/"
run wire_single.peec
L1=$(awk -F, 'NR==2{printf "%.9e", $4/(2*3.14159265358979*$2)}' "$csv")
R1=$(awk -F, 'NR==2{print $3}' "$csv")
chk "wire L (ndiv=1)" "$L1" 1.32038e-6 0.01
chk "wire Rin" "$R1" 5.48810e-3 0.01

# ndiv=8 でも同じ L になること (同一直線相互インダクタンスの回帰)
sed 's/5.8e7 1$/5.8e7 8/' "$SRC/wire_single.peec" > "$WORK/wire_ndiv8.peec"
run wire_ndiv8.peec
L8=$(awk -F, 'NR==2{printf "%.9e", $4/(2*3.14159265358979*$2)}' "$csv")
chk "wire L (ndiv=8)" "$L8" 1.32038e-6 0.01

# ── (c) 正方ループ : L=7.24689e-7 H ─────────────────────────────────
cp "$SRC/loop_square.peec" "$WORK/"
run loop_square.peec
LL=$(awk -F, 'NR==2{printf "%.9e", $4/(2*3.14159265358979*$2)}' "$csv")
chk "loop L" "$LL" 7.24689e-7 0.02

if [ "$status" -ne 0 ]; then
	echo "*** PEEC validation FAILED" >&2
else
	echo "PEEC validation passed"
fi
exit $status
