#!/bin/sh
# TEMPORARY diagnostic for the GCC 15 -O3 failure on PR #24 (to be removed)
set -u
R=$PWD; S=$R/data/sample; W=/tmp/diag; rm -rf $W; mkdir -p $W/obj
gcc --version | head -1
gcc -Q --help=target 2>/dev/null | grep -E "^ +-march=|^ +-mtune=" | head -2
echo | gcc -dM -E - | grep -E "__FMA__|__AVX2__|__x86_64__" | tr '\n' ' '; echo

# probe <binary> <label>
probe() {
	B=$1; L=$2
	d=$W/run_$L; rm -rf $d; mkdir -p $d; cd $d || return
	cp $S/pw_loop_rcs.peec .; "$B" -n 2 pw_loop_rcs.peec >/dev/null 2>&1
	y=$(awk -F, 'NR > 1 && ($2 + 0 == 90) && ($3 + 0 == 90) {print $8}' farpw.csv 2>/dev/null)
	voc1=$(sed -n 2p pw.csv | cut -d, -f3,4)
	grep -v "^farfield" pw_loop_rcs.peec > nf.peec; "$B" -n 2 nf.peec >/dev/null 2>&1
	voc2=$(sed -n 2p pw.csv | cut -d, -f3,4)
	cp $S/bent_pw.peec .; "$B" -n 2 bent_pw.peec >/dev/null 2>&1; cp far.csv far_bent.csv
	out=""
	for spec in "30 45 1" "90 30 2"; do
		sed "s/^planewave = .*/planewave = $spec 1.0/" bent_pw.peec > t.peec; "$B" -n 2 t.peec >/dev/null 2>&1
		set -- $spec
		dd=$(awk -F, -v th="$1" -v ph="$2" -v pol="$3" 'BEGIN {c = 2.99792458e8; lam = c / 1.5e8; eta = 4 * atan2(0, -1) * 1e-7 * c; fac = 2 * lam / eta}
		NR == FNR && FNR > 1 && (($2 + 0) == th) && (($3 + 0) == ph) { if (pol == 1) {a = $4; b = $5} else {a = $6; b = $7}; got = 1; next }
		FNR == 2 && got { er = fac * b; ei = -fac * a; dr = $3 - er; di = $4 - ei; printf "%.3e(Voc=%s,%s exp=%.6e,%.6e)", sqrt((dr * dr) + (di * di)) / (sqrt((er * er) + (ei * ei)) + 1e-12), $3, $4, er, ei }' far_bent.csv pw.csv)
		out="$out $dd"
	done
	grep -v "^farfield" t.peec > t2.peec; "$B" -n 2 t2.peec >/dev/null 2>&1
	vnf=$(sed -n 2p pw.csv | cut -d, -f3,4)
	cd $R
	echo "[$L] loop sig(y)=$y (exp 1.1814e-09) Voc(ff)=$voc1 Voc(noff)=$voc2"
	echo "[$L]   recip 30/45/1, 90/30/2 :$out  ; 90/30/2 Voc without farfield=$vnf"
}

# build <label> <default flags> [file=flags ...]
build() {
	lb=$1; df=$2; shift 2
	rm -f $W/obj/*.o
	for f in src/*.c; do
		fl="$df"
		for ov in "$@"; do
			case "$ov" in "$(basename $f)="*) fl="${ov#*=}";; esac
		done
		gcc -std=gnu11 -DNDEBUG $fl -Iinclude -c $f -o $W/obj/$(basename $f .c).o || echo "compile fail $f"
	done
	gcc $df $W/obj/*.o -o $W/peec_$lb -lm || echo "link fail $lb"
	probe $W/peec_$lb $lb
}

probe $R/bin/peec cmake
build O3 "-O3 -fopenmp"
build O2 "-O2 -fopenmp"
build O3noomp "-O3"
build O3novec "-O3 -fopenmp -fno-tree-vectorize"
build O3noipa "-O3 -fopenmp -fno-ipa-sra -fno-ipa-cp -fno-ipa-modref"
build O3noinl "-O3 -fopenmp -fno-inline"
build O3nofma "-O3 -fopenmp -ffp-contract=off"
build O3solve0 "-O3 -fopenmp" "solve.c=-O0 -fopenmp"
build O3far0 "-O3 -fopenmp" "farfield.c=-O0 -fopenmp"
build O3mna0 "-O3 -fopenmp" "mna.c=-O0 -fopenmp"
build O3main0 "-O3 -fopenmp" "main.c=-O0 -fopenmp"
build O3lu0 "-O3 -fopenmp" "lu.c=-O0 -fopenmp"
build O3part0 "-O3 -fopenmp" "partial.c=-O0 -fopenmp"
