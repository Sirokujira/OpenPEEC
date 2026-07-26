#!/bin/sh
# check-portability.sh — 編集後に走る移植性チェック (PostToolUse hook)
#
# Windows CI で実際に踏んだ落とし穴だけを、grep で確実に検出できる形で見る。
# 引数・標準入力には依存せず src/ 全体を走査するので、3 OS どこでも動く。
# 違反を見つけたら stderr に出して exit 2 (Claude にフィードバックされる)。

root="$(cd "$(dirname "$0")/../.." && pwd)"
status=0

# (1) MSVC の OpenMP は 2.0 相当 : #pragma omp parallel for の直後の for 文で
#     インデックスを宣言できない (error C3015)。ループ変数は事前宣言する。
hit=$(grep -n -A2 '^[[:space:]]*#pragma omp parallel for' "$root"/src/*.c 2>/dev/null \
	| grep -E 'for[[:space:]]*\([[:space:]]*(int|long|size_t|unsigned)[[:space:]]')
if [ -n "$hit" ]; then
	echo "[hook] MSVC C3015: '#pragma omp parallel for' の直後の for 文で" >&2
	echo "       ループ変数を宣言しています。'int i;' を前置して 'for (i = ...)' にしてください。" >&2
	echo "$hit" >&2
	status=2
fi

# (2) C99 VLA 禁止 (MSVC C2057/C2466)。malloc + フラット配列を使う。
#     配列サイズが定数でも大文字マクロでもないものを拾う。
hit=$(grep -nE '^[[:space:]]*(const[[:space:]]+)?(char|int|long|float|double|size_t)[[:space:]]+[a-z_][a-zA-Z0-9_]*\[[a-z_][a-zA-Z0-9_]*\][[:space:]]*;' \
	"$root"/src/*.c "$root"/include/*.h 2>/dev/null)
if [ -n "$hit" ]; then
	echo "[hook] C99 VLA の可能性 (MSVC C2057/C2466)。malloc + 明示インデックスの" >&2
	echo "       フラット配列に置き換えてください (サイズが定数式なら無視して構いません)。" >&2
	echo "$hit" >&2
	status=2
fi

# (3) C99 <complex.h> は MSVC 非準拠。include/complex.h の d_complex_t を使う。
hit=$(grep -n '#include[[:space:]]*<complex\.h>' "$root"/src/*.c "$root"/include/*.h 2>/dev/null)
if [ -n "$hit" ]; then
	echo "[hook] C99 <complex.h> は MSVC で使えません。include/complex.h の d_complex_t を使ってください。" >&2
	echo "$hit" >&2
	status=2
fi

exit $status
