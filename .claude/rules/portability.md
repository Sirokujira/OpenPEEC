---
paths:
  - "src/**/*.c"
  - "include/**/*.h"
  - "CMakeLists.txt"
---

# 移植性の絶対規則 (Windows CI で実際に踏んだもの)

Linux/macOS では通るが Windows (MSVC) で落ちるものだけを挙げる。
`.claude/hooks/check-portability.sh` が編集のたびに (1)(2)(3) を自動検査する。

1. **C99 VLA 禁止** (MSVC C2057/C2466)。`malloc` + 明示インデックスの
   フラット配列を使う。`a[i * n + j]` の形にする。
2. **OpenMP for のインデックスは事前宣言する** (MSVC C3015)。MSVC の OpenMP は
   2.0 相当で `#pragma omp parallel for` の直後に `for (int i = ...)` と
   書けない。`int i;` を前置して `for (i = ...)` にする。
3. **C99 `<complex.h>` は使わない** (MSVC 非準拠)。`include/complex.h` の
   `d_complex_t` と `d_add` / `d_mul` / `d_div` / `d_inv` 等を使う。
4. **float\*/double\* の取り違え禁止**。配列の実型と読み出しポインタ型の
   不一致は Windows で 0xC0000005 クラッシュする (glibc は偶然耐える)。
5. libm リンクは CMake の `MATH_LIB` 変数経由 (Windows には m.lib が無い)。
6. MSVC フラグは CMakeLists の既存ブロックに従う
   (`/utf-8`, `_USE_MATH_DEFINES`, `_CRT_SECURE_NO_WARNINGS`)。
7. 数学・物理定数は `peec.h` の自前マクロ (`PI` / `EPS0` / `MU0` / `C0`) を使う。
8. 外部ライブラリ (LAPACK/BLAS/HDF5 等) を追加しない。OpenMP のみ任意依存で、
   `#ifdef _OPENMP` でガードする。

## 検査

```bash
# 3 OS 共通の警告掃討 (CI は -Wall なしだが、ここで 0 件にしておく)
for f in src/*.c; do gcc -std=c11 -Wall -Wextra -fopenmp -Iinclude -c "$f" -o /dev/null; done
```
