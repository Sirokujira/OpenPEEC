# OpenPEEC

準静的 PEEC (部分要素等価回路法) 電気回路ソルバー (C)。OpenFDTD の
姉妹プロジェクトで、ビルド規約・移植性規則を共有する。
集中定数 MNA + ワイヤ形状からの部分要素抽出 (インダクタンス L / 電位係数 P /
抵抗 R) → Zin(f)。P (`capacitance`) と表皮効果 (`skineffect`) は既定で無効。

## ビルド / テスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 回帰 : 解析解の全ケース (MNA / 部分 L / 表皮効果 / 容量性 PEEC) と比較
sh data/sample/peec_check.sh bin/peec /tmp/peec-check
```

## 移植性の絶対規則 (OpenFDTD の Windows CI で実際に踏んだもの)

- **C99 VLA 禁止** (MSVC C2057/C2466)。`malloc` + 明示インデックスの
  フラット配列を使う。
- **float\*/double\* の取り違え禁止**: 配列の実型と読み出しポインタ型の
  不一致は Windows で 0xC0000005 クラッシュ (glibc は偶然耐える)。
- **OpenMP for のインデックスは事前宣言する** (MSVC C3015)。MSVC の OpenMP は
  2.0 相当で、`#pragma omp parallel for` の直後に `for (int i = ...)` と
  書けない。`int i;` を前置して `for (i = ...)` にする (実際に踏んだ)。
- libm リンクは CMake の `MATH_LIB` 変数経由 (Windows には m.lib が無い)。
- MSVC フラグは CMakeLists の既存ブロックに従う
  (`/utf-8`, `_USE_MATH_DEFINES`, `_CRT_SECURE_NO_WARNINGS`)。
- 数学定数は `PI` / `EPS0` / `MU0` 等の peec.h の自前マクロを使う。
- C99 `<complex.h>` は使わない (MSVC 非準拠)。include/complex.h の
  `d_complex_t` を使う。

## 設計の規則

- グローバル変数は使わない。状態は `peec_t` コンテキスト構造体 1 個を
  main で確保して関数に渡す。
- 入力キー追加は `src/input_data.c` に、既定値は「キー省略時に従来動作と
  完全一致」になるよう初期化する (後方互換)。未知キーは無視 (前方互換)。
- 新機能には data/sample/ の解析解付き検証ケースを追加し、
  `peec_check.sh` に判定を足す (CI 3 OS で自動実行される)。
- OpenMP は任意依存。`#ifdef _OPENMP` でガードし、スレッド数によらず
  出力がビット一致することを確認する (現状 lp_fill / pot_fill が並列)。
  並列ループ内で共有配列に `+=` しない (pot_fill は積分結果を一時配列に
  出してから直列で集約している)。
- 外部ライブラリ (LAPACK/BLAS/HDF5 等) を追加しない。

### 幾何カーネルの不変条件 (壊さないこと)

`neumann_pair()` / `neumann_self()` は同一の細線カーネル
(距離を導体半径 a で下限打ち切り) を使う。この統一により
`Σᵢⱼ Iᵢⱼ = I_self(全長)` が厳密に成立し、

- 直線ワイヤの合計インダクタンスが分割数に依存しない
- 対無限遠容量が平均電位法の解析解と厳密に一致する

という 2 つの性質が保証される。`peec_check.sh` の
"wire L (ndiv=1)/(ndiv=8)" と "wire Ctotal" が許容 0.1% でこれを守っている。
片方の項だけ正則化を変えると両方が静かに劣化するので注意。

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows (MSVC + Ninja)。
検証スクリプトは 3 OS とも同一の `data/sample/peec_check.sh` を
`shell: bash` (Windows は Git Bash) で実行する。タグ `v*` push で
Release にバイナリ添付。
