# OpenPEEC

準静的〜フルウェーブ PEEC (部分要素等価回路法) 電気回路ソルバー (C11)。
OpenFDTD の姉妹プロジェクトで、ビルド規約・移植性規則を共有する。

集中定数 MNA + 導体形状からの部分要素抽出 (インダクタンス L / 電位係数 P /
抵抗 R) → 入力インピーダンス Zin(f)。
導体は丸線 (`wire`) / 角線 (`bar`) / 面導体 (`plate` = 矩形、
`quad` = 凸四辺形、`disk` = 円板)。ほかに無限 PEC 地板 (`groundplane`、
鏡像法)、誘電体ブリック (`dielectric`、Ruehli の過剰容量)、
遠方界後処理 (`farfield` → `far.csv`)、平面波入射 (`planewave` → `pw.csv`、
EMC イミュニティ)、過渡応答 (`transient` → `tran.csv`、掃引の逆 FFT)、
対数掃引 (`frequency ... log`)、縁寄せ格子 (`grading`)、
誘電損 tanδ / 単極 Debye 分散 (`dielectric` の省略可能引数)。
`capacitance` / `skineffect` / `retardation` / `groundplane` / `farfield` /
`planewave` / `transient` は既定で無効 (キー省略時は従来動作と完全一致)。

CSV → HDF5 の変換は `tools/peec2h5.py` (numpy + h5py)。**ソルバー本体は
外部ライブラリに依存しない**規約なので、HDF5 はこのスクリプト側に置く。

外部ライブラリに依存しない (C11 + CMake、OpenMP のみ任意)。

**`AGENTS.md`** に同じ規約を単独で読める形でまとめてある (Codex 等、
`CLAUDE.md` / `.claude/rules/` を読まないエージェント向け)。
**規約を変えたら両方直すこと。**

## ビルド / テスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 回帰 : 解析解・文献値との比較 (MNA / 部分 L / 表皮効果 / 容量 / 遅延 / 角線 /
#        面導体 / 体積セル / パネル / 地板 / 遠方界 / 誘電体 / 対数掃引 /
#        平面波入射 / 過渡応答 / 縁寄せ格子 / Debye 分散)
sh data/sample/peec_check.sh "$PWD/bin/peec" /tmp/peec-check
```

`/check` (ビルド+検証)、`/preflight` (push 前の一括点検)、
`/add-feature` (検証ケース込みの機能追加) のスラッシュコマンドがある。

## ソース構成

| ファイル | 役割 |
|---|---|
| `src/input_data.c` | `.peec` パーサ (key = value) |
| `src/wire.c` | 幾何段 : 電流セル (`p->seg`) と電荷セル (`p->chg`) を作る |
| `src/partial.c` | 幾何二重積分と部分インダクタンス (細線) |
| `src/surface.c` | 面セル (リボン) の幾何二重積分 |
| `src/volume.c` | 体積セル (矩形バー) の Hoer–Love 閉形式 (plate の厚み分割) |
| `src/polygon.c` | 多角形セル (quad/disk パネル) の幾何二重積分 (多角形閉形式) |
| `src/potential.c` | 電位係数 P と節点容量行列 C = P⁻¹ (圧縮時は P の要素評価のみ) |
| `src/skin.c` | 表皮効果 (丸線は Bessel、角線は合成式) |
| `src/mna.c` | MNA 番号付けとスタンプ |
| `src/lu.c` | 複素 LU 分解 (部分ピボット) |
| `src/iterative.c` | GMRES (acceleration = 1 の掃引 LU 再利用と compression = 1 の行列フリー) |
| `src/hmatrix.c` | Lp / P の H 行列圧縮 (クラスタツリー + ACA、compression = 1) |
| `src/precond.c` | 葉ブロック消去 + 回路 Schur 補元の前処理 (compression = 1) |
| `src/solve.c` | 周波数掃引 |
| `src/output.c` | `peec.log` の表、`zin.csv`、Touchstone `peec.sNp`、`dist.csv` |
| `src/farfield.c` | 遠方界後処理 (`farfield` → `far.csv`、D / G / 放射効率) |
| `src/transient.c` | 過渡応答 (`transient` → `tran.csv`、掃引の逆フーリエ変換) |
| `tools/peec2h5.py` | CSV → HDF5 変換 (本体の依存を増やさないための外付け) |

## 詳細な規則

作業対象のファイルに応じて `.claude/rules/` が自動で読み込まれる。
先に目を通しておくべきもの:

- `.claude/rules/portability.md` — MSVC で実際に踏んだ落とし穴
  (VLA 禁止 / OpenMP インデックス事前宣言 / `<complex.h>` 不可 など)。
  編集のたびに `.claude/hooks/check-portability.sh` が自動検査する。
- `.claude/rules/physics-invariants.md` — **壊すと結果が静かに狂う 7 つの
  不変条件**と、その番人になっている検証判定の対応。幾何積分・セル構成・
  MNA に触る前に必読。
- `.claude/rules/validation.md` — 入力キーの後方互換規則と、検証ケースの
  作り方 (期待値はコードと独立な出所にすること)。

## 設計の規則

- グローバル変数は使わない。状態は `peec_t` コンテキスト構造体 1 個を
  main で確保して関数に渡す。
- 新機能には data/sample/ の解析解付き検証ケースを追加し、
  `peec_check.sh` に判定を足す (CI 3 OS で自動実行される)。
- 新たに踏んだ移植性の落とし穴は `.claude/rules/portability.md` に、
  新たな不変条件は `.claude/rules/physics-invariants.md` に追記する。

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows (MSVC + Ninja) の 3 OS +
`sanitize` (Linux, ASan + UBSan)。検証スクリプトは全ジョブとも同一の
`data/sample/peec_check.sh` を `shell: bash` (Windows は Git Bash) で実行する。
タグ `v*` push で Release にバイナリ添付。

`sanitize` ジョブは LeakSanitizer を有効にして同じ検証を流す。**確保した
メモリは `peec_free()` (src/input_data.c) で必ず解放すること** — `peec_t` に
新しい動的配列メンバを足したらここにも足す。忘れるとこのジョブが落ちる。
密行列 (`lp` / `cmat`) が規模に比例して増えるので実用上も効く。

## 多ポート解析

`port` を複数書くと Z 行列 (ポート j に 1A 注入、他ポート開放) と S 行列
(電力波の定義、`solve.c` の `z_to_s()`) を求め、`peec.log` の表と
Touchstone `peec.sNp` に出力する。基準抵抗はポートごとに指定できるが、
Touchstone 1.1 は 1 個しか記録できないので port#1 の値を書いて注記する。
2 ポートだけ Touchstone の列順が S11 S21 S12 S22 と転置になる (仕様)。

## 並列化とスレッド数不変性

OpenMP で並列化しているのは `lp_fill` (partial.c)、`pot_fill` (potential.c)、
`lu_decomp` の残余行列更新 (lu.c)、GMRES の行列ベクトル積 (iterative.c)、
H 行列のブロック充填と matvec (hmatrix.c)、前処理の葉 LU (precond.c) の
6 箇所。**いずれも要素・行・ブロックごとに独立で、順序依存の加算
(リダクション) を持たない** (GMRES の内積・Gram-Schmidt は直列。H 行列の
matvec は葉行クラスタ = 出力の互いに素な区間ごとに並列化し、区間内の
ブロック加算順は固定) ため、スレッド数を変えても結果はビット単位で一致する。
`peec_check.sh` が `-n 1` と `-n 4` の `zin.csv` 完全一致を判定しているので、
リダクションを持つ並列化を足すとここが落ちる。
その場合は「一致する」という README の主張ごと見直すこと。

**ブロックの行区間は木の階層をまたいで重なる**ことに注意 : 行クラスタ単位で
matvec を並列化すると複数スレッドが y の同じ要素を read-modify-write して
静かに寄与を落とす (実際に踏んだ : 相対誤差 ~0.5 が実行ごとに変わる)。
並列の単位は必ず葉 ([0, n) を重複なく分割する) にすること。

MSVC の OpenMP 2.0 はループ内でのインデックス宣言を許さない (C3015) ので、
新しい `#pragma omp parallel for` を書くときはループ変数を事前宣言する。
