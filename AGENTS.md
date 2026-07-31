# AGENTS.md — OpenPEEC

準静的〜フルウェーブ PEEC (部分要素等価回路法) 電気回路ソルバー (C11)。
[OpenFDTD](https://github.com/Sirokujira/OpenFDTD) の姉妹プロジェクトで、
ビルド規約・移植性規則を共有する。**外部数値ライブラリに依存しない**
(C11 + CMake、OpenMP のみ任意)。

集中定数 MNA + 導体形状からの部分要素抽出 (インダクタンス L / 電位係数 P /
抵抗 R) → 入力インピーダンス Zin(f)、Z / S 行列、電流・電荷分布。
導体は丸線 (`wire`) / 角線 (`bar`) / 面導体 (`plate`)。`plate` は `分割数t`
(省略時 1) を 2 以上にすると厚み方向にも分割し、体積セル (矩形バー) として
厚み方向の電流分布 (表皮効果) を解く。

> このファイルは Claude Code 用の `CLAUDE.md` / `.claude/rules/*.md` と同じ内容を
> 単独で読めるようまとめたもの。**片方だけ直すと食い違うので、規約を変えたら
> 両方直すこと。**

## ビルドとテスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 検証 (解析解・文献値との比較、50 件)
sh data/sample/peec_check.sh "$PWD/bin/peec" /tmp/peec-check

# メモリ健全性 (CI の sanitize ジョブと同じ)
cmake -B build-san -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_FLAGS="-g -fsanitize=address,undefined -fno-sanitize-recover=all" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-san -j"$(nproc)"
ASAN_OPTIONS=detect_leaks=1 sh data/sample/peec_check.sh "$PWD/bin/peec" /tmp/peec-san
```

**変更したら必ず `peec_check.sh` を通すこと。** 50 件すべて OK でなければ
完了ではない。「実行が通る」だけでは不十分で、この検証群が物理の番人になっている。

計測時の注意 : ソルバーが失敗すると `zin.csv` / `peec.log` は**前回の実行結果が
残る**。自前で値を読むときは先に `grep "normal end" peec.log` を確認すること
(`peec_check.sh` の `run()` は確認している)。これを怠って誤った結論を出した実例がある。

## ソース構成

| ファイル | 役割 |
|---|---|
| `src/input_data.c` | `.peec` パーサ (`key = value`) と `peec_free()` |
| `src/wire.c` | 幾何段 : 電流セル (`p->seg`) と電荷セル (`p->chg`) を作る |
| `src/partial.c` | 幾何二重積分と部分インダクタンス (細線) |
| `src/surface.c` | 面セル (リボン) の幾何二重積分 |
| `src/volume.c` | 体積セル (矩形バー) の Hoer–Love 閉形式 (plate の厚み分割) |
| `src/potential.c` | 電位係数 P と節点容量行列 C = P⁻¹ |
| `src/skin.c` | 表皮効果 (丸線は Bessel、角線は合成式) |
| `src/mna.c` | MNA 番号付けとスタンプ |
| `src/lu.c` | 複素 LU 分解 (部分ピボット) |
| `src/iterative.c` | 掃引 LU 再利用の GMRES (acceleration = 1) |
| `src/solve.c` | 周波数掃引、Z → S 変換 |
| `src/output.c` | `peec.log` の表、`zin.csv`、Touchstone `peec.sNp`、`dist.csv` |

状態はグローバル変数ではなく `peec_t` コンテキスト構造体 1 個を main で
確保して関数に渡す。

## 移植性の絶対規則 (Windows/MSVC で実際に踏んだもの)

Linux/macOS では通るが Windows で落ちるものだけを挙げる。

1. **C99 VLA 禁止** (MSVC C2057/C2466)。`malloc` + `a[i * n + j]` のフラット配列。
2. **OpenMP for のインデックスは事前宣言する** (MSVC C3015)。MSVC の OpenMP は
   2.0 相当で `#pragma omp parallel for` の直後に `for (int i = ...)` と書けない。
   `int i;` を前置して `for (i = ...)` にする。
3. **C99 `<complex.h>` は使わない** (MSVC 非準拠)。`include/complex.h` の
   `d_complex_t` と `d_add` / `d_mul` / `d_div` / `d_inv` を使う。
4. **float\* / double\* の取り違え禁止**。配列の実型と読み出しポインタ型の不一致は
   Windows で 0xC0000005 クラッシュする (glibc は偶然耐えるので Linux では出ない)。
5. libm リンクは CMake の `MATH_LIB` 変数経由 (Windows に m.lib は無い)。
6. MSVC フラグは CMakeLists の既存ブロックに従う
   (`/utf-8`, `_USE_MATH_DEFINES`, `_CRT_SECURE_NO_WARNINGS`)。
7. 数学・物理定数は `peec.h` の自前マクロ (`PI` / `EPS0` / `MU0` / `C0`) を使う。
8. **外部ライブラリを追加しない**。OpenMP のみ任意依存で `#ifdef _OPENMP` でガードする。

`.claude/hooks/check-portability.sh` が (1)(2)(3) を機械的に検査する。
Codex から使う場合も `sh .claude/hooks/check-portability.sh` を直接叩ける。

## 壊してはいけない不変条件

どれも「壊れても実行は通るが結果が静かに狂う」種類のもの。括弧内が番人の判定名。

1. **幾何カーネルの統一** — `neumann_pair()` / `neumann_self()` は同一の細線縮約
   カーネル `R = sqrt(|dr|² + a²)` を使う。これにより `Σᵢⱼ Iᵢⱼ = I_self(全長)` が
   厳密に成立し、合計インダクタンスが分割数に依存しなくなる。片方の項だけ
   正則化を変えると両方が静かに劣化する。
   (`wire L (ndiv=1)` / `wire L (ndiv=8)` / `wire Ctotal`)
2. **静的部と遅延補正部で評価法を揃える** — 放射抵抗は L 項と P 項の大きな値
   どうしの差として現れる。片方だけ中点近似に切り替えるとその差が汚染される。
   (`short dipole R_rad` / `dipole Rin at res.`)
3. **容量ありのとき実ノードを接地しない** — 容量の電位基準は無限遠。実ノードを
   接地すると電荷が無限遠へ逃げる経路ができ、遅延ありで放射抵抗の打ち消し
   (`Σq = 0`) が壊れて **Rin が負になる** (実際に踏んだ: 微小ダイポールで −23 Ω)。
   `mna_numbering()` は `capacitance = 1` のとき基準を必ずノード 0 にする。
   (`short dipole R_rad`)
4. **面セルの幅による規格化** — 面導体セルは `Î = (1/(w₁w₂))∬∬dS dS′/R` を返す。
   w→0 で細線の式に一致するので `Lp` と `P` の式を細線と共有できる。
   (`strip L (surface)`)
5. **電荷セルは幾何段で明示的に作る** — 面導体で「電流セルの半分」から電荷セルを
   導出すると面内 2 方向ぶん二重計上になる。(`plate C (extrapolated)`)
6. **体積カーネルのゼロ引数ガードと加法性** — Hoer–Love の原始関数 `hl_F(x,y,z)`
   (volume.c) では、引数が 0 のとき**係数が消える項だけ**を 0 にする。`F(0,y,z)` は
   x を因子に持つ項が消えるだけで残りは有限に残る。「引数に 0 があるから」と
   F 全体を 0 にすると自己項が約 38 倍ずれ、重なり配置への移行が不連続になる
   (実際に踏んだ)。対数項は `asinh` 形で書く (`log((x+R)/√(y²+z²))` は x < 0 で
   桁落ち)。正しく実装すると体積分割の部分和が全体に一致する加法性が機械精度で
   成り立つ。また体積セルの厚み方向再配分は Lp 行列が陽に解くので、`zint_seg()`
   は体積セルに表皮効果の合成式を適用しない (二重計上)。
   (`thick plate L (Grover)` / `thick plate R (ndivt=1/2)` / `slab skin ratio (extr.)`)

## メモリ

**`peec_t` に動的配列メンバを足したら `peec_free()` (src/input_data.c) にも足す。**
忘れると CI の `sanitize` ジョブが落ちる。密行列 (`lp` / `cmat`) が規模に比例して
増えるので実用上も効く。

## 並列化とスレッド数不変性

OpenMP で並列化しているのは `lp_fill` (partial.c)、`pot_fill` (potential.c)、
`lu_decomp` の残余行列更新 (lu.c)、GMRES の行列ベクトル積 (iterative.c) の
4 箇所。**いずれも要素・行ごとに独立で、順序依存の加算 (リダクション) を
持たない** (GMRES の内積・Gram-Schmidt は直列)。スレッド数を変えても結果は
ビット単位で一致する。

- 並列ループ内で共有配列に `+=` しない (`pot_fill` は一時配列に出してから直列で集約)。
- `peec_check.sh` が `-n 1` と `-n 4` の `zin.csv` 完全一致を判定している。
  リダクションを持つ並列化を足すとここが落ちる。その場合は「一致する」という
  README の主張ごと見直すこと。

## 入力キーを足すとき

- パースは `src/input_data.c` の `else if` 連鎖の末尾に足す。
- **既定値は「キー省略時に従来動作と完全一致」**にする (後方互換)。
  `memset(p, 0, sizeof(peec_t))` の後に既定値を代入する。
- 未知のキーは無視する (前方互換)。末尾に足すだけでこの性質は保たれる。
- README の入力キー表に 1 行足す。

## 検証ケースを足すとき (新機能には必須)

1. `data/sample/<name>.peec` を作り、**先頭コメントに解析解の導出を書く** —
   使った公式、代入した数値、期待値、許容誤差の根拠。後から検算できることが目的。
2. `data/sample/peec_check.sh` に判定を足す。`chk <label> <actual> <expected> <tol>`
   を使い、値は `zin.csv` か `peec.log` から awk で取る。
   ヘルパ : `getL` (= Xin/2πf)、`getR` (= Rin)、`getC` (= 総容量)、`getS` (= S 成分)。
3. **期待値はコードとは独立な出所にする**。教科書の公式、文献値、または別経路の
   実装との相互検証。コード自身の出力を期待値にすると回帰テストにしかならない。
4. スクリプトは POSIX sh + awk/grep/sed のみ。CI は Windows でも同じスクリプトを
   Git Bash (`shell: bash`) で走らせる。

許容誤差の目安 :

| 種類 | 目安 |
|---|---|
| 解析式を厳密に再現するはずのもの | 0.1% |
| 近似式・級数展開との比較 | 0.5% |
| 漸近極限との比較 (補正項が残る) | 1% |
| 低次基底の離散化誤差を含むもの | 2〜10% |

許容を緩めるときは、なぜ緩いのかを .peec のコメントに書く。

## 多ポート解析

`port` を複数書くと Z 行列 (ポート j に 1A 注入、他ポート開放) と S 行列
(電力波の定義、`solve.c` の `z_to_s()`) を求め、`peec.log` の表と Touchstone
`peec.sNp` に出力する。基準抵抗はポートごとに指定できるが、Touchstone 1.1 は
1 個しか記録できないので port#1 の値を書いて注記する。
**2 ポートだけ Touchstone の列順が S11 S21 S12 S22 と転置になる** (仕様)。

## CI

`.github/workflows/ci.yml` の 4 ジョブ :

| ジョブ | 内容 |
|---|---|
| `build-linux` / `build-macos` / `build-windows` | ビルド + `peec_check.sh` (全 OS 同一スクリプト) |
| `sanitize` (Linux) | ASan + UBSan + LeakSanitizer で同じ検証を実行 |

タグ `v*` push で Release にバイナリを添付する。

## 既知の制限

- 厚み方向の電流分布 (体積セル) は `plate` の `分割数t` で対応済み。ただし
  Hoer–Love 閉形式が使えるのは軸・断面がそろった平行バー対のみで、それ以外
  (遠方・非平行) は中立面リボンの面積分に落ちる。`wire` / `bar` の断面分割は
  未対応。層間は各格子点で並列接続 (共有ノード、VFI 近似) であり、厚み方向の
  伝導は陽には解かない。
- 曲面・非矩形面は未対応 (矩形/折れ線の集合で近似)。
- 電流は区間ごと一定、電荷はセルごと一定の低次基底。共振近傍の精度は分割数依存。
- 行列は密。O(N²+M²) の積分と O((N+M)³) の LU (LU は並列化済みだがオーダーは不変)。
  `acceleration = 1` で掃引の LU 回数を減らせる (前処理 GMRES、密 LU と実質同一の
  結果、収束しなければ自動で LU に戻る)。fill の O(N²) と行列の密メモリは不変
  (ACA 圧縮が次の段階)。
