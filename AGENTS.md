# AGENTS.md — OpenPEEC

準静的〜フルウェーブ PEEC (部分要素等価回路法) 電気回路ソルバー (C11)。
[OpenFDTD](https://github.com/Sirokujira/OpenFDTD) の姉妹プロジェクトで、
ビルド規約・移植性規則を共有する。**外部数値ライブラリに依存しない**
(C11 + CMake、OpenMP のみ任意)。

集中定数 MNA + 導体形状からの部分要素抽出 (インダクタンス L / 電位係数 P /
抵抗 R) → 入力インピーダンス Zin(f)、Z / S 行列、電流・電荷分布、遠方界。
集中定数には伝送線路 (`tline`、TEM 線路を進行波の関係で組む) も含む。
導体は丸線 (`wire`) / 角線 (`bar`) / 面導体 (`plate` = 矩形、`quad` = 凸四辺形、
`disk` = 円板)、らせん / 円形ループ (`helix`、直線 `wire` 列に展開)。`plate` は `分割数t` (省略時 1) を 2 以上にすると厚み方向にも
分割し、体積セル (矩形バー) として厚み方向の電流分布 (表皮効果) を解く。
ほかに無限 PEC 地板 (`groundplane`、鏡像法)、誘電体ブリック (`dielectric`、
Ruehli の過剰容量)、遠方界後処理 (`farfield` → `far.csv`、平面波併用時は散乱界・RCS `farpw.csv`)、近傍界後処理
(`nearfield` → `near.csv`、平面波併用時は全界 `nearpw.csv`)、平面波入射
(`planewave` → `pw.csv`、EMC イミュニティ)、過渡応答 (`transient` →
`tran.csv`、掃引の逆フーリエ変換)、対数掃引 (`frequency ... log`)。
いずれもキー省略時は無効 (従来動作と完全一致)。

CSV → HDF5 の変換は `tools/peec2h5.py` (numpy + h5py)。**ソルバー本体は
外部ライブラリに依存しない**規約 (下記 8) を守るため、HDF5 はスクリプト側に
置いてある。ビルド・CI には一切影響しない。

> このファイルは Claude Code 用の `CLAUDE.md` / `.claude/rules/*.md` と同じ内容を
> 単独で読めるようまとめたもの。**片方だけ直すと食い違うので、規約を変えたら
> 両方直すこと。**

## ビルドとテスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 検証 (解析解・文献値との比較、265 判定)
sh data/sample/peec_check.sh "$PWD/bin/peec" /tmp/peec-check

# メモリ健全性 (CI の sanitize ジョブと同じ)
cmake -B build-san -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_FLAGS="-g -fsanitize=address,undefined -fno-sanitize-recover=all" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-san -j"$(nproc)"
ASAN_OPTIONS=detect_leaks=1 sh data/sample/peec_check.sh "$PWD/bin/peec" /tmp/peec-san
```

**変更したら必ず `peec_check.sh` を通すこと。** 265 判定すべて OK (h5py が
無い環境の `peec2h5 converter` だけは SKIP) でなければ完了ではない。
「実行が通る」だけでは不十分で、この検証群が物理の番人になっている。

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
| `src/polygon.c` | 多角形セル (quad/disk パネル) の幾何二重積分 (多角形閉形式) |
| `src/potential.c` | 電位係数 P と節点容量行列 C = P⁻¹ |
| `src/basis.c` | 区分線形の電流基底 (basis = 1) の結合 K / X、平面波の右辺と遠方界の傾斜関数の項 |
| `src/spice.c` | 準静的な PEEC 等価回路の SPICE ネットリスト出力 (spice = 1 → peec.cir) |
| `src/skin.c` | 表皮効果 (丸線は Bessel、角線は合成式) |
| `src/mna.c` | MNA 番号付けとスタンプ |
| `src/lu.c` | 複素 LU 分解 (部分ピボット) |
| `src/sparse.c` | 疎行列 LU (最小次数順序 + Gilbert-Peierls、前処理の Schur 補元用) |
| `src/iterative.c` | GMRES (acceleration = 1 の掃引 LU 再利用と compression = 1 の行列フリー) |
| `src/hmatrix.c` | Lp / P の H 行列圧縮 (クラスタツリー + ACA、compression = 1) |
| `src/precond.c` | 葉ブロック消去 + 回路 Schur 補元の前処理 (compression = 1) |
| `src/solve.c` | 周波数掃引、Z → S 変換 |
| `src/output.c` | `peec.log` の表、`zin.csv`、Touchstone `peec.sNp` / 2.0 `peec.ts`、混合モード `mixed.csv`、`dist.csv` |
| `src/farfield.c` | 遠方界後処理 (`farfield` → `far.csv`、D / G / 放射効率。平面波の散乱界・RCS・光学定理 → `farpw.csv`) |
| `src/nearfield.c` | 近傍界後処理 (`nearfield` → `near.csv` / 平面波の全界 `nearpw.csv`、E / H の遅延ポテンシャル) |
| `src/transient.c` | 過渡応答 (`transient` → `tran.csv`、掃引の逆フーリエ変換) |
| `tools/peec2h5.py` | CSV → HDF5 変換 (本体の依存を増やさないための外付け) |
| `tools/basis_study.py` | 電流基底の次数の検討 (純 Python の独立実装。basis = 1 の検証の参照値) |

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
   これは**ソルバー本体 (`src/`, `include/`, `CMakeLists.txt`) の規則**。
   `tools/` の後処理スクリプトは対象外で、そこでなら numpy / h5py 等を使ってよい
   (`tools/peec2h5.py` の HDF5 出力がこれ)。ビルドと CI に影響しないことが条件。

GCC の最適化で踏んだもの (Linux CI) :

9. **FMA 縮約は切る (`-ffp-contract=off`、CMakeLists の GNU ブロック)**。
   ubuntu-26.04 ランナーの GCC 15.2 は既定で `-march=x86-64-v3` (AVX2 + FMA)
   になり、`-O2` 以上の SLP ベクトル化が farfield.c の `ff_nvec()` を誤って
   コンパイルした (far.csv / farpw.csv が丸め誤差ではなく狂う。Debug の
   sanitize・macOS・Windows・ローカル GCC 13/14・clang はすべて通る)。
   `-ffp-contract=off` / `-fno-tree-slp-vectorize` / FMA なしの `-march` の
   どれでも直る。縮約を切ると結果が FMA の有無に依存しなくなり、FMA の無い
   既定の x86-64 では出力がバイト一致のまま。Linux だけ落ちたら既定の -march
   (`gcc -Q --help=target | grep march`) を確かめ、CI 上でファイル単位の `-O0`
   とフラグ単位の切り分けを回す。

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
   例外は地板 (`groundplane`) があるとき : 鏡像電荷 −q が総電荷の保存を自動で
   満たすため、ノード 0 (= 地板電位) にポートを繋いでよい (モノポール給電)。
   (`short dipole R_rad` / `monopole Rin at res.`)
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
   bar の断面分割 (`分割数w 分割数t`) も同じ体積セルで、区間の両端ノードを共有する。
   丸線の断面分割 (`wire` の `断面分割数` n) も同じで、n × n 格子のうち中心が円内の
   セルを残し、総面積が πa² になるよう格子間隔を伸縮する (DC 抵抗が厳密)。
   Hoer–Love の交代和は条件数 ~ (区間長/セル幅)³ なので、区間長はセル幅の数百倍
   以下にする (座標の丸めの差 1e-16 が L に 1e-7 で出る)。
   (`thick plate L (Grover)` / `thick plate R (ndivt=1/2)` / `slab skin ratio (extr.)` /
   `bar sub L (Grover)` / `bar sub L (10x4 = 5x2)` / `bar slab = plate slab` /
   `round sub skin c_iso` / `round sub prox c_prox` / `round sub skin x=2` /
   `round sub L internal` / `round sub image = pair`)
7. **多角形セルの規格化と桁落ちフリーの対数項** — パネル (quad/disk) のセルは
   `w̄ = 面積/len` で規格化する ((4) の一般化。`len·w̄ = 面積` なので P の式は
   共通のまま)。`poly_potential()` の対数項は l < 0 の側を恒等式
   `R + l = R0²/(R − l)` で計算する。素朴な和は隣接セルの求積点が共有辺の
   延長線近傍に乗ると 0 になり、log が inf → MNA で NaN になる (実際に踏んだ)。
   DC 抵抗は構造四辺形格子の双対幅から取る (重なりのない三角形辺基底は一様流の
   散逸を約 2 倍に過大評価するので使わない)。断面が plate と同じ幾何量
   (双対幅 × 厚さ) なので、**表皮効果も plate と同じ合成式を適用する** —
   合成式の DC 極限は `s->res` に厳密に一致するので不連続は生じない。適用
   しないでいると厚み方向の表皮効果が丸ごと落ちる (1 GHz・厚さ 0.1 mm の
   銅シートで −24% だった)。
   (`quad sheet Rin (DC)` / `quad sheet L vs plate` / `quad sheet R vs plate (HF)` /
   `quad sheet R (skin at DC)` / `disk C (extrapolated)` / `annulus Rin (DC)`)
8. **鏡像の符号規約 (groundplane)** — 鏡像セルは幾何鏡像 (`seg_image()`,
   partial.c)。Lp は `(t₁·t₂′) I(s₁, M s₂)` を**減算** (= 鏡像電流の
   面に平行な成分反転・垂直成分保存)、P は鏡像電荷 −q の相互項を減算、
   **自己項 (対角) も自分の鏡像との相互項を引く**。遠方界・近傍界・平面波の
   反射波も同じ規約 (鏡像位置 + t′ + 電流 −I)。M は等長変換かつ対合なので
   行列の対称性は保たれる。直交する複数の面 (3 枚まで) では面の部分集合
   m ごとに鏡像が 1 つで、符号は (−1)^|m| (`gp_sign()`。2 回反射した鏡像は
   **加算**)。平行な 2 面は無限級数になるのでパーサが拒否する。
   (`wire-gp L (images)` / `wire-gp C (to plane)` / `monopole = dipole/2` —
   モノポール = ダイポール/2 が全符号を同時に判定する。複数面は
   `corner3 L (7 images)` / `corner3 C (2 cells)` / `x-plane dipole far` /
   `monopole x-plane = z` / コーナーの `recip`)
9. **誘電体は εr* − 1 だけを枝に載せる (二重計上禁止)** — 真空 (ε0) は電位
   係数 P が受け持ち、誘電体の枝は複素比誘電率 εr* = εr(1 − j tanδ) の
   **過剰分だけ**を持つ : `Y = jωε0(εr*−1)A/len = ω(gexc + j cexc)`、
   枝は Z = 1/Y (mna.c が `seg->diel` を見る。skin は適用しない)。
   εr = 1 のブリックは**セルを作らない** (作ると 1/(jω·0) で NaN。
   スキップすれば真空と bit 一致)。tanδ = 0 なら Z = −j/(ωC_e) で従来と一致。
   周波数分散 (多極 Debye / Lorentz、npole > 0) も同じ原則で、`zint_diel()`
   (mna.c) が周波数ごとの εr*(f) から過剰分だけを組む。定数 tanδ と分散の
   併用は損失の二重計上なので入力検証で拒否する。Lorentz 極
   Δε f0²/(f0²−f²+jfg) は分母の虚部が +jfg (e^{jωt} 規約で因果的・損失側)。
   導電性の Drude 項 −fp²/(f²−jfg) と導電率 −jσ/(ωε0) は直流で発散するので
   直流極限の和 ε∞ + ΣΔε = εr に入れない (Y は直流で σA/len に収まる)。
   導電性があれば εr = 1 でもセルを作る。
   (`dielectric dC (pp)` / `dielectric epsr=1 noop` / `dielectric G (tand)` /
   `dielectric tand=0 noop` / `debye dY (3 freqs)` / `lorentz dY (5 freqs)` /
   `lorentz+debye dY` / `drude dY (5 freqs)` / `sigma dY (5 freqs)` /
   `sigma brick (epsr=1) G`)
10. **遠方界の規格化とポインティング整合** — r E = −j(ωμ0/4π)[N − (N·r̂)r̂]、
    U = |rE|²/(2η0)。係数を触るとパターン (D) は変わらず**効率だけが静かに
    狂う**ため、Prad = ∮U dΩ が Pin = Re(Zin)/2 と一致することが番人になる。
    (`dipole ff efficiency` / `short dipole D` / `monopole D = 2 x dipole`)
11. **送信と受信は相反定理で結ばれている** — 平面波入射 (`mna_rhs_planewave`)
    と遠方界 (`farfield.c`) は独立な実装だが、離散化した系でも
    `l_eff = rE·2λ/(−jη0)`、`Voc = −E0(ê·l_eff)` が厳密に成り立つ。符号 (−) は
    ポート規約 (n1 に +1A 注入 vs Voc = v(n1)−v(n2)) によるもので、**角度・
    偏波によらず厳密に −1 倍**になる (バグなら一定倍率にならない)。両者は
    位相規約・sinc 積分・偏波ベクトル・鏡像規約を共有しているので、片方だけ
    触ると必ず落ちる。(`recip ...` 6 通り / `dipole pw leff` / `loop pw |Voc|`)
12. **過渡応答は exp(+jωt) 規約と等間隔掃引に依存する** — `transient.c` の
    合成はコードベース全体の exp(+jωt) 規約 (inductor が +jωL) と一致して
    いる必要がある。**純抵抗の判定は S11 が実数なので共役の取り違えを検出
    できない**ため、周波数依存のあるケースが要る (微小ループは Voc ∝ jω な
    ので時間波形が励振の微分になり、共役を誤ると符号が反転する)。掃引は
    f_k = k·df でなければ時間軸が定義できず、`tran_check()` が検査して
    対数掃引を拒否する。DC 項は実数条件のもとで下 2 点から線形補外。
    (`loop transient (d/dt)` / `tdr R=...` 4 通り / `transient sweep guard`)
13. **非一様格子の幾何は格子座標の隣接中点 (双対区間) から導く** — 縁寄せ
    格子 (`grading = 1`) の plate はセル幅・中心・電荷セルを格子座標配列の
    隣接中点 [lo, hi] から導く (Σwid = 全幅が厳密で DC 抵抗の厳密性が保た
    れる)。quad/disk は `panel_node()` のパラメータ変換だけで全幾何が整合。
    `node =` で格子点に繋ぐ入力は縁寄せ後の座標に置くこと — 座標が食い違うと
    nodetol マージが外れて**接続が静かに切れる** (quad 1 スクエアで +46%、
    実際に踏んだ)。(`graded quad Rin (DC)` / `graded quad L vs plate` /
    `graded thick plate R` / `graded8 beats uniform16`)
14. **圧縮経路 (`compression = 1`) は密経路と同じ要素評価を共有する** —
    Lp は `lp_entry()`、P は `pot_entry()` の 1 箇所だけで評価し、求積の
    次数が引数順に依存する (不変条件 13) ため添字を上三角に正規化する。
    H 行列の matvec は**葉**単位で並列化する ([0,n) を重複なく分割する)。
    木の階層をまたぐ行区間で並列化すると複数スレッドが y の同じ要素を
    read-modify-write して静かに寄与を落とす (実際に踏んだ)。
    (`compress wire L` / `compress plate Rin` / `compress thread invar.`)
15. **圧縮の容量性 PEEC は電荷を陽な未知数にする** — 節点容量 C = P⁻¹ は
    作らず、「KCL に +jωq」「電荷行に Pq − v = 0」と組む。電荷行から
    q = P⁻¹v なので密経路と厳密に同じ系。符号を取り違えると Xin が反転する。
    近似が入らない規模で両経路の Zin がビット一致することが最強の判定。
    (`compress cap Rin` / `compress cap L` / `compress cap exact`)
16. **前処理の Schur 補元は疎に持ち、厳密に分解する** — 疎にするのは格納と
    分解だけで前処理そのものは変えない (近似 ILU ではない) ので、GMRES の
    反復数も収束先も変わらない。最小次数の消去グラフを作るとき **A + A^T の
    重複辺を必ず潰す** — (i,j) と (j,i) が両方あると次数が 2 倍に膨らみ、
    次数バケット `head[deg]` (長さ n+1) を突き抜けて隣の確保領域を壊す
    (実際に踏んだ : n=18 で `adj[0].len = 34`、ASan では再現せず)。
    MNA には対角が構造的にゼロの行 (電圧源の枝電流) があるので、閾値
    ピボットは非対角へ落ちられること。
    (`compress schur Rin` / `compress schur L` / `compress cap exact`)
17. **近傍界は解と同じ規約・同じ求積点で評価する** — 位相規約は
    exp(+jωt) / g = exp(-jkR)/R で farfield.c と共通、求積点は
    `cell_qpts()` (polygon.c) の 1 箇所から取る (セル形状の振り分けを二重に
    持たない)。電流セルの幅規格化は不変条件 4、電荷のサブセル重みは
    `pot_entry()` と同じ len/carea、鏡像規約は不変条件 8。`retardation = 0`
    では k = 0 (部分要素が静的なのに場だけ遅延させない)。地板面上で接線 E と
    法線 H が機械精度で 0 になることが鏡像符号の番人。平面波の全界
    (`nearpw.csv`) の入射界は右辺と同じ `pw_waves()` (mna.c) から組み、
    retardation に依らず k = ω/c。平面波の誘起電流・電荷は補償定理で
    ポート経路と突き合わせる。
    (`near loop Hz` / `near disk Ez` / `near vs far` / `near gp boundary` /
    `pw standing` / `loop shield` / `pw compensation` / `nearpw gp boundary`)
18. **混合モード S は単端 S の直交な組替えで作る** — `S_mm = M S Mᵀ`、M は
    差動対ごとに a_d = (a_p − a_n)/√2、a_c = (a_p + a_n)/√2 (対に入らない
    ポートは 1)。√2 を落とすと |S| が 2 倍、符号を落とすと d と c が入れ替わる。
    基準は差動 2·z0、同相 z0/2 で、対の z0 が等しいときだけ成り立つので
    入力検査で強制する。並びは d1..dP, c1..cP, 単端 (番号順)。
    (`mixed-mode (2 pairs)` / `mixed-mode Zref (2Z0,Z0/2)` /
    `twinline no conversion` / `twinline norm preserved`)
19. **散乱遠方界は送信と同じ放射ベクトルで作り、光学定理は回路と恒等的** —
    `farpw.csv` の散乱界は `segipw` (誘起電流) を `ff_nvec()` (far.csv と共通)
    に通すだけ。消光断面積 σ_e = −(4π/k) Im[ê·rE_s(d)/E₀] (exp(+jωt) 規約。
    e^{−iωt} の教科書式の複素共役で符号が逆) は、離散系では回路が受け取る電力
    ½ Re Σ emf conj(I) と恒等的に等しい。したがってエネルギー保存
    σ_e = σ_s + σ_a は誘起電流に放射の反作用がある retardation = 1 でだけ
    成り立つ (retardation = 0 では σ_e は吸収分のみ)。
    (`loop RCS sigma(y^)` / `loop sigma_ext (optical)` / `dipole scat energy` /
    `farpw compensation`)
20. **伝送線路は進行波の関係で組み、T = a e^{−jω td} の偏角が遅延** —
    端子電流 i1 / i2 (どちらも線路へ流れ込む向き) を未知数に足し、
    `v1 − Z0 i1 = T (v2 + Z0 i2)`、`v2 − Z0 i2 = T (v1 + Z0 i1)` の 2 行を
    スタンプする (mna.c)。Z / Y 行列の形は cot / csc を含み ω td = nπ で
    特異になるが、この形は対角が常に −Z0。exp(+jωt) 規約で e^{−jωtd} が遅延
    なので、符号を取り違えると反射が前進し S21 の位相が逆に回る。端子対ごとに
    電流は行き帰りで等しい (同相の経路を作らない) ので、容量性 PEEC の総電荷の
    保存 (3) を壊さず浮いた端子対で構造に繋げる。損失 dB は無ひずみ線路の
    片道減衰で、Z0 は実数のまま。`rlgc` の形 (R + Rs√f、G + ωC tanδ) は複素 Z0 と
    T = e^{−γ·長さ} を同じスタンプに入れるだけ (主値の平方根 Re ≥ 0 で受動側)。
    未知数ブロックは offS より前 (回路側) なので圧縮経路では前処理の Schur 補元に入る。
    (`tline Zin (lossless)` / `tline Zin (2 dB loss)` / `tline matched S (3 dB)` /
    `tline TDR delay` / `tline-fed dipole` / `tline compress = dense` /
    `tline rlgc Zin` / `tline rlgc vs ladder` / `tline rlgc = Z0/td`)
21. **励振の重みはポート解の重ね合わせ、Pin は wᴴZw** — `excitation` の重み
    w_j で電流 (`segi`) と**電荷 (`cellq`) の両方**を重ね合わせる (output.c の
    `exc_mix()`。近傍界の −∇φ は電荷から来る)。入力電力は相互結合を含む
    Pin = ½Re(Σ conj(w_i) Z_ij w_j) (`exc_pin()`)。行が無ければ base をそのまま
    返して従来とビット一致。重みの位相は e^{+j 位相} (電源と同じ規約)。
    (`array antiphase null` / `array near Ez (x=0)` / `array Pin = Re(wHZw)/2` /
    `array ff eff (antiphase)` / `array ff eff (in-phase)` /
    `excitation superposition` / `excitation default`)
22. **区分線形の電流基底 (basis = 1) は電荷で駆動する傾斜関数** — J = Σ I_k p_k +
    Σ c_n χ_n、c_n = jω q_n、χ_n は節点 n の電荷セルの半区間上の (L_h − s) ô_h/L_n
    (ô_h は節点から外向き、節点側の端は seg_t の `nend`)。χ_n で試験すると
    <χ_n, −∇φ> = φ(節点) − φ̄_n なので節点電位は点電位で、電荷行に
    jω(Kᵀ I + X c)、KVL 行に + ω² K q を足す (電荷を陽な未知数にする並びが必須)。
    平面波の右辺 (電荷行に + <χ_n, E_inc>) と遠方界 (c_n χ_n の放射) も同じ関数で
    足さないと送受の相反が崩れる。ω → 0 で従来に戻り、basis = 0 はビット一致。
    (`basis1 = study (4/arm)` / `basis1 coarse accuracy` / `basis1 ff efficiency` /
    `basis1 monopole = dip/2` / `basis1 2-port symmetry` / `basis1 low-f = basis0` /
    basis 系の `recip`)
23. **SPICE 出力 (spice = 1) は f = 0 の部分要素を solve() の前に書く** — 遅延ありでは
    solve() が Lp / C を周波数ごとに作り直すので、main の静的な充填の直後に書く。
    セルは n1 —R— X_k —L— n2 (誘電体は R の代わりに C_e)、L の dot は n1 側、
    K = Lp_kl/√(Lp_kk Lp_ll) は逆向きで負。容量は Maxwell の行列を節点 - 基準
    (行和) と節点間 (−C_ij) に分け、負の枝も省かない (`capacitor` は負を受け付ける)。
    数値は %.16e。読み戻しは元の PEEC と丸めの差だけで一致する。
    (`spice reimport (mixed)` / `spice reimport (gp 4-port)` / `spice reimport (diel)` /
    `spice netlist integrity` / `spice kmin threshold`)

## メモリ

**`peec_t` に動的配列メンバを足したら `peec_free()` (src/input_data.c) にも足す。**
忘れると CI の `sanitize` ジョブが落ちる。密行列 (`lp` / `cmat`) が規模に比例して
増えるので実用上も効く。

## 並列化とスレッド数不変性

OpenMP で並列化しているのは `lp_fill` (partial.c)、`pot_fill` (potential.c)、
`lu_decomp` の残余行列更新 (lu.c)、GMRES の行列ベクトル積 (iterative.c)、
H 行列のブロック充填と matvec (hmatrix.c)、前処理の葉 LU (precond.c)、
近傍界の観測点ループ (nearfield.c) の 7 箇所。**いずれも要素・行・ブロックごとに独立で、順序依存の加算
(リダクション) を持たない** (GMRES の内積・Gram-Schmidt は直列。H 行列の
matvec は葉行クラスタ = 出力の互いに素な区間ごとに並列化する。ブロックの
行区間は木の階層をまたいで重なるので、行クラスタ単位の並列化は競合する —
実際に踏んだ)。スレッド数を変えても結果はビット単位で一致する。

- 並列ループ内で共有配列に `+=` しない (`pot_fill` は一時配列に出してから直列で集約)。
- `peec_check.sh` が `-n 1` と `-n 4` の `zin.csv` / `near.csv` / `nearpw.csv` 完全一致を
  判定している。
  リダクションを持つ並列化を足すとここが落ちる。その場合は「一致する」という
  README の主張ごと見直すこと。
- 疎行列 LU (sparse.c、前処理の Schur 補元) は**全体が直列**。GMRES の前処理
  として固定の順序で呼ばれるだけなので、スレッド数不変性の心配が無い。

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
5. **awk は実装差に注意 (CI は gawk、ローカルは mawk のことがある)**。gawk は
   0 除算で停止し、mawk は inf / NaN を黙って返す。NaN は `if (e > mx) mx = e`
   の最大値から**黙って抜け落ちる**ので、その点の比較が無かったことになる
   (実際に踏んだ : tan の発散する θ = π/2 の行がローカルでは比較されず、CI では
   停止した)。期待値の式は発散点の無い形 (tan より sin/cos、tanh より
   cosh/sinh) で書き、誤差の最大値は `if (!(e <= mx)) mx = e` と書いて NaN を
   判定の失敗に変える。点数 (`n == 5` など) も必ず確かめる。
6. **走っている `peec_check.sh` を編集しない**。sh はスクリプトを逐次読みするので、
   実行中 (sanitize ビルドだと 1 時間近くかかる) に行を足すと読み位置がずれ、
   無関係な行で構文エラーや不可解なコマンド (実際に踏んだ : `cp: invalid option
   -- '8'`) になって止まる。長い検証は写しを `data/sample/` に置いて走らせる
   (`SRC` はスクリプトの場所から決まる) か、編集を終えてから走らせる。

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
1 個しか記録できないので port#1 の値を書いて注記する。`touchstone = 2` で
Touchstone 2.0 の `peec.ts` も書き、全ポートの基準抵抗を `[Reference]` に記録する
(`output_touchstone2()`、データは行優先で 2 ポートは `[Two-Port Data Order] 12_21`)。
**2 ポートだけ Touchstone の列順が S11 S21 S12 S22 と転置になる** (仕様)。

`mixedmode = p1 n1 [p2 n2 ...]` で混合モード S_mm = M S Mᵀ (差動 / 同相) を
`mixed.csv` と `peec.log` に出す (`output_mixed()`)。対の 2 ポートの z0 は
等しいこと (差動 2·z0 / 同相 z0/2 の基準はその場合だけ定義できる)。

## CI

`.github/workflows/ci.yml` の 4 ジョブ :

| ジョブ | 内容 |
|---|---|
| `build-linux` / `build-macos` / `build-windows` | ビルド + `peec_check.sh` (全 OS 同一スクリプト) |
| `sanitize` (Linux) | ASan + UBSan + LeakSanitizer で同じ検証を実行 |

タグ `v*` push で Release にバイナリを添付する。

## 既知の制限

- 厚み方向の電流分布 (体積セル) は `plate` の `分割数t` と `bar` の
  `分割数w 分割数t` (断面分割) で対応済み。ただし Hoer–Love 閉形式が使えるのは
  軸・断面がそろった平行バー対のみで、それ以外 (遠方・非平行) は中立面リボンの
  面積分に落ちる。丸線 (`wire`) は `断面分割数` で断面を等面積の階段近似に
  分割できる (収束は O(セル幅)、n = 10 で抵抗の増分 −2〜3%)。分割しない丸線の
  表皮効果は Bessel の内部インピーダンス (近接効果は扱わない)。層間・断面内は各ノードで並列接続
  (共有ノード、VFI 近似) であり、厚み方向・断面内の横方向の伝導は陽には解かない。
- 面導体は平面パネル (矩形 `plate` / 凸四辺形 `quad` / 円板 `disk`) のみ。
  3 次元曲面 (円筒など) は平面パネルのファセット化で近似する。
  `quad` / `disk` の厚み方向分割 (`分割数t`) は層を厚みのない多角形の面積分で
  結合する (plate の層は Hoer–Love) ので少ない層数での誤差が大きい。パネルの
  DC 抵抗は格子の双対幅で決まり、矩形格子では厳密、歪んだ格子・極格子
  では O(格子幅²) の系統誤差を含む (弦近似)。
- 既定 (basis = 0) は電流が区間ごと一定、電荷がセルごと一定の低次基底。共振近傍の
  精度は分割数依存 (互い違い格子で Rin は 2 次以上の収束)。`basis = 1` は細線の
  電流を連続な区分線形にする高次基底 (節点電荷で駆動する傾斜関数の Galerkin、
  basis.c) で、粗い格子で正確 (Rin 誤差が片腕 4 分割で 1/3、8 分割で 1/2) だが、
  電荷が区分一定のままなので漸近的には 1 次程度で従来より遅い (電荷も高次に
  すると節点電位 = セル平均の回路解釈が成り立たない)。面導体・体積セル・誘電体・
  圧縮とは併用不可 (近傍界は傾斜関数の電流も電荷サブセルの求積点で積分する)。独立実装の `tools/basis_study.py` で検討を再現できる。
  Xin はデルタギャップのギャップ容量で対数的に増えるので収束判定に使わない。
- 地板 (`groundplane`) は軸に垂直な無限 PEC 面で、互いに直交する 3 枚まで
  (平行な 2 面・斜めの面は不可)。積分回数は 2^面数 倍になるが未知数は
  増えない。導体は各面の + 側に置く (はみ出しはエラー)。x / y 面の遠方界は
  面の手前の半空間を φ の台形則 (面上は重み 1/2) で積分するので、φ 分割数は
  4 の倍数にして面上を格子点に載せる。
- 誘電体 (`dielectric`) は直交直方体のみ。損失は一定 tanδ または多極分散
  (Debye 極 `f de`、Lorentz 極 `lorentz f0 de g`、導電性の Drude 項
  `drude fp g` と導電率 `sigma s` の任意の組合せ) のどちらか (併用は損失の
  二重計上なので入力検証で拒否)。Drude と導電率は直流で発散するので直流極限の
  和 ε∞ + ΣΔε = εr には入れない (枝のアドミタンスは直流で σA/len に収まる)。
- 面格子の縁寄せは `grading = 1` (plate/quad 余弦・disk 正弦)。plate の
  メッシュ生成は非一様間隔対応 : 幅・中心・電荷セルは格子座標配列の
  隣接中点 (双対区間) から導く。等間隔時の式は従来と同一で省略時は
  完全後方互換。`node =` で格子点に繋ぐ入力は縁寄せ後の座標に置くこと
  (座標がずれると nodetol マージが外れて接続が切れる)。
  節点の束縛電荷は a–b 面内の双対矩形パネルで表す近似 (表面電荷が a–b 面に
  支配的な配置 — 基板・平行平板 — で正確)。`capacitance = 1` 必須。
- 遠方界 (`farfield`) の `far.csv` は port #1 の 1 A 励振、または `excitation` の
  電流励振の重ね合わせ (Pin は相互結合込みの Re(wᴴZw)/2) に対する値。効率が意味を持つのは
  `retardation = 1` のとき (準静的電流は放射を含まない)。平面波併用時の
  `farpw.csv` は誘起電流の散乱界と RCS。光学定理の消光断面積は地板なしのときのみ。
- 近傍界 (`nearfield`) の `near.csv` は port #1 の 1 A 励振 (または `excitation`)、`nearpw.csv` は
  平面波入射時の全界 (入射 + 反射 + 散乱)。入射界は `retardation` に依らず
  k = ω/c で評価する。観測点はセル寸法
  程度以上離すこと (セル上の積分は数点の求積)。`capacitance = 0` では電荷が
  無いので E は −jωA の項だけ、`retardation = 0` では部分要素と揃えて
  静的 (k = 0) に評価する。
- 平面波 (`planewave`) は 1 行が 1 波。複数行 (PW_MAX = 8 行まで) は同時入射の
  重ね合わせ (起電力・入射界・散乱界とも)。leff / RCS は 1 行目の振幅で規格化、
  光学定理の消光断面積は単一波・地板なしのときだけ。`pw.csv` はポート間に素子が
  無いときだけ Voc として読める。
- 過渡応答 (`transient`) は掃引の逆フーリエ変換なので、応答が `1/Δf` 以内に
  減衰しないと巻き込み (時間領域エイリアス) が起きる。Q の高い共振では
  分割数を増やす。時間分解の下限は Δt = 1/(2f_max) で、既定の −40 dB では
  ガウス幅が σ ≈ Δt (帯域を使い切った状態) になる。
- 面セル・多角形セルの静的部は外側求積の次数を距離で落としている
  (`ribbon_nq()` / `POLY_R7`、不変条件 13)。充填が平板で 3.1 倍・パネルで
  1.5 倍速くなり、Zin の変化は 1e-8 以下。次数を落としてよいのは静的部だけで、
  自己項・隣接項は最高次を使う。
- 行列は密。O(N²+M²) の積分と O((N+M)³) の LU (LU は並列化済みだがオーダーは不変)。
  実測では未知数 ~7500 までは充填が支配的で、それを超えると LU が追い越す。
  `acceleration = 1` で掃引の LU 回数を減らせる (前処理 GMRES、密 LU と実質同一の
  結果、収束しなければ自動で LU に戻る)。fill の O(N²) と行列の密メモリは不変
  (ACA 圧縮が次の段階)。
