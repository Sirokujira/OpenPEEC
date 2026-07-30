# OpenPEEC

A quasi-static PEEC (Partial Element Equivalent Circuit) electric-circuit
solver in C. Companion project to
[OpenFDTD](https://github.com/Sirokujira/OpenFDTD) — same build conventions,
same portability rules, no external numerical libraries.

準静的 PEEC (部分要素等価回路法) による電気回路解析ソルバー (C 言語)。
[OpenFDTD](https://github.com/Sirokujira/OpenFDTD) の姉妹プロジェクトで、
ビルド規約・移植性規則を共有します。外部数値ライブラリ (LAPACK/BLAS/HDF5)
には依存しません。

## Theory / 理論概要

- **MNA (修正節点解析)** : 集中定数素子 R / L / C / 相互インダクタンス K /
  独立電源 (AC) をスタンプした複素連立一次方程式を、周波数ごとに自前の
  複素 LU 分解 (部分ピボット) で解きます。
- **部分インダクタンス (L)** : 直線導体 (`wire` / `bar` キー) を区間分割し、
  幾何二重積分 `I = ∬ dl₁dl₂/R` から `Lp = (μ0/4π)(t̂₁·t̂₂) I` を求めて
  区間枝電流を MNA に組み込みます。ポートは 1 A 電流源励振で Zin(f) を
  計算します。
  - 自己項・平行/同一直線の相互項は解析式、一般配置は 8 点
    Gauss–Legendre の複合則 (近接時は 4×4 分割、遠方は中点近似)。
  - 細線縮約カーネル `R = √(|Δr|² + a²)` を全項で統一しているため、
    **直線導体の合計インダクタンスは分割数に依存しません**
    (`Σᵢⱼ Iᵢⱼ = I_self(全長)` が厳密に成立)。
- **導体形状** : 3 種類に対応します。
  - `wire` : 丸線 (半径 a) — 細線フィラメントとして扱う
  - `bar` : 角線 (幅 w × 厚さ t) — 等価半径で扱う。インダクタンスは
    幾何平均距離 `a_L = 0.2235(w+t)` (Grover の角線式を再現)、容量は
    薄板↔丸線の等価 `a_P = (w+t)/4`、抵抗は断面積 `R = l/(σwt)`。
  - `plate` : 平面矩形の面導体 — 格子に分割し、隣接ノード間を幅をもつ
    面セル (リボン)、各ノードの双対矩形を電荷セルとする標準的な面 PEEC。
    電流は面内 2 方向に流れます。
- **面セルの積分** : 幅で規格化した `Î = (1/w₁w₂)∬∬dS dS′/R` を使うことで、
  `Lp = (μ0/4π)(t̂₁·t̂₂)Î`, `P = Î/(4πε0 L₁L₂)` と細線とまったく同じ式に
  なります (w→0 で細線の式に一致)。内側の面積分は矩形の閉形式で厳密に
  評価するため、自己項・隣接項も安定です。
- **電位係数 (P) / 容量性 PEEC** — `capacitance = 1` で有効:
  容量セルを幾何ノードに割り当て (各区間の半分ずつを両端ノードが持つ)、
  同じ幾何二重積分から `P = I/(4πε0 L₁L₂)` を作り、節点容量行列
  `C = P⁻¹` を MNA の節点ブロックに `jωC` として加えます。
  対無限遠の総容量 `ΣᵢⱼCᵢⱼ` をログに出力します。
- **表皮効果 (R(f) + 内部 L)** — `skineffect = 1` で有効:
  丸線は内部インピーダンス `Z_int = k/(2πaσ)·I₀(ka)/I₁(ka)`
  (`k = (1+j)/δ`, `δ = √(2/ωμ0σ)`) を各区間に適用します。
  `I₀/I₁` は `|ka| < 15` で級数、以上で漸近展開。
  低周波では `R_dc + jωμ0l/(8π)`、高周波では `R = X = l/(2πaδσ)` に
  漸近します。角線は閉形式が無いため、両極限 (`R_dc` と周長 P の表皮層
  `l/(σδP)`) が厳密になる合成式で近似します。
- **遅延 (フルウェーブ PEEC)** — `retardation = 1` で有効:
  カーネルを `e^{-jkR}/R` に置き換え、L と P を周波数依存の複素行列として
  毎周波数で組み直します。特異性抽出
  `∬e^{-jkR}/R = ∬1/R + ∬(e^{-jkR}−1)/R` により、解析式を保ったまま
  第 2 項 (R→0 で −jk に収束する正則関数) だけを求積します。
  これにより放射抵抗が正しく得られます (半波長ダイポールで Rin ≈ 73 Ω)。

  遅延を使う場合、**構造ノードを接地してはいけません**。容量の電位基準は
  無限遠であり、実ノードを接地するとそこから無限遠へ電荷が逃げて構造の
  総電荷が保存せず、放射抵抗の打ち消しが壊れて Rin が負になります。
  `capacitance = 1` のときは基準ノードを常にノード 0 (= 無限遠) とし、
  ノード 0 が未使用ならどの実ノードも消去しません (容量が電位を確定させる)。

### 制限 / Limitations

- `skineffect` / `capacitance` / `retardation` は既定で無効です
  (キー省略時は従来動作と完全に一致)。
- `retardation = 0` は準静的モデルです。寸法 ≪ λ/10 で正確ですが、
  共振周波数は半波長共振より 15〜20% 高く出ます (`dipole_full.peec` 参照)。
  共振・放射を扱うなら `retardation = 1` を使ってください。
- 電流は区間ごとに一定、電荷はセルごとに一定の低次基底です。共振近傍の
  精度は分割数に依存します (半波長ダイポールで片腕 16 分割 → Rin 誤差 1% 程度)。
- 角線は等価半径による近似です。表皮効果の遷移域 (δ ~ 断面寸法) は
  両極限を繋ぐ合成式なので、その領域の精度が要る用途では丸線を使ってください。
- 面導体は平面矩形のみ (曲面・非矩形は矩形の集合で近似)。厚さ方向は
  1 層 (体積導体・厚み方向の電流分布は未対応)。表皮効果は角線と同じ
  合成式を面セルの断面 (幅 × 厚さ) に適用します。
- 導体は直線区間・平面矩形のみ (曲線は折れ線で近似)。
- 行列は密です。計算量は電流セル数 N と電荷セル数 M に対して
  O(N²+M²) の積分と O((N+M)³) の LU。遅延ありのときはこれが周波数ごとに
  必要になります。

## Build / ビルド

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 検証 (解析解と比較、許容誤差は各ケース 0.5〜2%)
sh data/sample/peec_check.sh bin/peec /tmp/peec-check

# メモリ健全性の検証 (AddressSanitizer + UndefinedBehaviorSanitizer)
cmake -B build-san -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_FLAGS="-g -fsanitize=address,undefined -fno-sanitize-recover=all" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address,undefined"
cmake --build build-san -j"$(nproc)"
ASAN_OPTIONS=detect_leaks=1 sh data/sample/peec_check.sh bin/peec /tmp/peec-san
```

- Linux / macOS / Windows (MSVC + Ninja) の 3 OS を CI で検証しています。
  あわせて Linux で ASan + UBSan (LeakSanitizer 有効) を掛けた `sanitize`
  ジョブを実行し、値の正しさとメモリ健全性を同時に判定しています。
- OpenMP は任意です (無くてもビルド・実行可能)。並列化しているのは
  **部分インダクタンス行列の充填・電位係数行列の充填・LU 分解の残余行列更新**
  の 3 箇所で、いずれも要素ごとに独立 (順序依存の加算が無い) なため
  **スレッド数によらず出力はビット単位で一致します**
  (`peec_check.sh` が `-n 1` と `-n 4` の完全一致を判定)。
  LU 支配的なケース (未知数 645、21 周波数掃引) で 4 スレッド 3.6 倍。

## Usage / 実行

```bash
peec [-n <threads>] input.peec
```

出力 :

| ファイル | 内容 |
|---|---|
| `peec.log` | 入力インピーダンス表、S パラメータ表、`=== normal end ===` |
| `zin.csv` | 各ポートの Zin (機械読み取り用) |
| `peec.sNp` | S パラメータの Touchstone 1.1 形式 (N = ポート数) |
| `dist.csv` | 電流・電荷分布 (`distribution = 1` のときのみ) |

**多ポート解析** : `port` を複数書くと、ポート j に 1A を注入し (他ポートは
開放) 全ポートの端子電圧を読むことで Z 行列を求め、電力波の定義 (Kurokawa)

```
S = F (Z − Z0) (Z + Z0)⁻¹ F⁻¹,   F = diag(1 / (2√z0ᵢ))
  → Sᵢⱼ = √(z0ⱼ / z0ᵢ) · [(Z − Z0)(Z + Z0)⁻¹]ᵢⱼ
```

で S 行列に変換します。基準抵抗 `z0` は `port` キーでポートごとに指定でき、
全ポートで等しいときは教科書どおりの `(Z − Z0)(Z + Z0)⁻¹` に一致します。
Touchstone 1.1 は基準抵抗を 1 個しか記録できないため、ポートごとに異なる
`z0` を与えた場合は port#1 の値を書き、その旨をコメント行に残します。

## Input format / 入力形式

テキスト形式。1 行目は `OpenPEEC 1 0`、最終行は `end`。`#` はコメント。
ノード id は非負整数で **0 = GND (基準)**。単位はすべて SI。

| キー | 書式 | 意味 |
|---|---|---|
| `title` | `title = 任意文字列` | タイトル |
| `resistor` | `resistor = n1 n2 R` | 抵抗 [Ω] |
| `capacitor` | `capacitor = n1 n2 C` | キャパシタ [F] |
| `inductor` | `inductor = 名前 n1 n2 L` | インダクタ [H] (名前は `mutual` 参照用) |
| `mutual` | `mutual = 名前1 名前2 k` | 結合係数 k (0 < k ≤ 1)、M = k√(L1L2) |
| `vsource` | `vsource = n1 n2 振幅 [位相deg]` | AC 電圧源 (+ 端子 = n1) |
| `isource` | `isource = n1 n2 振幅 [位相deg]` | AC 電流源 (電源内部を n1→n2、SPICE 規約) |
| `node` | `node = n x y z` | 回路ノード n を 3D 座標に束縛 (ワイヤ接続点) |
| `wire` | `wire = x1 y1 z1 x2 y2 z2 半径 導電率 分割数` | 丸線 (PEEC 抽出対象) |
| `bar` | `bar = x1 y1 z1 x2 y2 z2 幅 厚さ 導電率 分割数` | 角線 (矩形断面) |
| `plate` | `plate = ox oy oz ax ay az bx by bz 厚さ 導電率 分割数a 分割数b` | 平面矩形の面導体 (2 辺ベクトルは直交) |
| `port` | `port = n1 n2 Z0` | ポート (Zin / S パラメータの基準抵抗 Z0 [Ω])。複数書くと多ポート解析になる |
| `frequency` | `frequency = f開始 f終了 分割数` | 周波数掃引 (分割数+1 点) |
| `nodetol` | `nodetol = 1e-8` | 座標マージ許容 [m] (省略時 1e-8) |
| `gmin` | `gmin = 0` | 全ノード対地コンダクタンス [S] (省略時 0) |
| `skineffect` | `skineffect = 1` | 1 で表皮効果 + 内部インダクタンス (省略時 0 = DC 抵抗) |
| `capacitance` | `capacitance = 1` | 1 で容量性 PEEC (電位係数) を有効化 (省略時 0) |
| `retardation` | `retardation = 1` | 1 で遅延 (フルウェーブ PEEC) を有効化 (省略時 0) |
| `distribution` | `distribution = 1` | 1 で電流・電荷分布を `dist.csv` に出力 (省略時 0) |

- 導体と回路素子の接続は `node` キーで行います。導体端点 (面導体は格子点)
  の座標が `nodetol` 以内で一致するとそのノード id に接続されます。
- `plate` は原点 `o` と 2 辺ベクトル `a`, `b` で矩形 `o + s·a + t·b`
  (s,t ∈ [0,1]) を表します。格子点は `o + (i/分割数a)·a + (j/分割数b)·b`。
- ワイヤ内部の分割点は最大ノード id + 1 から記載順に自動採番されます
  (決定的)。
- 基準ノードはノード 0。`capacitance = 0` でノード 0 が使われていない場合は
  port #1 の n2 を基準にします (浮遊回路の特異回避)。`capacitance = 1` では
  常にノード 0 (= 無限遠) が基準です (上記「遅延」参照)。ログに表示されます。
- 未知のキーは無視されます (前方互換)。キー省略時の既定値は従来動作と
  一致します (後方互換)。

## Validation / 検証 (`data/sample/`)

`sh data/sample/peec_check.sh bin/peec /tmp/peec-check` が全ケースを実行し、
解析解と比較します (CI で 3 OS 自動実行)。

| ケース | 解析解 | 許容 |
|---|---|---|
| `rlc_series.peec` 直列 RLC | Zin(f0) = 50 + j0、Xin(2f0) = 47.434 Ω | 0.5% |
| `wire_single.peec` 単線 1 m / a = 1 mm | L = 1.32038 µH (ndiv = 1 と 8 で同一)、Rin = 5.4881 mΩ | 0.1% |
| `loop_square.peec` 正方ループ s = 0.2 m | L = (2μ0s/π)[ln(s/a)+a/s−0.774013] = 724.69 nH | 2% |
| `wire_skin.peec` 表皮効果 (低周波 1 kHz) | Rin = R_dc(1+q⁴/48)、内部 L → μ0l/(8π) → L = 1.37035 µH | 0.5% |
| `wire_skin_hf.peec` 表皮効果 (高周波 100 MHz) | Rin → l/(2πaδσ) = 41.52 mΩ | 1% |
| `wire_cap.peec` 単線の対無限遠容量 | C = 4πε0l²/I_self = 8.42673 pF | 0.1% |
| `twowire_cap.peec` 平行 2 線の容量 (相互項) | C = 8πε0l²/(I_self+I_par) = 10.1982 pF | 0.1% |
| `loop_square.peec` + `capacitance=1` | 電気的に小さいループでは Zin が変わらないこと | 2% |
| `dipole_full.peec` フル PEEC デモ | 直列共振 (Xin の符号反転) が存在すること | — |
| `dipole_halfwave.peec` 半波長ダイポール (遅延) | 共振長 l = 0.478λ、共振時 Rin = 73.1 Ω | 5% / 10% |
| `dipole_short.peec` 微小ダイポール (遅延) | R_rad = 20π²(l/λ)² = 0.12354 Ω | 5% |
| `loop_square.peec` + `retardation=1` | 低周波では準静的結果に一致すること | 2% |
| `bar_single.peec` 角線 1 m / 10×1 mm | L = (μ0l/2π)[ln(2l/(w+t))+0.5+0.2235(w+t)/l] = 1.14109 µH、R_dc = 1.72414 mΩ | 0.5% |
| `bar_strip_cap.peec` 薄板の容量 | 半径 w/4 の丸線と等価 → C = 9.78221 pF | 0.5% |
| `bar_skin.peec` 角線の表皮効果 | R → l/(σδP) = 0.118589 Ω (P = 2(w+t)) | 1% |
| `plate_cap.peec` 正方形平板の容量 | 8/16 分割から Richardson 補外 → 0.3667892·4πε0a = 40.810 pF | 1% |
| `tnetwork_spara.peec` 抵抗性 T 型 2 ポート | Z = [[75,50],[50,75]], Z0 = 50 → S11 = 1/21、S21 = 8/21 (実数)。相反性 S12 = S21、対称性 S11 = S22 | 0.1% |
| `tnetwork_l_spara.peec` 直列 L 入り T 型 2 ポート | ωL = 50 Ω で S11 = 0.101124+j0.561798、S21 = 0.359551−j0.224719、S22 = 0.056180+j0.089888。非対称だが相反 | 0.1% |
| `wire_dist.peec` 電流分布 (容量なし) | 電荷が溜まらないので全 8 区間に同じ 1 A (キルヒホッフの電流則、厳密) | 1e-9 |
| `wire_dist.peec` + `capacitance=1` | ポートは 1 A を入れて 1 A を出すので構造の総電荷は 0 | max\|q\| 比 1e-6 |
| `plate_cap.peec` スレッド数不変性 | `-n 1` と `-n 4` の `zin.csv` が完全一致 (未知数 225 で LU 並列経路も通る) | 完全一致 |
| `plate_strip.peec` 帯のインダクタンス | 面積分 vs GMD 近似 → (μ0l/2π)[ln(2l/w)+0.5] = 1.15966 µH | 1% |

## License

MIT
