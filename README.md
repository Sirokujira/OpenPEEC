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
- **導体断面** : 丸線 (`wire`, 半径 a) と角線 (`bar`, 幅 w × 厚さ t) に対応。
  角線は等価半径で扱います — インダクタンスは幾何平均距離
  `a_L = 0.2235(w+t)` (Grover の角線式を再現)、容量は薄板↔丸線の等価
  `a_P = (w+t)/4`、抵抗は断面積 `R = l/(σwt)`。
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
- 導体は直線区間のみ (曲線は折れ線で近似)。面・体積導体は未対応。
- 遅延ありのときは周波数ごとに O(N²) 個の積分と O(N³) の逆行列が必要です。

## Build / ビルド

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 検証 (解析解と比較、許容誤差は各ケース 0.5〜2%)
sh data/sample/peec_check.sh bin/peec /tmp/peec-check
```

- Linux / macOS / Windows (MSVC + Ninja) の 3 OS を CI で検証しています。
- OpenMP は任意です (見つかれば部分インダクタンス・電位係数行列の充填を
  並列化。無くてもビルド・実行可能)。スレッド数によらず出力は一致します。

## Usage / 実行

```bash
peec [-n <threads>] input.peec
```

出力 : `peec.log` (入力インピーダンス表 + `=== normal end ===`)、
`zin.csv` (機械読み取り用)。

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
| `port` | `port = n1 n2 Z0` | Zin 計算ポート |
| `frequency` | `frequency = f開始 f終了 分割数` | 周波数掃引 (分割数+1 点) |
| `nodetol` | `nodetol = 1e-8` | 座標マージ許容 [m] (省略時 1e-8) |
| `gmin` | `gmin = 0` | 全ノード対地コンダクタンス [S] (省略時 0) |
| `skineffect` | `skineffect = 1` | 1 で表皮効果 + 内部インダクタンス (省略時 0 = DC 抵抗) |
| `capacitance` | `capacitance = 1` | 1 で容量性 PEEC (電位係数) を有効化 (省略時 0) |
| `retardation` | `retardation = 1` | 1 で遅延 (フルウェーブ PEEC) を有効化 (省略時 0) |

- 導体と回路素子の接続は `node` キーで行います。導体端点の座標が
  `nodetol` 以内で一致するとそのノード id に接続されます。
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

## License

MIT
