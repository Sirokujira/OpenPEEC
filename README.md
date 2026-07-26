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
- **PEEC 抽出** : 直線ワイヤ (`wire` キー) を区間分割し、各区間の
  - 自己部分インダクタンス : ロッドの閉形式
    `Lp = (μ0 l/2π)[asinh(l/a) − √(1+(a/l)²) + a/l]`
  - 相互部分インダクタンス : 平行/同一直線は解析式、一般配置は
    Neumann 二重積分の 8 点 Gauss–Legendre 数値積分 (遠方は中点近似)
  - 直流抵抗 : `R = l/(σπa²)`

  を求め、区間枝電流を MNA に組み込みます。ポートは 1 A 電流源励振で
  入力インピーダンス Zin(f) を計算します。
- **v1 の制限** : 容量性結合 (電位係数) は未実装の準静的 L+R モデルです。
  電気的に小さい (寸法 ≪ λ/10) 閉電流路で有効です。表皮効果は未実装
  (抵抗は DC 値)。

## Build / ビルド

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 検証 (解析解と比較、許容誤差は各ケース 0.5〜2%)
sh data/sample/peec_check.sh bin/peec /tmp/peec-check
```

- Linux / macOS / Windows (MSVC + Ninja) の 3 OS を CI で検証しています。
- OpenMP は任意です (見つかれば部分インダクタンス行列の充填を並列化。
  無くてもビルド・実行可能)。

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
| `wire` | `wire = x1 y1 z1 x2 y2 z2 半径 導電率 分割数` | 直線ワイヤ (PEEC 抽出対象) |
| `port` | `port = n1 n2 Z0` | Zin 計算ポート |
| `frequency` | `frequency = f開始 f終了 分割数` | 周波数掃引 (分割数+1 点) |
| `nodetol` | `nodetol = 1e-8` | 座標マージ許容 [m] (省略時 1e-8) |
| `gmin` | `gmin = 0` | 全ノード対地コンダクタンス [S] (省略時 0) |

- ワイヤと回路素子の接続は `node` キーで行います。ワイヤ端点の座標が
  `nodetol` 以内で一致するとそのノード id に接続されます。
- ワイヤ内部の分割点は最大ノード id + 1 から記載順に自動採番されます
  (決定的)。
- 基準ノードはノード 0。ノード 0 が使われていない場合は port #1 の n2 を
  基準にします (ログに表示)。
- 未知のキーは無視されます (前方互換)。キー省略時の既定値は従来動作と
  一致します (後方互換)。

## Validation / 検証 (`data/sample/`)

| ケース | 解析解 | 許容 |
|---|---|---|
| `rlc_series.peec` 直列 RLC | Zin(f0) = 50 + j0、Xin(2f0) = 47.434 Ω | 0.5% |
| `wire_single.peec` 単線 1 m / a = 1 mm | L = 1.32038 µH、Rin = 5.4881 mΩ (ndiv = 1, 8) | 1% |
| `loop_square.peec` 正方ループ s = 0.2 m | L = (2μ0s/π)[ln(s/a)+a/s−0.774013] = 724.69 nH | 2% |

## License

MIT
