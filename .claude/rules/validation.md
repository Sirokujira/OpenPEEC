---
paths:
  - "data/sample/**"
  - "src/input_data.c"
  - ".github/workflows/*.yml"
---

# 検証ケースと入力キーの規則

## 入力キーを足すとき

- パースは `src/input_data.c` に追加する。
- **既定値は「キー省略時に従来動作と完全一致」**になるよう初期化する
  (後方互換)。`memset(p, 0, sizeof(peec_t))` の後に既定値を代入する。
- 未知のキーは無視する (前方互換)。既存の `else if` 連鎖の末尾に足すだけで
  この性質は保たれる。
- README の入力キー表に 1 行足す。

## 検証ケースを足すとき (新機能には必須)

1. `data/sample/<name>.peec` を作る。**先頭コメントに解析解の導出を書く** —
   使った公式、代入した数値、期待値、許容誤差の根拠。後から誰かが期待値を
   検算できることが目的。
2. `data/sample/peec_check.sh` に判定を足す。`chk <label> <actual> <expected> <tol>`
   を使う。値は `zin.csv` (機械読み取り用) か `peec.log` から awk で取る。
   ヘルパは `getL` (= Xin/2πf)、`getR` (= Rin)、`getC` (= 総容量) がある。
3. 期待値は**コードとは独立な出所**にする。教科書の公式、文献値、または
   別経路の実装 (面積分 vs 細線近似のような) との相互検証が望ましい。
   コード自身の出力を期待値にすると回帰テストにしかならない。
4. スクリプトは POSIX sh + awk/grep/sed のみ。CI は Windows でも同じ
   スクリプトを Git Bash (`shell: bash`) で走らせる。

## 許容誤差の付け方

| 種類 | 目安 |
|---|---|
| 解析式を厳密に再現するはずのもの | 0.1% |
| 近似式・級数展開との比較 | 0.5% |
| 漸近極限との比較 (補正項が残る) | 1% |
| 低次基底の離散化誤差を含むもの | 2〜10% |

許容を緩めるときは、なぜ緩いのかを .peec のコメントに書く。

## CI

`.github/workflows/ci.yml` は Linux / macOS / Windows (MSVC + Ninja) の 3 本 +
`sanitize` (Linux, ASan + UBSan + LeakSanitizer) の計 4 本。検証は全ジョブとも
同一の `peec_check.sh` を実行する。依存ライブラリが無いので vcpkg 等の手順は
不要。タグ `v*` push で Release にバイナリを添付する。
