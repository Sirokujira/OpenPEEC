---
description: push 前の一括点検 (クリーンビルド・警告掃討・検証・スレッド一致・サンプル整合)
allowed-tools: Bash(cmake:*), Bash(gcc:*), Bash(sh data/sample/peec_check.sh:*), Bash(./bin/peec:*), Bash(grep:*), Bash(ls:*), Bash(git status:*), Bash(git diff:*)
---

push 前の点検を順に実行し、結果を表 1 つにまとめて報告してください。
途中で落ちたら、そこで止めて原因を示してください。

1. **クリーンビルド** — `rm -rf build && cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build`
   でエラー・警告が 0 件であること。
2. **警告掃討** — `for f in src/*.c; do gcc -std=c11 -Wall -Wextra -fopenmp -Iinclude -c "$f" -o /dev/null; done`
   が無出力であること。
3. **移植性チェック** — `sh .claude/hooks/check-portability.sh` が exit 0 であること。
4. **検証** — `sh data/sample/peec_check.sh "$PWD/bin/peec" /tmp/peec-preflight` が全判定 OK であること。
5. **スレッド非依存** — 代表ケースを `-n 1` と `-n 4` で実行し `zin.csv` が
   バイト一致すること (OpenMP を使う `plate_cap.peec` と `wire_single.peec` で確認)。
6. **サンプル整合** — `data/sample/*.peec` がすべて `peec_check.sh` から
   参照されていること (検証ケースの付け忘れ検出)。
7. **差分の確認** — `git status --short` と `git diff --stat` で、意図しない
   ファイル (build/, bin/, *.log, zin.csv) が混ざっていないこと。

報告フォーマット:

| 点検 | 結果 |
|---|---|
| クリーンビルド | OK / NG (理由) |
| ... | ... |

すべて OK なら最後に「push 可」とだけ書いてください。コミットや push は
明示的に頼まれるまで実行しないこと。
