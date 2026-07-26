---
description: ビルドして解析解との検証 (peec_check.sh) を全ケース走らせる
allowed-tools: Bash(cmake:*), Bash(sh data/sample/peec_check.sh:*), Bash(ls:*)
---

ビルドして検証を走らせ、結果を報告してください。

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
sh data/sample/peec_check.sh "$PWD/bin/peec" /tmp/peec-check
```

報告は簡潔に:

- 全判定 OK なら「全 N 判定 OK」の 1 行 + 誤差が大きい上位 2〜3 件だけ
- NG があれば、その判定名・実測値・期待値・誤差と、
  `.claude/rules/physics-invariants.md` のどの不変条件の番人かを対応づける
- ビルドが失敗したら検証は走らせず、エラー箇所だけを示す

長い出力をそのまま貼らないこと。
