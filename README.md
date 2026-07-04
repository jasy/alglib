![Action Status](https://github.com/jasy/alglib/workflows/C++%20CI/badge.svg)

alglib
======

algorithm libraries for programming contest

概要
----
アルゴリズムプログラミングコンテストでの覚書的なライブラリ。
主にTopCoderが対象。で他で使いにくいようなのも入るかもね。
Genericに書くことと速度はそれなりに気にすることにしておく。

環境
----
C++11を使っているので全機能対応の以下あたり推奨
- Clang 3.3
- GCC 4.8.1

ビルド / テスト
----
このリポジトリは Google Test 1.12.1 を submodule として使用しています。
Google Test のライセンスは `googletest/LICENSE` で BSD 3-clause ライセンスとして提示しています。

ビルド手順:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

カバレッジビルド:

```sh
cmake -S . -B build-coverage -DENABLE_COVERAGE=ON
cmake --build build-coverage
ctest --test-dir build-coverage
```

HTML レポートを生成するには `gcovr` が必要です。
```sh
gcovr -r . --exclude googletest --html --html-details -o coverage/index.html
```
