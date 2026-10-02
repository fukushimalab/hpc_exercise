# Makefileについて
配布されている`Makefile`をそのまま使用する．CSEサーバーの動作確認環境はg++ 13.3.0である．

このディレクトリで次のコマンドを実行する．

```shell
make
```

コンパイラは`CXX`，コンパイルオプションは`CXXFLAGS`に指定されている．最適化レベルを変更した場合は，`make -B`で全オブジェクトを再コンパイルすること．例えば，ファイルを編集せずに`-O3`を指定するには次のように実行する．

```shell
make -B CXXFLAGS='-std=c++0x -fopenmp -Wno-unused-result -march=native -O3'
```

`Makefile_for_CSE`には実行環境の確認表示が追加されている．通常のビルドにファイルのリネームは必要ない．

WindowsでLinux用のMakefileを使用する場合はWSLを利用する．AVX命令を扱う演習には，AVX対応のx86環境を使用すること．
