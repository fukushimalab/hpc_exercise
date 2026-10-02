# ネットワーク系演習II：ハイパフォーマンスコンピューティング
名古屋工業大学 情報工学科 ネットワーク系分野 3年後期 ネットワーク系演習II：ハイパフォーマンスコンピューティングの演習資料です．詳しい説明は下記リンクから．

* [ドキュメントURL](https://fukushimalab.github.io/hpc_exercise/)
* [Intel Intrinsics](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html)
* [2部実験など短い演習用ページ](./simple.md)

# Todo
## 2023
* [ ] 画像表示関数の作成（tmp.pngを書きだして，popenあたりで，開くだけ．）

# ディレクトリ構成
```
/                       ドキュメントファイル
/docfig                 ドキュメントファイルに必要な画像ファイル
/src                    演習用のファイル一覧
/src/hpc-exercise       最初の演習用のファイル
/src/image-processing   画像処理の演習用のファイル
/src/boxfilter          画像処理の選択課題用のファイル
/report_sample_md       レポート提出に使うサンプル用のマークダウンファイル（要pdf変換）．もちろんtexでもwordでもpdfに変換して提出してもらえれば，書くツールはなんでもよい．
```

src内の各プロジェクトは，Makefileでコンパイルできるようになっています．   

Windows用にはVisual Studioのプロジェクトファイル（`*.sln`，`*.vcxproj*`）も提供しています．本演習の動作確認環境はCSEのg++です．ファイルはLF改行で管理しています．

# 動作確認環境
2026年10月2日にCSEサーバーで確認した環境です．

|項目|環境|
|---|---|
|OS|Linux（x86_64）|
|コンパイラ|g++ 13.3.0（Ubuntuパッケージ：13.3.0-6ubuntu2~24.04.1）|
|CPU|AMD EPYC 9745 × 2|
|コア数|256物理コア・512論理CPU|
|ビルド|各プロジェクトのMakefileを使用|

プロジェクトのルートから，次のようにコンパイルします．

```shell
make -C src/hpc-exercise
```

最適化オプションを変更した場合は，`make -B`で全オブジェクトを再コンパイルしてください．CSEのg++ 13.3では，`-O2`でも自動ベクトル化が有効です．オプションの説明とアセンブリの観察方法は[演習資料](index.md)に記載しています．

WindowsでLinux環境を利用する場合はWSLを使用できます．AVX命令を扱う演習は，AVXに対応するx86 CPUが必要です．Apple SiliconなどARM環境を使用している場合は，CSEサーバーで実行してください．
