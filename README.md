# ネットワーク系演習II：ハイパフォーマンスコンピューティング
名古屋工業大学 情報工学科 ネットワーク系分野 3年後期 ネットワーク系演習II：ハイパフォーマンスコンピューティングの演習資料です．詳しい説明は下記リンクから．

* [ドキュメントURL](https://fukushimalab.github.io/hpc_exercise/)
* [Intel Intrinsics](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html)
* [2部実験など短い演習用ページ] (simple.md)

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

また，Visual Studio 2019でもコンパイルできるようにしています．
下記のファイルはVisual Studio 2019用のファイルです．必要ない人は無視してください．  
`*.sln` `*.vcxproj*`

ただしファイルの文字コードや改行コードには注意すること．
Linux用に改行コードはLFになっていますが，Visual Stdio2019上で実行するには，CRLFになっていないといけません．

# 動作確認
2020/5/26現在の動作確認

|OS等 |コンパイラ|備考|
|---|---------|---|
|20号館ローカル（名工大）|g++|〇：make※１|
|CSE@384コア（名工大）|g++|〇：make※２|
|CSE@384コア（名工大）|icc|〇：make※２|
|Linux (Ubuntu)|g++|〇：make|
|Linux (Ubuntu)|clang++|〇：make|
|Windows|Visual Studio2019|〇：slnを開く|
|Windows|g++on WSL|〇：make|
|Windows|MinGW|△：※３|
|Windows|clang+VS2019|×：※４|
|Mac|clang++|△：※５|

* ※１：普段はこれがデフォルト
* ※２：普通の時間にこれで回すと，1～4年生の全ユーザのCPU資源が枯渇するので，他の授業に影響がでないように使うなら深夜．
* ※３：インストーラでデフォルトではついてこないpthreadを必ずチェック．getclock_timeがないので`mat_util.h`の`＃USU_TIME_CHRONO`をコメントアウトを戻す．そのあとmake．この場合，1ms以下の精度がないので，タイマーの測り方を変更すること．
* ※４：Visual StudioのMSBuild用のclangのOpenMPが有効化できずに動作していない．頑張ったら動くはず．普通のLLVM+clangなら動くはず（検証していない）
* ※５：デフォルトのclangはOpenMPに対応していない可能性が高いので，OpenMPに対応したg++に変更する．場合によっては，Makefileのg++のところはclang++に．[インストール用の参考ページ](https://mem-archive.com/2019/08/17/post-2038/)．

AVX命令が前提なので，IntelかAMDのCPUが必要です．ARMのCPUでは動きません．
Apple M1を使っている人は，どうやっても対応できないのでCSEの計算機を使ってください．

CSEの情報（2023年～）：
* 各自端末 [Intel Core-i5 11500T](https://www.intel.co.jp/content/www/jp/ja/products/sku/212272/intel-core-i511500t-processor-12m-cache-up-to-3-90-ghz/specifications.html)（1.5GHz，6コア，12スレッド）2号館
* 各自端末 [Intel Core-i5 11500](https://www.intel.co.jp/content/www/jp/ja/products/sku/212277/intel-core-i511500-processor-12m-cache-up-to-4-60-ghz/specifications.html)（2.7GHz，6コア，12スレッド）20号館
* サーバー [AMD EPYC 7763](https://www.amd.com/ja/products/cpu/amd-epyc-7763) 64コア128スレッドx2
