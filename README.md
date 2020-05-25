# ネットワーク系演習II：ハイパフォーマンスコンピューティング
名古屋工業大学 情報工学科 ネットワーク系分野 3年後期 ネットワーク系演習II：ハイパフォーマンスコンピューティングの演習資料です．

[ドキュメントURL](https://fukushimalab.github.io/hpc_exercise/)

各プロジェクトは，Makefileでコンパイルできるようになっています．   

また，Visual Studio 2019でもコンパイルできるようにしています．
下記のファイルはVisual Studio 2019用のファイルです．必要ない人は無視してください．  
`*.sln` `*.vcxproj*`

ただしファイルの文字コードや改行コードには注意すること．
Linux用に改行コードはLFになっていますが，Visual Stdio2019上で実行するには，CRLFになっていないといけません．

# 動作確認
2020/5/26現在の動作確認

|OS等 |コンパイラ|備考|
|---|---------|---|
|CSE@384コア（名工大）|g++|〇：make|
|CSE@384コア（名工大）|icc|〇：make|
|20号館ローカル（名工大）|g++|〇：make|
|Linux|clang++|〇：make|
|Linux|g++|〇：make|
|Linux|clang++|〇：make|
|Windows|Visual Studio2019|〇：slnを開く|
|Windows|g++on WSL|〇：make|
|Windows|MinGW|△：※１|
|Windows|clang+VS2019|×：※２|
|Mac|clang++|△：※３|

* ※１インストーラでデフォルトではついてこないpthreadを必ずチェック．getclock_timeがないので`mat_util.h`の`＃USU_TIME_CHRONO`をコメントアウトを戻す．そのあとmake．
* ※２clangのOpenMPが有効化できずに動作していない．頑張ったら動くはず．
* ※３デフォルトのclangはOpenMPに対応していない可能性があるので，OpenMPに対応したg++に変更する．場合によっては，Makefileのg++のところはclang++に．



