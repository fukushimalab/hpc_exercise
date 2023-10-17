# ネットワーク系演習II 高能率計算　1回目レポート
## 学籍番号：01234567

## 名前：名工 大

## 使用計算機環境情報
* CPU: [Intel Core i7 8650U 1.9GHz （ターボブースト4.2GHz） 4コア8スレッド](https://en.wikichip.org/wiki/intel/core_i7/i7-8650u)
* メモリ: 16GB
* OS: Windows 10 with Windows subsystem for Linux (WSL1)
* コンパイラ: g++

CPUの詳細情報は以下等で調べてください．
https://en.wikichip.org/wiki/WikiChip

# 課題１
ソースコードの計測部分前後にタイマー関数を追記して実行速度を計測した．

なんどか実行した結果の実行時間の例を以下に示す．

exercise 1: loop = 10, size = 30

|time    |ms  |
|--------|----|
|time   0|0.0084|
|time   1|0.0018|
|time   2|0.0018|
|time   3|0.0015|
|time   4|0.0016|
|time   5|0.0016|
|time   6|0.0018|
|time   7|0.0016|
|time   8|0.0016|
|time   9|0.0015|
|time avg|0.00164444|

考察：計算時間にはばらつきが観測できたが，おおむね一回目の計測時間が遅かった．

# 課題２
コンパイラオプションをO0~Ofastに変えてプログラムをコンパイルした．

結果は以下のようになった．

exercise 2: loop = 100, size = 1024

|option|time [ms]|
|------|---------|
|-O0   |2321.71  |
|-O1   |xxxxx    |
|・・・   |・・・      |
|-Ofast|92.4779 |

またコンパイル時間を下記に示す．
```shell-session 
fukushima@darkroom:~/hpc$ time make
g++ -std=c++0x -fopenmp -Wno-unused-result -O0 -march=native -mtune=native   -c -o utils/mat.o utils/mat.cpp
g++ -std=c++0x -fopenmp -Wno-unused-result -O0 -march=native -mtune=native   -c -o utils/mat_util.o utils/mat_util.cpp
g++ -std=c++0x -fopenmp -Wno-unused-result -O0 -march=native -mtune=native   -c -o main.o main.cpp
g++ -std=c++0x -fopenmp -Wno-unused-result -O0 -march=native -mtune=native -o hpc_exercise utils/mat.o utils/mat_util.o main.o

real    0m2.558s
user    0m1.484s
sys     0m1.016s

fukushima@darkroom:~/hpc$ time make
g++ -std=c++0x -fopenmp -Wno-unused-result -Ofast -march=native -mtune=native   -c -o utils/mat.o utils/mat.cpp
g++ -std=c++0x -fopenmp -Wno-unused-result -Ofast -march=native -mtune=native   -c -o utils/mat_util.o utils/mat_util.cpp
g++ -std=c++0x -fopenmp -Wno-unused-result -Ofast -march=native -mtune=native   -c -o main.o main.cpp
g++ -std=c++0x -fopenmp -Wno-unused-result -Ofast -march=native -mtune=native -o hpc_exercise utils/mat.o utils/mat_util.o main.o

real    0m3.573s
user    0m2.563s
sys     0m0.969s
```

考察：Ofastが最も高速に動作した．またコンパイル時間も最も長かった．O0とO1の実行速度の差は顕著だが，O1～3の差はごくわずかであり，優位な差を計測するにはより大きなサイズの行列で計測する必要があると考えられる．
そのほか，`-march=native -mtune=native`のオプションについても検証した結果，`-march=native`は指定しないとコンパイルエラーが出た．`-mtune=native`については，．．．

# 課題３
以下のようにコードを書き換えた．

書き換え前
```cpp
const float v = x.data[i];
ret.data[i] = 3.f * v * v * v * v * v * v
    + 3.f * v * v * v * v * v
    + 3.f * v * v * v * v
    + 3.f;

```
書き換え後
```cpp
const float v = x.data[i];
//ret.data[i] = XXXXXXX;
```

考察：ｘｘｘｘｘｘｘｘｘｘｘ



# 課題４

表１：表のサンプル
|項目|時間|
|------|-----|
|展開前|xx ms|
|展開後|xx ms|

．．．


# 課題２２

画像の張り付けサンプル

<img src="loofline.png" alt="ルーフライン" width="600px">
図ｘ：シングルスレッドとマルチスレッドのルーフライン．

参考までに，CSEは1.3TFLOPSくらい出ます（試すのは絶対に夜中で．）．

．．．

# 課題２９





# 画像処理課題１

# 画像処理課題２

1回目の課題はここまで．
画像処理の共通の課題である上記１，２を忘れずにやること．
これ以降の課題は2回目のレポートです．

