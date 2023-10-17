# 入門用 
二部実験など，短い時間で演習を行う時用．

1. プログラムをコンパイルする
テキストファイルに，hello worldを表示するプログラムを書き，コンパイルして実行する．
例えば，main.cというファイルを作成して下記のように書く．
```cpp
#include <stdio.h>
int main()
{
  printf("hello world");
  return 0;
}
```

そのあと，そのファイルが存在するディレクトリでターミナルを立ち上げ，下記コマンドを打つと，テキストファイルがコンパイルされ，C言語が機械語に変換される．
```console
gcc main.c
```
そうするとa.outというファイルができるので，
```console
./a.out
```
とするとターミナルでプログラムが実行できる．

2. 配列をmallocで確保し，配列の中身を0から順番に埋めよ．
malloc関数はstdlibヘッダがないと使えないためインクルードが必要である．

```cpp
  #include <stdio.h>
  #include <stdlib.h>
  int main()
  {
    const int size = 100000000;
    int* data = (int*)malloc(sizeof(int)*size);
    for(int i=0;i<size;i++)
    {
      data[i]=i;
    }
    return 0;
  }
```

3. この埋める作業の計算時間を計測せよ．
linuxには，プログラムの実行時間を計測するコマンドtimeがある．
下記で実行するとtimeに続くプログラムの実行時間を計測できる．
```
time ./a.out
```

この時，2番で用意した配列のサイズを変えて，大きいほど時間がかかることを確認せよ．

4. dataと同じサイズの配列data2を作り，そこに，逆順（i-size）でデータを用意せよ．

5. dataと同じサイズの配列data3を作り，dataとdata2を足した結果をdata3に格納し，計算時間を計測せよ．

6. data3は作らず，計算結果をdata2に上書きして入れて，計算時間を計測せよ．
これは，メモリの書き込みにキャッシュをうまく使うため，5よりも高速化するはずである．

7. 今使っているプログラムはマルチコアCPUである．複数のCPUを使って計算せよ．
C言語は，普通に書いただけでは，CPU一つしか使えない．
これを複数のCPUで計算することを簡単にできるようにするものにOpenMPがある．
例えば以下のようにするとforループをコア数の数だけに分割して，CPUの数だけ分業してくれる．
これを，足し算のプログラムに適用し，速度がどうなるか確認せよ．
ただし，コンパイルは，下記で行う．-fopenmpでOpenMPが使えるようになる．
```console
gcc main.c -fopenmp
```

```cpp
  int main()
  {
    const int size = 100000000;
    int* data = (int*)malloc(sizeof(int)*size);
  #pragma omp parallel for
     for(int i=0;i<size;i++)
     {
      data[i]=i;
    }
    return 0;
  }
```
８. 今使っているプログラムはSIMD演算が使用できる．SIMDとはの資料を読み，下記プログラムを実行し高速化したことを確認せよ．
なお，サイズは必ず８の倍数にすること．
これは，課題２のプログラムをSIMDによるベクトル演算化したものである．

```cpp
int main()
{
	const int size = 80000000;//must be 8 multple
	int* data = (int*)malloc(sizeof(int) * size);
	__m256i offset = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

	for (int i = 0; i < size; i += 8)
	{
		__m256i a = _mm256_add_epi32(_mm256_set1_epi32(i), offset);
		_mm256_storeu_si256((__m256i*)(data + i), a);
	}
	return 0;
}
```

コンパイルの方法は以下である．
```console
gcc main.c -mavx2
```

９．
下記コードは全ての値を総和するプログラムである．
このプログラムをどんな方法でもよいので高速化せよ．
ただし，数列を使って数式を展開する方法はやってもよいが，別途すべて計算する方法も作ること．
ヒントとして，初期化部分は高速化している．
また，コンパイルオプションを変えてみてもよい．
```cpp
#include <stdio.h>
#include <stdlib.h>
int main()
{
	const int size = 100000000;
	const int size32 = (size / 32) * 32;//必ず32の倍数にする
	int* data = (int*)malloc(sizeof(int) * size32);
	__m256i mstep = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);

#pragma omp parallel for schedule (dynamic)
	for (int i = 0; i < size32; i += 32)
	{
		_mm256_storeu_si256((__m256i*)(data + i + 0), _mm256_add_epi32(_mm256_set1_epi32(i + 0), mstep));
		_mm256_storeu_si256((__m256i*)(data + i + 8), _mm256_add_epi32(_mm256_set1_epi32(i + 8), mstep));
		_mm256_storeu_si256((__m256i*)(data + i + 16), _mm256_add_epi32(_mm256_set1_epi32(i + 16), mstep));
		_mm256_storeu_si256((__m256i*)(data + i + 24), _mm256_add_epi32(_mm256_set1_epi32(i + 24), mstep));
	}
	/*
	//上と同じコード
#pragma omp parallel for
	for (int i = 0; i < size32; i+=8)
	{
		_mm256_store_si256((__m256i*)(data + i), _mm256_add_epi32(_mm256_set1_epi32(i), mstep));
	}
	*/
	/*
	//上と同じコード
#pragma omp parallel for
	for (int i = 0; i < size32; i++)
	{
		data[i] = i;
	}
	*/

	int sum = 0;
	for (int i = 0; i < size32; i++)
	{
		sum += data[i];
	}
	printf("sum: %d", sum);// 887459712

	return 0;
}
```
