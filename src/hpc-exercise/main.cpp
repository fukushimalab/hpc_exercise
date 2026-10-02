#if 1
#define EX1
#define EX2
#define EX3
#define EX4
#define EX5
#define EX6
#define EX65
#define EX7
#define EX8
#define EX9
#define EX10
#define EX11
#define EX12
#define EX13
#define EX14
#define EX15
#define EX16
#define EX17
#define EX18
#define EX19
#define EX20
#define EX21
#define EX22
#define EX23
#define EX24
#define EX25
#define EX26
#define EX27
#define EX28
#define EX29
#else
//一番先頭の#if 1を#if 0に変えると，ここで定義した課題だけがコンパイルされるようになる．
//対象外の課題をコンパイルから除外し，アセンブリ出力を絞るために使用する．

#define EX9

#endif

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include "utils/simd_util.h"
#include "utils/mat_util.h"
#include <iostream>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <omp.h>
#ifdef EX22
#include "utils/loofline_performace.h"
#endif
#ifdef EX27
void inline _mm256_transpose_8x8_ps(__m256* dst, const __m256* src);
#endif
#ifdef EX12
//課題12用：回転を繰り返す段数．変更したら再ビルドすること．
#define ROT_STAGES 8
//関数のinline展開を禁止する指定．GCC/clangとVisual Studio（MSVC）で書き方が異なる．
#ifdef _MSC_VER
#define NOINLINE __declspec(noinline)
#else
#define NOINLINE __attribute__((noinline))
#endif
NOINLINE void rot_noinline(double a, double b, double& x, double& y, double c, double s);
void rot_function(double a, double b, double& x, double& y, double c, double s);
inline void rot_inline(double a, double b, double& x, double& y, double c, double s);
#endif
#ifdef EX2
void GEMM(Mat_32F& a, Mat_32F& b, Mat_32F& c);
#endif


int main(const int argc, const char** argv)
{
	//ヘルプの表示
	if (argc == 1)
	{
		std::cout << "Usage: hpc_exercise exercise_nums loop_nums mat_size" << std::endl;
		std::cout << "ex: run exercise 1 with default parameters" << std::endl;
		std::cout << "./hpc_exercise 1" << std::endl;
		std::cout << "ex: run exercise 2 with 10 iterations and default matrix size" << std::endl;
		std::cout << "./hpc_exercise 2 10" << std::endl;
		std::cout << "ex: run exercise 3 with 10 iterations and matrix size is 32 x 32" << std::endl;
		std::cout << "./hpc_exercise 3 10 32" << std::endl;
		std::cout << std::endl;
	}
	//ヘルプの表示ここまで．この間のcoutは削除可能．

	const int exercise = (argc < 2) ? 0 : atoi(argv[1]);
	const int arg_loop = (argc < 3) ? 0 : atoi(argv[2]);
	const int arg_size = (argc < 4) ? 0 : atoi(argv[3]);

#ifdef EX1
	//課題1
	//(1) 行列の"要素"の積和演算AX + Bを計算するプログラムにおいて，行列積と和それぞれの実行時間をタイマーを挟むことで測定せよ．
	//なお，デフォルトのサンプルコードには行列積の後，和をとるプログラムがサンプルとして書かれており，タイマーもすでに記述済みである．
	//(2) また，１度目の計算時間は，２回目以降よりも統計的に遅いことを確認せよ．
	if (exercise == 1)
	{
		//慣れるためのの問題．timerをコメントアウトしてみたり行列の表示をしたりすることを習得する．
		const int default_loop = 10;
		const int default_size = 3;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int row = (arg_size == 0) ? default_size : arg_size;
		const int col = row;
		std::cout << "exercise 1: loop = " << loop << ", size = " << row << std::endl;

		Mat_32F a(row, col);
		mat_rand(a, 0.f, 100.f);
		Mat_32F x(row, col);
		mat_rand(x, 0.f, 100.f);
		Mat_32F b(row, col);
		mat_rand(b, 0.f, 100.f);
		Mat_32F ret(row, col);
		mat_zero(ret);

		CalcTime t;
		std::cout << "|time    |ms  |" << std::endl;
		std::cout << "|--------|----|" << std::endl;
		for (int i = 0; i < loop; i++)
		{
			t.start();	//時間計測開始
			ret = mat_add(mat_mul(a, x), b);
			t.end();	//時間計測終了
			std::cout << "|time   " << i << "|" << t.getLastTime() << "|" << std::endl;
		}
		std::cout << "|time avg|" << t.getAvgTime() << "|" << std::endl << std::endl;

		//行列の中身を表示．邪魔になったらコメントアウトすること．
		a.show();
		x.show();
		b.show();
		ret.show();

		return 0;
	}

#endif
#ifdef EX2
	//課題2
	//行列積を計算するプログラムにおいて，コンパイラオプションを変えて計算速度の計測し，その違いを観測せよ．
	//なお，行列の要素積ではなく行列積であり，呼ぶ関数はGEMMである．
	if (exercise == 2)
	{
		//コンパイルオプションを変えることに慣れる問題．空欄埋めはない．
		const int default_loop = 100;
		const int default_size = 512;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 2: loop = " << loop << ", size = " << size << std::endl << std::endl;

		Mat_32F a(size, size);
		mat_rand(a, 0.f, 100.f);
		Mat_32F b(size, size);
		mat_rand(b, 0.f, 100.f);
		Mat_32F c(size, size);
		mat_zero(c);

		CalcTime t;

		for (int i = 0; i < loop; i++)
		{
			mat_zero(c);
			t.start();
			GEMM(a, b, c);
			t.end();
			//std::cout << "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|option|time [ms]|" << std::endl;
		std::cout << "|------|---------|" << std::endl;
		std::cout << "|-Oxxxx|" << t.getAvgTime() << "|" << std::endl;

		return 0;
	}

#endif
#ifdef EX3
	//課題3
	//小さな行列に対して，各要素を
	//3x^6+3x^5+3*x^4+3*x^3+3
	//する計算するプログラムを作成し，乗算回数を削減する前と後で計算速度を比較せよ．
	if (exercise == 3)
	{
		//空欄埋め問題
		const int default_loop = 100000;
		const int default_size = 64;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 3: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;
		Mat_32F x(row, col);
		mat_rand(x, 0.f, 1.f);
		Mat_32F ans(row, col);
		Mat_32F y(row, col);
		mat_zero(ans);
		mat_zero(y);

		//before
		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			//そのまま
			t.start();
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算 ansに書き込み
				const float v = x.data[i];
				ans.data[i] = 3.f * v * v * v * v * v * v
					+ 3.f * v * v * v * v * v
					+ 3.f * v * v * v * v
					+ 3.f * v * v * v
					+ 3.f;
			}
			t.end();
			//std::cout << "before: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|method|time [ms]|" << std::endl;
		std::cout << "|------|---------|" << std::endl;
		std::cout << "|before|" << t.getAvgTime() << "|" << std::endl;

		//after
		for (int k = 0; k < loop; k++)
		{
			//ホーナー法つかった場合
			//3 * ((x * x * x * (x * (x * (x + 1) + 1) + 1)) +1);
			t.start();
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算 yに書き込み
				const float v = x.data[i];
				//y.data[i] = XXXXXXXX;
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|after |" << t.getAvgTime() << "|" << std::endl;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;

		std::cout << "diff from ans: " << mat_diff(ans, y) / double(ans.cols * ans.rows) << std::endl << std::endl;

		return 0;
	}

#endif
#ifdef EX4
	//課題4
	//小さな行列に対して，各要素を下記の定数倍するプログラムを作成し，数式の展開前後で計算速度を比較せよ．
	//(2π+sqrt(5)+0.5^2)x
	if (exercise == 4)
	{
		//空欄埋め問題
		const int default_loop = 10000;
		const int default_size = 64;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 4: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;
		Mat_32F x(row, col);
		mat_rand(x, 0.f, 100.f);
		Mat_32F ans(row, col);
		Mat_32F y(row, col);
		mat_zero(ans);
		mat_zero(y);

		CalcTime t;

		//before
		for (int k = 0; k < loop; k++)
		{
			//毎回計算する場合
			t.start();
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算 ansに書き込み
				//ans.data[i] = XXXXXXXX
			}

			t.end();
			//std::cout << "before: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|method |time [ms]|" << std::endl;
		std::cout << "|-------|---------|" << std::endl;
		std::cout << "|inline |" << t.getAvgTime() << "|" << std::endl;

		//after
		for (int k = 0; k < loop; k++)
		{
			//先に計算する場合
			t.start();
			const int s = x.rows * x.cols;
			//XXXXXXXX //定数値を計算
			for (int i = 0; i < s; i++)
			{
				//計算 yに書き込み
				//y.data[i] = XXXXXXXX;
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|precomp|" << t.getAvgTime() << "|" << std::endl;


		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;

		std::cout << "diff from ans: " << mat_diff(ans, y) / double(ans.cols * ans.rows) << std::endl << std::endl;

		return 0;
	}

#endif
#ifdef EX5
	//課題5
	//小さな行列に対して，各要素を3.141592で除算する計算するプログラムを作成し，除算を削減する前と後で計算速度を比較せよ．
	// 本課題は行列を2重ループで処理している．これをループをつぶした場合に，効果がどれくらいはっきりするか比較せよ． なお，このループ構造については，演習用コードと行列演算ライブラリの説明におけるチュートリアル課題3で示してある通り．
	//大きな行列で行うと，効果が少ない可能性があるため注意すること．
	//-O2と-O3の両方で計測し，課題4.5の方法でアセンブリを観察して考察せよ．
	if (exercise == 5)
	{
		//空欄埋め問題
		const int default_loop = 10000;
		const int default_size = 128;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 5: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;
		Mat_32F x(row, col);
		mat_rand(x, 1.f, 100.f);
		Mat_32F ans(row, col);
		Mat_32F y(row, col);
		mat_zero(ans);
		mat_zero(y);

		CalcTime t;
		//before
		for (int k = 0; k < loop; k++)
		{
			//除算の場合
			t.start();
			for (int j = 0; j < x.rows; j++)
			{
				for (int i = 0; i < x.cols; i++)
				{
					//計算 ansに書き込み
					//ans.data[XXXXXXXX]=xxxx
				}
			}
			t.end();
			//std::cout << "div: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|method |time [ms]|" << std::endl;
		std::cout << "|-------|---------|" << std::endl;
		std::cout << "|div 2lp|" << t.getAvgTime() << "|" << std::endl;

		//after
		for (int k = 0; k < loop; k++)
		{
			//乗算にした場合
			t.start();
			for (int j = 0; j < x.rows; j++)
			{
				for (int i = 0; i < x.cols; i++)
				{
					//計算
					//y.data[XXXXXXXX]=xxxxx
				}
			}
			t.end();
			//std::cout << "mul: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|mul 2lp|" << t.getAvgTime() << "|" << std::endl;

		//1重ループの場合
		for (int k = 0; k < loop; k++)
		{
			//1重ループの除算の場合
			t.start();
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算
				//ans.data[XXXXXXXX]=xxxx
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|div 1lp|" << t.getAvgTime() << "|" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//1重ループの乗算にした場合
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				//y.data[XXXXXXXX]=xxxxx
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|mul 1lp|" << t.getAvgTime() << "|" << std::endl;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;

		std::cout << "diff from ans: " << mat_diff(ans, y) / double(ans.cols * ans.rows) << std::endl << std::endl;

		return 0;
	}

#endif
#ifdef EX6
	//課題6
	//小さな4つのfloat型の行列A, B, C, D（要素は1〜100）について，同じ位置の要素から下記の式を計算せよ．
	//a / b / c / d / b / c / d / b / c / d / b / c / d
	//除算は左から順に行う．aをb，c，dで順に割る処理を4回重ねる．
	//式を変形して除算回数を削減し，変形前後の計算速度を-O0と-O3で比較せよ．
	//課題4.5の方法で生成されたアセンブリを観察して結果を考察せよ．
	if (exercise == 6)
	{
		//空欄埋め問題
		const int default_loop = 10000;
		const int default_size = 64;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 6: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;

		Mat_32F a(row, col);
		mat_rand(a, 1.f, 100.f);
		Mat_32F b(row, col);
		mat_rand(b, 1.f, 100.f);
		Mat_32F c(row, col);
		mat_rand(c, 1.f, 100.f);
		Mat_32F d(row, col);
		mat_rand(d, 1.f, 100.f);
		Mat_32F ans(row, col);
		Mat_32F ret(row, col);
		mat_zero(ans);
		mat_zero(ret);

		CalcTime t;
		//before
		for (int k = 0; k < loop; k++)
		{
			//普通に計算
			t.start();
			for (int j = 0; j < ans.rows; j++)
			{
				for (int i = 0; i < ans.cols; i++)
				{
					//計算 ansに書き込み
					//ans.data[XXXXXXXX]=xxxx
				}
			}
			t.end();
			//std::cout << "before: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|method|time [ms]|" << std::endl;
		std::cout << "|------|---------|" << std::endl;
		std::cout << "|before|" << t.getAvgTime() << "|" << std::endl;

		//after
		for (int k = 0; k < loop; k++)
		{
			//除算を削減した場合
			t.start();
			for (int j = 0; j < ret.rows; j++)
			{
				for (int i = 0; i < ret.cols; i++)
				{
					//計算 retに書き込み
					//ret.data[XXXXXXXX]=xxxx
				}
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|after |" << t.getAvgTime() << "|" << std::endl;


		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		std::cout << "diff from ans: " << mat_diff(ans, ret) / double(ans.cols * ans.rows) << std::endl << std::endl;
		return 0;
	}

#endif
#ifdef EX65
	//課題6.5：データIOと演算コスト
	//定数のコピー，入力のコピー，単純な加算，入力値の繰り返し加算を比較せよ．
	//-O3で行列サイズと加算回数を変え，生成されたアセンブリと実行時間を考察せよ．
	//繰り返し加算の回数Nは下記のADDITION_COUNTで指定する（4，8，16）．変更したら再ビルドすること．
	//コンパイル時にNが決まっていると，コンパイラが加算のループを展開し，要素ループのSIMD化を保てる．
#define ADDITION_COUNT 4
	if (exercise == 65)
	{
		const int loop = (arg_loop == 0) ? 1001 : arg_loop;
		const int size = (arg_size == 0) ? 64 : arg_size;
		const int addition_count = ADDITION_COUNT;
		if (loop < 2 || size < 1 || addition_count < 1)
		{
			std::cout << "exercise 65: loop >= 2, size >= 1, additions >= 1" << std::endl;
			return 1;
		}
		std::cout << "exercise 65: loop = " << loop << ", size = " << size
			<< ", additions = " << addition_count << std::endl << std::endl;
		Mat_32F x(size, size);
		Mat_32F y(size, size);
		mat_rand(x, 1.f, 100.f);
		mat_zero(y);
		const int elements = x.rows * x.cols;
		CalcTime t;
		std::cout << "|method|time [ms]|" << std::endl;
		std::cout << "|------|---------|" << std::endl;

		//定数のコピー
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int i = 0; i < elements; i++)
			{
				//XXXXXXXX y.data[i]に1.fを書き込む．
			}
			t.end();
		}
		std::cout << "|constant copy|" << t.getAvgTime() << "|" << std::endl;
		std::cout << "result sample: " << y.data[0] << std::endl;

		//入力のコピー
		mat_zero(y);//前の方式の結果を消す（計測外）
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int i = 0; i < elements; i++)
			{
				//XXXXXXXX x.data[i]をy.data[i]にコピーする．
			}
			t.end();
		}
		std::cout << "|input copy|" << t.getAvgTime() << "|" << std::endl;
		std::cout << "result sample: " << y.data[0] << std::endl;

		//単純な加算
		mat_zero(y);//前の方式の結果を消す（計測外）
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int i = 0; i < elements; i++)
			{
				//XXXXXXXX x.data[i]に1.fを加え，y.data[i]に書き込む．
			}
			t.end();
		}
		std::cout << "|add|" << t.getAvgTime() << "|" << std::endl;
		std::cout << "result sample: " << y.data[0] << std::endl;

		//非定数の繰り返し加算
		mat_zero(y);//前の方式の結果を消す（計測外）
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int i = 0; i < elements; i++)
			{
				//XXXXXXXX x.data[i]を一度読み出し，その値を0.fにaddition_count回加える．
				//最終結果だけをy.data[i]に書き込む．
			}
			t.end();
		}
		std::cout << "|repeated add|" << t.getAvgTime() << "|" << std::endl;
		std::cout << "result sample: " << y.data[0] << std::endl;
		return 0;
	}
#endif
#ifdef EX7
	//課題7
	//行列の各要素を2乗，3乗，4乗．．．n乗としたときに，mulで作ったものとpowで作ったものの速度を比較せよ．
	//また，nがいくつの時にmulのほうが速くなるのか（それとも常時powのほうが遅い・速いのか）比較せよ．
	//ただし，累乗をすると値が大きくなるため，浮動小数点の最大値を超える可能性がある．その時はrandの初期化の値を小さくせよ．
	//なお，この課題はコンパイルオプションによって結果が大きく変わる．
	//ansにpowの結果を，retにmulの結果を入れること．
	if (exercise == 7)
	{
		//変数を書き換える問題．空欄埋めはない．
		const int default_loop = 10000;
		const int default_size = 64;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 7: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;

		Mat_64F x(row, col);
		mat_rand(x, 1.0, 1.01);
		Mat_64F ans(row, col);
		mat_zero(ans);
		Mat_64F ret(row, col);
		mat_zero(ret);

		CalcTime t;

		//2乗をpow計算
		int pow_n = 2;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			//powで計算
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算（ansに入れる）
				ans.data[i] = pow(x.data[i], pow_n);
			}
			t.end();
		}
		std::cout << "|method|time [ms]|" << std::endl;
		std::cout << "|------|---------|" << std::endl;
		std::cout << "|pow  2|" << t.getAvgTime() << "|" << std::endl;

		//2乗をmulで計算
		for (int j = 0; j < loop; j++)
		{
			t.start();
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算（retに入れる）
				//ret.data[i] = xxxxx;
			}
			t.end();
		}
		std::cout << "|mul  2|" << t.getAvgTime() << "|" << std::endl;


		//3乗をpow計算
		pow_n = 3;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			//powで計算
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算
				ans.data[i] = pow(x.data[i], pow_n);
			}
			t.end();
		}
		std::cout << "|pow  3|" << t.getAvgTime() << "|" << std::endl;

		//3乗をmulで計算
		for (int j = 0; j < loop; j++)
		{
			t.start();
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算
				//ret.data[i] = xxxx;
			}
			t.end();
		}
		std::cout << "|mul  3|" << t.getAvgTime() << "|" << std::endl;


		//4乗をpow計算
		pow_n = 4;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			//powで計算
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算
				ans.data[i] = pow(x.data[i], pow_n);
			}
			t.end();
		}
		std::cout << "|pow  4|" << t.getAvgTime() << "|" << std::endl;

		//4乗をmulで計算
		for (int j = 0; j < loop; j++)
		{
			t.start();
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算
				//ret.data[i] = xxxx;
			}
			t.end();
		}
		std::cout << "|mul  4|" << t.getAvgTime() << "|" << std::endl;
		//std::cout << "diff from ans: " << mat_diff(ans, ret) / double(ans.cols * ans.rows) << std::endl << std::endl;


		//n乗をpow計算(nを変える問題．powは32をサンプルとして入力している．)
		pow_n = 32;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			//powで計算
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//計算
				ans.data[i] = pow(x.data[i], pow_n);
			}
			t.end();
		}
		std::cout << "|pow " << pow_n << "|" << t.getAvgTime() << "|" << std::endl;

		//n乗をmulで計算（nは任意）
		for (int j = 0; j < loop; j++)
		{
			t.start();
			const int s = x.rows * x.cols;
			for (int i = 0; i < s; i++)
			{
				//ret.data[i] = xxxx
			}
			t.end();
		}
		std::cout << "|mul " << pow_n << "|" << t.getAvgTime() << "|" << std::endl;


		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		std::cout << "diff from ans: " << mat_diff(ans, ret) / double(ans.cols * ans.rows) << std::endl << std::endl;
		return 0;
	}

#endif
#ifdef EX8
	//課題8
	//2つの行列の和を`char, unsigned char, short, int, float, double`で計算しそれぞれ比較せよ．
	//なお，大きい行列サイズでないと，効果がでない場合がある．
	if (exercise == 8)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 1024;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 8: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;

		CalcTime t;

		Mat_8U a_8u(row, col);
		Mat_8U b_8u(row, col);
		Mat_8U ret_8u(row, col);
		mat_rand(a_8u, 0, 63);//64は64+64=128でcharがオーバーフロー
		mat_rand(b_8u, 0, 63);//64は64+64=128でcharがオーバーフロー

		Mat_8S a_8s(a_8u);
		Mat_8S b_8s(b_8u);
		Mat_8S ret_8s(row, col);

		Mat_16S a_16s(a_8u);
		Mat_16S b_16s(b_8u);
		Mat_16S ret_16s(row, col);

		Mat_32S a_32s(a_8u);
		Mat_32S b_32s(b_8u);
		Mat_32S ret_32s(row, col);

		Mat_32F a_32f(a_8u);
		Mat_32F b_32f(b_8u);
		Mat_32F ret_32f(row, col);

		Mat_64F a_64f(a_8u);
		Mat_64F b_64f(b_8u);
		Mat_64F ret_64f(row, col);

		//unsigend char
		for (int i = 0; i < loop; i++)
		{
			t.start();
			//unsigned char
			//参考に埋めた．ret_8u = mat_add(a_8u, b_8u)のスタイルだと，内部でメモリを確保するため挙動が異なる．
			mat_add(a_8u, b_8u, ret_8u);
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|method|time [ms]|" << std::endl;
		std::cout << "|------|---------|" << std::endl;
		std::cout << "|8U    |" << t.getAvgTime() << "|" << std::endl;

		//char
		for (int i = 0; i < loop; i++)
		{
			t.start();
			// char
			//XXXXXXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|8S    |" << t.getAvgTime() << "|" << std::endl;

		//short
		for (int i = 0; i < loop; i++)
		{
			t.start();
			//short
			//XXXXXXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|16S   |" << t.getAvgTime() << "|" << std::endl;


		//int
		for (int i = 0; i < loop; i++)
		{
			t.start();
			//int
			//XXXXXXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|32S   |" << t.getAvgTime() << "|" << std::endl;


		//float
		for (int i = 0; i < loop; i++)
		{
			t.start();
			//float
			//XXXXXXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|32F   |" << t.getAvgTime() << "|" << std::endl;


		//double
		for (int i = 0; i < loop; i++)
		{
			t.start();
			//double
			//XXXXXXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "|64F   |" << t.getAvgTime() << "|" << std::endl;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;

		Mat_64F ans(ret_8u);
		Mat_64F temp = Mat_64F(ret_8s);
		std::cout << "diff  8S from 8U: " << mat_diff(temp, ans) / double(ans.cols * ans.rows) << std::endl;
		temp = Mat_64F(ret_16s);
		std::cout << "diff 16S from 8U: " << mat_diff(temp, ans) / double(ans.cols * ans.rows) << std::endl;
		Mat_64F temp2 = Mat_64F(ret_32s);
		std::cout << "diff 32S from 8U: " << mat_diff(temp2, ans) / double(ans.cols * ans.rows) << std::endl;
		temp = Mat_64F(ret_32f);
		std::cout << "diff 32F from 8U: " << mat_diff(temp, ans) / double(ans.cols * ans.rows) << std::endl;
		temp = Mat_64F(ret_64f);
		std::cout << "diff 64F from 8U: " << mat_diff(temp, ans) / double(ans.cols * ans.rows) << std::endl << std::endl;
		return 0;
	}

#endif
#ifdef EX9
	//課題9
	//intの行列を整数で2倍，浮動小数点で2.f倍,整数を１ビットだけビットシフトすることで2倍する場合の計算速度を比較せよ．
	//また，intの行列を整数で2で除算する場合，浮動小数点で2で除算する場合，浮動小数点の0.5で乗算する場合，１ビットだけビットシフトすることで1/2倍する場合の速度を比較せよ．
	//加えて，floatの行列で，2.0で除算する場合と0.5で乗算する場合を比較せよ．
	//なお，浮動小数点で乗算する場合は整数の場合よりも遅いことが多いが，新しい計算機ではそれほど違わないことも多い．
	//また，小さい行列サイズでないと，効果がでない場合がある．
	//差が出なかった場合は，課題4.5の方法でアセンブリを観察し，各方式が同じ命令になっているかを確認せよ．
	if (exercise == 9)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 1024;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 9: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;

		//int
		Mat_32S x_32s(row, col);
		Mat_32S ret_32s(row, col);


		CalcTime t;
		//2x mul
		for (int k = 0; k < loop; k++)
		{
			//2倍 乗算
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32s.data[i]=xxxx
			}
			t.end();
		}
		std::cout << "|method      |time [ms]|" << std::endl;
		std::cout << "|------------|---------|" << std::endl;
		std::cout << "|2 mul  (32S)|" << t.getAvgTime() << "|" << std::endl;

		//2.0x mul
		for (int k = 0; k < loop; k++)
		{
			//2.f倍 乗算(float)
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32s.data[i] =xxx
			}
			t.end();
		}
		std::cout << "|2 mul  (32F)|" << t.getAvgTime() << "|" << std::endl;

		//2x bit shift
		for (int k = 0; k < loop; k++)
		{
			//2倍 ビットシフト
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32s.data[i] =xxx
			}
			t.end();
		}
		std::cout << "|1 shift(32S)|" << t.getAvgTime() << "|" << std::endl;
		std::cout << "|------------|---------|" << std::endl;


		//1/2 div int
		for (int k = 0; k < loop; k++)
		{
			//1/2 除算
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32s.data[i] =xxx
			}
			t.end();
		}
		std::cout << "|2 div  (32S)|" << t.getAvgTime() << "|" << std::endl;

		//1/2 div double
		for (int k = 0; k < loop; k++)
		{
			//1/2.0 除算
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32s.data[i] =xxx
			}
			t.end();
		}
		std::cout << "|2 div  (64F)|" << t.getAvgTime() << "|" << std::endl;

		//1/2 -> mul 0.5
		for (int k = 0; k < loop; k++)
		{
			//1/2 0.5乗算で実現
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32s.data[i] =xxx
			}
			t.end();
		}
		std::cout << "|0.5 mul(64F)|" << t.getAvgTime() << "|" << std::endl;

		//1/2->bit shift
		for (int k = 0; k < loop; k++)
		{
			//1/2 ビットシフト
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32s.data[i] =xxx
			}
			t.end();
		}
		std::cout << "|1 shift(32S)|" << t.getAvgTime() << "|" << std::endl;
		std::cout << "|------------|---------|" << std::endl;


		//float
		Mat_32F x_32f(row, col);
		mat_rand(x_32f, 0, 100);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);

		//1/2 div
		for (int k = 0; k < loop; k++)
		{
			//1/2.f 除算
			t.start();
			const int size = ret_32f.cols * ret_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|2 div  (32F)|" << t.getAvgTime() << "|" << std::endl;

		//1/2 -> mul 0.5
		for (int k = 0; k < loop; k++)
		{
			//1/2 0.5f乗算
			t.start();
			const int size = ret_32f.cols * ret_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|0.5 mul(32F)|" << t.getAvgTime() << "|" << std::endl;


		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif 
#ifdef EX10
	//課題10
	//floatの行列を3.141倍する場合と，intの行列を3.141倍する場合と，intの行列を3.141倍を固定小数点で行う場合とshortの行列を3.141倍を固定小数点で行う場合で計算し比較せよ．
	//ただし，shortの配列ではオーバーフローに注意せよ．
	//例えば，3.141を10ビットシフトで表現する場合，3.141 * 1024 = 3216であり，short maxは32768であるため，入力の値は最大10までしかとることができない．課題のコードでは0～100の乱数であるため，適宜シフトの量を工夫せよ．
	if (exercise == 10)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 1024;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 10: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;

		CalcTime t;

		Mat_32F x_32f(row, col);
		mat_rand(x_32f, 0, 100);
		Mat_32F ans_32f(row, col);
		Mat_32F ret_32f(row, col);
		//float
		for (int k = 0; k < loop; k++)
		{
			//浮動小数点
			t.start();
			const int size = ret_32f.cols * ret_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ans_32f.data[i] =
			}
			t.end();
		}
		std::cout << "|method     |time [ms]|" << std::endl;
		std::cout << "|-----------|---------|" << std::endl;
		std::cout << "|32F mul 32F|" << t.getAvgTime() << "|" << std::endl;


		Mat_32S x_32s = Mat_32S(x_32f);
		Mat_32S ret_32s(row, col);

		//int floating point mul
		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32s.data[i] =
			}
			t.end();
		}
		std::cout << "|32S mul 32F|" << t.getAvgTime() << "|" << std::endl;

		//int fix
		for (int k = 0; k < loop; k++)
		{
			//固定小数点
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32s.data[i]
			}
			t.end();
		}
		std::cout << "|32S mul fix|" << t.getAvgTime() << "|" << std::endl;


		Mat_16S x_16s = Mat_16S(x_32f);
		Mat_16S ret_16s(row, col);

		//short fix
		for (int k = 0; k < loop; k++)
		{
			//固定小数点
			t.start();
			const int size = ret_16s.cols * ret_16s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_16s.data[i] =
			}
			t.end();
		}
		std::cout << "|16S mul fix|" << t.getAvgTime() << "|" << std::endl;


		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		Mat_32F temp = Mat_32F(ret_32s);
		std::cout << "32S diff from ans: " << mat_diff(ans_32f, temp) / double(ans_32f.cols * ans_32f.rows) << std::endl;
		temp = Mat_32F(ret_16s);
		std::cout << "16S diff from ans: " << mat_diff(ans_32f, temp) / double(ans_32f.cols * ans_32f.rows) << std::endl << std::endl;
		return 0;
	}
#endif
#ifdef EX11
	//課題11
	//floatの行列への定数値(3.141592f)の四則演算(加算，乗算，除算)と，floatの行列の値自体に`sqrt, sin, cos, exp, log`関数の適用した場合の計算時間を比較せよ．
	//また，`sin, cos, exp, log, sqrt`計算はテーブル参照も作成した場合についても比較せよ．
	//なお，環境によっては，演算したほうが速い演算もある可能性がある．
	if (exercise == 11)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 512;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 11: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;

		CalcTime t;

		Mat_32F x_32f(row, col);
		mat_rand(x_32f, 1.f, 255.f);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);

		//四則演算
		for (int k = 0; k < loop; k++)
		{
			//加算
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX ret_32f.data[i] =
			}
			t.end();
		}
		std::cout << "|method |time [ms]|" << std::endl;
		std::cout << "|-------|---------|" << std::endl;
		std::cout << "|add    |" << t.getAvgTime() << "|" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//乗算
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|mul    |" << t.getAvgTime() << "|" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//除算
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|div    |" << t.getAvgTime() << "|" << std::endl;

		//sqrt
		for (int k = 0; k < loop; k++)
		{
			//sqrt関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|sqrt   |" << t.getAvgTime() << "|" << std::endl;

		//sin
		for (int k = 0; k < loop; k++)
		{
			//sin関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|sin    |" << t.getAvgTime() << "|" << std::endl;

		//cos
		for (int k = 0; k < loop; k++)
		{
			//cos関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|cos    |" << t.getAvgTime() << "|" << std::endl;

		//exp
		for (int k = 0; k < loop; k++)
		{
			//exp関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|exp    |" << t.getAvgTime() << "|" << std::endl;

		//log
		for (int k = 0; k < loop; k++)
		{
			//log関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|log    |" << t.getAvgTime() << "|" << std::endl;


		float LUT[256];

		//sqrt LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			LUT[i] = (float)sqrt(i);// 以下の課題用のヒントのために，埋めてあります．
		}
		for (int k = 0; k < loop; k++)
		{
			//sqrt LUT
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				ret_32f.data[i] = LUT[(int)x_32f.data[i]];// 以下の課題用のヒントのために，埋めてあります．
			}
			t.end();
		}
		std::cout << "|sqrtLUT|" << t.getAvgTime() << "|" << std::endl;

		//sin LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXXXXXX
		}
		for (int k = 0; k < loop; k++)
		{
			t.start();
			//sin LUT
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|sin LUT|" << t.getAvgTime() << "|" << std::endl;

		//cos LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXXXXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//cos LUT
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|cos LUT|" << t.getAvgTime() << "|" << std::endl;

		//exp LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXXXXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//exp LUT
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|exp LUT|" << t.getAvgTime() << "|" << std::endl;

		//log LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXXXXXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//log LUT
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|log LUT|" << t.getAvgTime() << "|" << std::endl;


		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX115
	//課題(省略可)
	//ループ構造変換の課題にあるので省略した
	if (exercise == 115)
	{
		std::cout << "exercise 11.5(省略可)" << std::endl;
		CalcTime t;
		const int loop = 1000;
		const int size = 1024;
		Mat_32F a(size, size);
		Mat_32F b(size, size);

		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < size; j++)
			{
				for (int i = 0; i < size; i++)
				{
					if (i == j)
						a.data[size * j + i] = 1.f;
					else
						a.data[size * j + i] = b.data[size * j + i];
				}
			}
			t.end();
		}
		std::cout << "simple loop: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < size; j++)
			{
				for (int i = 0; i < size; i++)
				{
					a.data[size * j + i] = b.data[size * j + i];
				}
			}
			for (int i = 0; i < size; i++)
			{
				a.data[size * i + i] = 1.f;
			}
			t.end();
		}

		std::cout << "loop unswitching: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		return 0;
	}
#endif
#ifdef EX12
	//課題12
	//小さな行列A,Bの各要素の組(a,b)を回転させてx,yに格納する処理を，次の4方式で記述し，計算速度を比較せよ．
	//(1) noinline指定の関数呼出し（記述済み），(2) inline指定のない通常の関数呼出し，(3) inline指定の関数呼出し，(4) べた書き
	//関数はrot_noinline，rot_function，rot_inlineで，同じ処理（回転をROT_STAGES回繰り返す）を実装している．
	//cos，sinは関数の外で一度だけ計算して引数で渡す．
	//-O0，-O2，-O3で計測し，課題4.5の方法でアセンブリを観察して，関数呼出しが残っているか，要素ループがSIMD化されているかを確認せよ．
	//ROT_STAGESを1に変えて再ビルドし，8の場合と比較せよ．
	if (exercise == 12)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 128;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 12: loop = " << loop << ", size = " << size << ", stages = " << ROT_STAGES << std::endl << std::endl;

		const int row = size;
		const int col = size;

		CalcTime t;

		Mat_64F a(row, col);
		mat_rand(a, 0.0, 100.0);
		Mat_64F b(row, col);
		mat_rand(b, 0.0, 100.0);
		Mat_64F ansx(row, col);
		Mat_64F ansy(row, col);
		Mat_64F x(row, col);
		Mat_64F y(row, col);

		const double radian = 2.2;
		//cos，sinは計測の外で一度だけ計算する．
		const double c = cos(radian);
		const double s = sin(radian);

		std::cout << "|method   |time [ms]|" << std::endl;
		std::cout << "|---------|---------|" << std::endl;

		//(1) noinline指定の関数呼出し（記述済み．結果を正解としてansx，ansyに格納）
		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = a.cols * a.rows;
			for (int i = 0; i < size; i++)
			{
				rot_noinline(a.data[i], b.data[i], ansx.data[i], ansy.data[i], c, s);
			}
			t.end();
		}
		std::cout << "|noinline |" << t.getAvgTime() << "|" << std::endl;

		//(2) inline指定のない通常の関数呼出し
		mat_zero(x);//前の方式の結果を消す（計測外）
		mat_zero(y);
		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = a.cols * a.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX call rot_function
			}
			t.end();
		}
		std::cout << "|function |" << t.getAvgTime() << "|" << std::endl;
		const double diff_function = mat_diff(ansx, x) + mat_diff(ansy, y);

		//(3) inline指定の関数呼出し
		mat_zero(x);//前の方式の結果を消す（計測外）
		mat_zero(y);
		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = a.cols * a.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX call rot_inline
			}
			t.end();
		}
		std::cout << "|inline   |" << t.getAvgTime() << "|" << std::endl;
		const double diff_inline = mat_diff(ansx, x) + mat_diff(ansy, y);

		//(4) べた書き（関数の中身をループ内に書く）
		mat_zero(x);//前の方式の結果を消す（計測外）
		mat_zero(y);
		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = a.cols * a.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX rot_functionと同じ計算を書き，結果をx.data[i]，y.data[i]に格納する
			}
			t.end();
		}
		std::cout << "|hardcode |" << t.getAvgTime() << "|" << std::endl;
		const double diff_hardcode = mat_diff(ansx, x) + mat_diff(ansy, y);

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		std::cout << "diff from noinline: function = " << diff_function << ", inline = " << diff_inline << ", hardcode = " << diff_hardcode << std::endl;
		return 0;
	}
#endif
#ifdef EX13
	//課題13
	//行列A，Bの各要素の乗算を行うときに，結果を行列Cに格納する場合と行列Aに上書きする場合との計算時間をせよ．
	if (exercise == 13)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 1024;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 13: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;

		//c = a x b
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				//C=A*B
				Mat_64F a(row, col);
				Mat_64F b(row, col);
				Mat_64F c(row, col);
				const int size = a.rows * a.cols;
				for (int i = 0; i < size; i++)
				{
					//XXXXXXXX
				}
				t.end();
			}
			std::cout << "|method |time [ms]|" << std::endl;
			std::cout << "|-------|---------|" << std::endl;
			std::cout << "|alloc  |" << t.getAvgTime() << "|" << std::endl;
		}

		//a = a x b
		{

			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				//A=AxB
				Mat_64F a(row, col);
				Mat_64F b(row, col);
				const int size = a.rows * a.cols;
				for (int i = 0; i < size; i++)
				{
					//XXXXXXXX
				}
				t.end();
			}
			std::cout << "|inplace|" << t.getAvgTime() << "|" << std::endl;
		}

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX14
	//課題14
	//上記の0で初期化するコードをループの順序を変えてどちらが速いか計測して検証せよ．
	//また，行列積のコードのループの順序を変えてどれが速いか計測して検証せよ．
	//なお，行列積のコードは記述済みであり，どれも違いがないことを正解との差分の二乗和（MSE）で評価している．
	//行列のサイズを64，128，256，512，．．．と大きくしていくと，ループ順序の効果が大きいことがよくわかる．
	//その時，ループ回数を調整して少ないループにしないと計算がなかなか終わらない．（差が大きいため少ないループ数でも十分効果がわかる．）
	if (exercise == 14)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 128;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 14: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int row = size;
		const int col = size;

		CalcTime t;

		Mat_32F x(row, col);

		//col, row
		for (int k = 0; k < loop; k++)
		{
			int i = 0, j = 0;
			t.start();
			for (int i = 0; i < row; ++i)
			{
				for (int j = 0; j < col; ++j)
				{
					x.data[col * i + j] = 0.f;
				}
			}
			t.end();
		}
		std::cout << "|method |time [ms]|" << std::endl;
		std::cout << "|-------|---------|" << std::endl;
		std::cout << "|col-row|" << t.getAvgTime() << "|" << std::endl;

		//row, col
		for (int k = 0; k < loop; k++)
		{
			int i = 0, j = 0;
			t.start();
			for (int j = 0; j < col; ++j)
			{
				for (int i = 0; i < row; ++i)
				{
					x.data[col * i + j] = 0.f;
				}
			}
			t.end();
		}
		std::cout << "|row-col|" << t.getAvgTime() << "|" << std::endl;
		std::cout << std::endl << std::endl;


		//行列積
		Mat_32F a(size, size);
		Mat_32F b(size, size);
		Mat_32F c(size, size);
		Mat_32F ans(size, size);
		mat_rand(a, 0, 10);
		mat_rand(b, 0, 10);

		//for anser
		mat_zero(ans);
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				for (int k = 0; k < size; ++k)
				{
					//ans[i][j] += a[i][k] * b[k][j];
					ans.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
				}
			}
		}

		//i, j, k
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int i = 0; i < size; ++i)
			{
				float* pc = c.data + i * size;
				float* pa = a.data + i * size;
				for (int j = 0; j < size; ++j)
				{
					float* pb = b.data + j;
					float* pcj = &pc[j];
					for (int k = 0; k < size; ++k)
					{
						//c[i][j] += a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						*pcj += pa[k] * pb[k * size];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "|method|time [ms]|" << std::endl;
		std::cout << "|------|---------|" << std::endl;
		std::cout << "|i-j-k |" << t.getAvgTime() << "|" << std::endl;

		if (mat_diff(ans, c) > 1)std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//i, k, j
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int i = 0; i < size; ++i)
			{
				float* pc = c.data + i * size;
				float* pa = a.data + i * size;
				for (int k = 0; k < size; ++k)
				{
					float* pb = b.data + k * size;
					const float pak = pa[k];
					for (int j = 0; j < size; ++j)
					{
						//c[i][j] += a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						//pc[j] += (pak * pb[j]);//コンパイラ最適化でfmaがかかる
						pc[j] += pa[k] * pb[j];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "|i-k-j |" << t.getAvgTime() << "|" << std::endl;

		if (mat_diff(ans, c) > 1)std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//j, i, k
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int j = 0; j < size; ++j)
			{
				float* pb = b.data + j;
				for (int i = 0; i < size; ++i)
				{
					float* pc = c.data + i * size;
					float* pa = a.data + i * size;
					for (int k = 0; k < size; ++k)
					{
						//c[i][j] = c[i][j] + a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						pc[j] += pa[k] * pb[k * size];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "|j-i-k |" << t.getAvgTime() << "|" << std::endl;

		if (mat_diff(ans, c) > 1)std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//j, k, i
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int j = 0; j < size; ++j)
			{
				float* pc = c.data + j;
				float* pb = b.data + j;
				for (int k = 0; k < size; ++k)
				{
					float* pa = a.data + k;
					const float pbj = pb[k * size];
					for (int i = 0; i < size; ++i)
					{
						//c[i][j] = c[i][j] + a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						pc[i * size] += pa[i * size] * pbj;
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "|j-k-i |" << t.getAvgTime() << "|" << std::endl;
		if (mat_diff(ans, c) > 1)std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//k, i, j
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int k = 0; k < size; ++k)
			{
				float* pb = b.data + k * size;
				for (int i = 0; i < size; ++i)
				{
					float* pc = c.data + i * size;
					float* pa = a.data + i * size;
					const float pak = pa[k];
					for (int j = 0; j < size; ++j)
					{
						//c[i][j] = c[i][j] + a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						//pc[j] += pak * pb[j];//コンパイル最適化でFMAがかかる
						pc[j] += pa[k] * pb[j];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "|k-i-j |" << t.getAvgTime() << "|" << std::endl;

		if (mat_diff(ans, c) > 1)std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//k, j, i
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int k = 0; k < size; ++k)
			{
				float* pa = a.data + k;
				float* pb = b.data + k * size;
				for (int j = 0; j < size; ++j)
				{
					float* pc = c.data + j;
					const float pbj = pb[j];
					for (int i = 0; i < size; ++i)
					{
						//c[i][j] = c[i][j] + a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						pc[i * size] += pa[i * size] * pbj;
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "|k-j-i |" << t.getAvgTime() << "|" << std::endl;

		if (mat_diff(ans, c) > 1)std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX15
	//課題15
	//アンローリングの段数を2,4,8,16,32,...と変更することで，速度がどのように変わるか計測せよ．
	if (exercise == 15)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 65536;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 15: loop = " << loop << ", size = " << size << std::endl << std::endl;

		CalcTime t;

		//行列の2重ループだとコードが大変なので1重の配列をnew
		float* x = new float[size];
		float* y = new float[size];
		const float a = 2.f;
		const float b = 1.f;

		for (int i = 0; i < size; i++)
		{
			x[i] = rand_32f(0.f, 100.f);
		}

		//unrolling 1
		for (int j = 0; j < loop; j++)
		{
			//unrolling 1
			t.start();
			for (int i = 0; i < size; i++)
			{
				y[i + 0] = a * x[i + 0] + b;
			}
			t.end();
		}
		std::cout << "|method   |time [ms]|" << std::endl;
		std::cout << "|---------|---------|" << std::endl;
		std::cout << "|no unroll|" << t.getAvgTime() << "|" << std::endl;

		//unrolling 2
		for (int j = 0; j < loop; j++)
		{
			//unrolling 2
			t.start();
			for (int i = 0; i < size; i += 2)
			{
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|unroll  2|" << t.getAvgTime() << "|" << std::endl;

		//unrolling 4
		for (int j = 0; j < loop; j++)
		{
			//unrolling 4
			t.start();
			for (int i = 0; i < size; i += 4)
			{
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|unroll  4|" << t.getAvgTime() << "|" << std::endl;

		//unrolling 8
		for (int j = 0; j < loop; j++)
		{
			//unrolling 8
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|unroll  8|" << t.getAvgTime() << "|" << std::endl;

		//unrolling 16
		for (int j = 0; j < loop; j++)
		{
			//unrolling 16
			t.start();
			for (int i = 0; i < size; i += 16)
			{
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|unroll 16|" << t.getAvgTime() << "|" << std::endl;

		//unrolling 32
		for (int j = 0; j < loop; j++)
		{
			//unrolling 32
			t.start();
			for (int i = 0; i < size; i += 32)
			{
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|unroll 32|" << t.getAvgTime() << "|" << std::endl;

		//unrolling 64
		for (int j = 0; j < loop; j++)
		{
			//unrolling 64
			t.start();
			for (int i = 0; i < size; i += 64)
			{
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|unroll 64|" << t.getAvgTime() << "|" << std::endl;


		delete[] x;
		delete[] y;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX16
	//課題16
	//上記のプログラムを実装し，ループピーリングの有無で速度がどのように変わるか計測せよ．
	if (exercise == 16)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 65535;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 16: loop = " << loop << ", size = " << size << std::endl << std::endl;

		float* x = new float[size];
		float* y = new float[size];

		CalcTime t;

		//init
		for (int i = 0; i < size; i++)
		{
			x[i] = rand_32f(0, 100);
			y[i] = rand_32f(0, 100);
		}

		//Original
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; ++i)
			{
				if (i == 0)
				{
					y[i] = (x[i] + x[i + 1]) / 2.f;
				}
				else if (i == size - 1)
				{
					y[i] = (x[i - 1] + x[i]) / 2.f;
				}
				else
				{
					y[i] = (x[i - 1] + x[i] + x[i + 1]) / 3.f;
				}
			}
			t.end();
		}
		std::cout << "|method   |time [ms]|" << std::endl;
		std::cout << "|---------|---------|" << std::endl;
		std::cout << "|original |" << t.getAvgTime() << "|" << std::endl;

		//Loop peeling
		for (int j = 0; j < loop; j++)
		{
			t.start();
			{
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|L-peeling|" << t.getAvgTime() << "|" << std::endl;


		delete[] x;
		delete[] y;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX17
	//課題17
	//Mat_32Sの配同士の加算(a=a+b)をループつぶしをするか否かで計算時間を比較せよ．
	if (exercise == 17)
	{
		//空欄埋め問題
		const int default_loop = 1000;
		const int default_size = 256;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 17: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int width = size;
		const int height = size;
		//heightとwidthを大きくする場合は，注意すること（今はコメントアウト済み）．
		//下記にa[height][width]という宣言があり，大きすぎる配列はヒープ領域からメモリを取ってこれないため．

		CalcTime t;

		Mat_32S data(height, width);
		mat_rand(data, 0, 255);

		Mat_32S ma_ans(data);
		Mat_32S mb(data);

		//original
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < height; j++)
			{
				for (int i = 0; i < width; i++)
				{
					ma_ans.data[j * width + i] = ma_ans.data[j * width + i] + mb.data[j * width + i];
				}
			}
			t.end();
		}
		std::cout << "|method    |time [ms]|" << std::endl;
		std::cout << "|----------|---------|" << std::endl;
		std::cout << "|original  |" << t.getAvgTime() << "|" << std::endl;

		Mat_32S ma(data);
		//loop collapse
		//下記の2次元配列の場合のサンプルを参考にするとよい．
		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = width * height;
			/*//XXXXXXXX// hint: int* pa = ma.data;
			//XXXXXXXX
			for (int i = 0; i < size; i++)
			{
				//XXXXXXXX
			}*/

			t.end();
		}
		std::cout << "|L-collapse|" << t.getAvgTime() << "|" << std::endl;

		/*
		//2次元配列の場合のsample code
		//以下はコード記述済み．コンパイル時最適化のせいで，おそらく違いはない．興味がある人はgcc main.cpp -sでアセンブラ出力して確認すること．
		int a[height][width];
		int b[height][width];
		// before
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < height; j++)
			{
				for (int i = 0; i < width; i++)
				{
					a[j][i] = a[j][i] + b[j][i];
				}
			}
			t.end();
		}
		std::cout << "array before: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//after
		const int size = width * height;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			int* pa = &a[0][0];
			int* pb = &b[0][0];
			for (int i = 0; i < size; i++)
			{
				*pa++ = *pa + *pb++;
			}
			t.end();
		}
		std::cout << "array after : time (avg): " << t.getAvgTime() << " ms" << std::endl;
		*/

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		std::cout << "diff from ans: " << mat_diff(ma_ans, ma) << std::endl;
		return 0;
	}
#endif
#ifdef EX18
	//課題18
	//上記のコードを実行し，並列に動作していることを確認せよ．
	//また，並列化を有効にする場合としない場合の計算時間を比較せよ．
	//目的は並列動作の確認である．表示はスレッド間の同期が必要なため，並列化しても必ずしも速くはならない．
	if (exercise == 18)
	{
		//空欄埋め問題
		const int default_loop = 64;
		const int loop = (arg_loop == 0) ? default_loop : arg_loop;

		std::cout << "exercise 18: loop = " << loop << std::endl;
		//XXXXXXXX hint: #pragma...
		for (int i = 0; i < 100; i++)
		{
			std::cout << i << std::endl; //並列化したい処理
		}

		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = 0" << std::endl;
		return 0;
	}
#endif
#ifdef EX19
	//課題19
	//総和を計算するコードで，reduction指定子を使用する場合としない場合で計算結果がどのようになるか比較せよ．
	//なお，違いを確認するには，loopサイズは10回程度である．
	if (exercise == 19)
	{
		//空欄埋め問題
		const int default_loop = 10;
		const int default_size = 4096;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 19: loop = " << loop << ", size = " << size << std::endl << std::endl;

		std::cout << "|method       |";
		for (int j = 0; j < loop; j++)
		{
			std::cout << "      " << j << "|";
		}
		std::cout << std::endl;
		std::cout << "|-------------|";
		for (int j = 0; j < loop; j++)
		{
			std::cout << "-------|";
		}
		std::cout << std::endl;
		//並列化なし（逐次実行：全て記入済み


		std::cout << "|serial       |";
		for (int j = 0; j < loop; j++)
		{
			int sum = 0;
			for (int i = 0; i < size; i++)
			{
				sum += i;
			}
			std::cout << sum << "|";
		}
		std::cout << std::endl;

		//omp使用：空欄埋め
		std::cout << "|parallel     |";
		for (int j = 0; j < loop; j++)
		{
			int sum = 0;
			//XXXXXXXX
			for (int i = 0; i < size; i++)
			{
				sum += i;
			}
			std::cout << sum << "|";
		}
		std::cout << std::endl;

		//omp reduction使用：空欄埋め
		std::cout << "|parallel red.|";
		for (int j = 0; j < loop; j++)
		{
			int sum = 0;
			//XXXXXXXX hint reduction指定
			for (int i = 0; i < size; i++)
			{
				sum += i;
			}
			std::cout << sum << "|";
		}
		std::cout << std::endl;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX20
	//課題20（ローカルPCでの実行を前提とする．CSEサーバーで行う場合は講義時間外など利用者の少ない時間帯を選ぶこと．）
	//二つの行列の各要素の積を計算するコードで，スレッド数を変更して，計算時間がどのように推移するのかを確認せよ．
	//なお，スレッド数は，計算機のコア数以上の物まで指定せよ．
	//8コア16スレッドのPCでは，16コアよりも大きいスレッド数（例えば32までなど）までを指定せよ．
	if (exercise == 20)
	{
		//空欄埋め問題
		const int default_loop = 10;
		const int default_size = 512;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 20: loop = " << loop << ", size = " << size << std::endl << std::endl;

		CalcTime t;

		Mat_32F a(size, size);
		mat_rand(a, 0, 100);
		Mat_32F b(size, size);
		mat_rand(b, 0, 100);
		Mat_32F c(size, size);

		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			//#pragma omp parallel for num_threads(n)で並列化，nに任意の整数を入れる
			//XXXXXXXX
			for (int i = 0; i < size; ++i)
			{
				float* pc = c.data + i * size;
				float* pa = a.data + i * size;
				for (int k = 0; k < size; ++k)
				{
					float* pb = b.data + k * size;
					const float pak = pa[k];
					for (int j = 0; j < size; ++j)
					{
						pc[j] += (pak * pb[j]);
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "|num_threads(8)|time [ms]|" << std::endl;
		std::cout << "|--------------|---------|" << std::endl;
		std::cout << "|xx            |" << t.getAvgTime() << "|" << std::endl;

		std::cout << std::endl;

		//上記はプログラムコンパイル時にスレッド数を決定する方法．
		//プログラム実行時にスレッド数を変更するには下記のset_num_threadを用いる．
		//詳細は，並列化用の関数群 omp.h　の章を参照のこと．
		//上記は，何度もコンパイルしないといけないが，下記は一度だけでよいのでif (isUseSetNumThread)のコメントアウトなどをうまく使うこと．
		bool isUseSetNumThread = false;
		//if (isUseSetNumThread)
		{
			const int threadMax = 32;
			std::cout << "|set_thread|time [ms]|" << std::endl;
			std::cout << "|----------|---------|" << std::endl;
			for (int nt = 1; nt < threadMax; nt++)
			{
				omp_set_num_threads(nt);
				for (int l = 0; l < loop; l++)
				{
					mat_zero(c);
					t.start();
#pragma omp parallel for
					for (int i = 0; i < size; ++i)
					{
						float* pc = c.data + i * size;
						float* pa = a.data + i * size;
						for (int k = 0; k < size; ++k)
						{
							float* pb = b.data + k * size;
							const float pak = pa[k];
							for (int j = 0; j < size; ++j)
							{
								pc[j] += (pak * pb[j]);
							}
						}
					}
					t.end();
				}
				std::cout << "|set      " << nt << "|" << t.getAvgTime() << "|" << std::endl;
			}
		}

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX21
	//課題21
	//四則演算のコードを書いてprintfデバッグで確認せよ．
	if (exercise == 21)
	{
		std::cout << "exercise 21" << std::endl;
		const __m256 a = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
		const __m256 b = _mm256_set_ps(15, 14, 13, 12, 11, 10, 9, 8);
		//逆順で入力(setr)こちらのほうが使いやすい
		const __m256 c = _mm256_setr_ps(0, 1, 2, 3, 4, 5, 6, 7);
		const __m256 d = _mm256_setr_ps(8, 9, 10, 11, 12, 13, 14, 15);
		__m256 e = _mm256_setzero_ps();

		//加算
		//e = a + b
		//XXXXXXXX
		std::cout << "add a b: ";
		print_m256(e);
		//e = c + d
		//XXXXXXXX
		std::cout << "add c d: ";
		print_m256(e);//出力は上と同じのはず．

		//表示の仕方色々．
		//print_m256の関数を使わない場合はいったんstoreして出力
		//_mm256_store_psは32バイト境界にアラインされたアドレスが必要なため，alignas(32)を付ける．
		//付けない場合はsegmentation faultする可能性がある．アラインを気にしない場合は_mm256_storeu_psを使う．
		alignas(32) float temp[8];
		_mm256_store_ps(&temp[0], e);
		std::cout << "cout ex1: ";
		for (int i = 0; i < 8; i++) std::cout << temp[i] << ", ";
		std::cout << std::endl;
		//下記のように直接参照も可能
		std::cout << "cout ex2: ";
		for (int i = 0; i < 8; i++) std::cout << ((float*)&e)[i] << ", ";
		std::cout << std::endl;
		//visual studioなら更にこのような出力も可能（gccはコンパイルが通らないのでmsc_verのマクロでコメントアウト中．）
#ifdef _MSC_VER
		std::cout << "cout ex3: ";
		for (int i = 0; i < 8; i++) std::cout << e.m256_f32[i] << ", ";
		std::cout << std::endl;
#endif

		std::cout << std::endl;
		//減算
		//e = a - b
		//XXXXXXXX
		std::cout << "sub a b: ";
		print_m256(e);
		//e = c - d
		//XXXXXXXX
		std::cout << "sub c d: ";
		print_m256(e);

		std::cout << std::endl;
		//乗算
		//e = a * b
		//XXXXXXXX
		std::cout << "mul a b: ";
		print_m256(e);
		//e = c * d
		//XXXXXXXX
		std::cout << "mul c d: ";
		print_m256(e);

		std::cout << std::endl;
		//除算
		//e = a / b
		//XXXXXXXX
		std::cout << "div a b: ";
		print_m256(e);
		//e = c / d
		//XXXXXXXX
		std::cout << "div c d: ";
		print_m256(e);

		return 0;
	}
#endif
#ifdef EX22
	//課題22
	//(1) 配列a,x,bに対して，`(((a*x+b)*x+b)*x+b)*x+b `の計算を配列ｃに格納するコードをmul/addで記述するものとFMAを使うもので記述し，計算速度と計算結果を比較せよ．
	//なお，上記の関数は以下に等しい．
	//a=_mm256_fmadd_ps(a,b,c);
	//a=_mm256_fmadd_ps(a,b,c);
	//a=_mm256_fmadd_ps(a,b,c);
	//a=_mm256_fmadd_ps(a,b,c);
	//これは，単純にFMAが1度だとメモリで律速するこのコードでは計算速度の差が出にくいためである．差が小さければ，より演算を増やせば良い．
	//なお，現在のg++では，-O2以上でmul/addの命令もFMAに自動的に変換される．-ffp-contract=offを付けると変換を抑制できる．

	//(2) また，GFLOPSと演算強度[FLOPS / BYTE]をFMA命令で計算する関数`loofline_test`を使って，ルールラインのグラフとして図示せよ．
	if (exercise == 22)
	{
		//空欄埋め問題
		const int default_loop = 1000000;
		const int default_size = 32;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 22: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int col = size;
		const int row = size;

		CalcTime t;

		Mat_32F a(row, col);
		Mat_32F x(row, col);
		Mat_32F b(row, col);
		Mat_32F ans(row, col);
		Mat_32F c(row, col);
		mat_rand(a, 0, 1);
		mat_rand(b, 0, 1);
		mat_rand(x, 0, 1);
		mat_zero(ans);
		mat_zero(c);

		//mul, add
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* pa = a.data;
			float* px = x.data;
			float* pb = b.data;
			float* pc = ans.data;
			const int s = col * row;
			for (int i = 0; i < s; i += 8)
			{
				const __m256 ma = _mm256_load_ps(pa + i);
				const __m256 mx = _mm256_load_ps(px + i);
				const __m256 mb = _mm256_load_ps(pb + i);

				//mul,addを使って（結果はtempに入れること）
				__m256 temp = _mm256_setzero_ps();
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				_mm256_store_ps(pc + i, temp);
			}
			t.end();
		}
		std::cout << "|method |time [ms]|" << std::endl;
		std::cout << "|-------|---------|" << std::endl;
		std::cout << "|mul-add|" << t.getAvgTime() << "|" << std::endl;

		//fma
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* pa = a.data;
			float* px = x.data;
			float* pb = b.data;
			float* pc = c.data;
			const int s = col * row;
			for (int i = 0; i < s; i += 8)
			{
				const __m256 ma = _mm256_load_ps(pa + i);
				const __m256 mx = _mm256_load_ps(px + i);
				const __m256 mb = _mm256_load_ps(pb + i);

				//fmaを使って（結果はtempに入れること）
				__m256 temp = _mm256_setzero_ps();
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				_mm256_store_ps(pc + i, temp);
			}
			t.end();
		}
		std::cout << "|fma    |" << t.getAvgTime() << "|" << std::endl;

		std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;


		//loofline test
		//上記の実験がしたいときは，下記フラグを制御するとよい．
		const bool isUseLoopLineTest = false;
		//if (isUseLoopLineTest)
		{
			std::cout << "plot loofline" << std::endl;

			const int loofline_size = 16 * 1024 / sizeof(float);//16KByte
			const int iteration = 1000000;

			std::cout << "single core performance" << std::endl;
			loofline_test<loofline_size>(iteration, 1);
			std::cout << std::endl;
			std::cout << "multi core performance" << std::endl;
			loofline_test<loofline_size>(iteration);

			//ベクトル化の性能をオートベクタライゼーションに頼るときはこっち
			//std::cout << "single core performance" << std::endl;
			//loofline_test_cpp<loofline_size>(iteration, 1);
			//std::cout << std::endl;
			//std::cout << "multi core performance" << std::endl;
			//loofline_test_cpp<loofline_size>(iteration);
		}

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX23
	//課題23
	//divとrcp,sqrtとrsqrtの実行速度を比較せよ．
	//また，絶対値の計算をandによる実装とmaxによる実装を比較せよ．
	//なお，mulとrcpの組み合わせは，計算機によってはdivよりも遅い可能性があり，また，ビット演算による絶対値計算も，演算強度とレジスタ数の関係からmaxとsubの演算よりも遅い場合も遅い可能性がある．
	if (exercise == 23)
	{
		//空欄埋め問題
		const int default_loop = 1000000;
		const int default_size = 64;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 23: loop = " << loop << ", size = " << size << std::endl << std::endl;

		Mat_32F a(size, size);
		Mat_32F b(size, size);
		Mat_32F c(size, size);
		mat_rand(a, 1.f, 100.f);
		mat_rand(b, 1.f, 100.f);

		CalcTime t;
		//div
		const int matsize = size * size;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < matsize; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a.data + i);
				const __m256 mb = _mm256_load_ps(b.data + i);

				__m256 temp = _mm256_setzero_ps();
				//divを使って
				//XXXXXXXX

				_mm256_store_ps(c.data + i, temp);
			}
			t.end();
		}
		std::cout << "|method    |time [ms]|" << std::endl;
		std::cout << "|----------|---------|" << std::endl;
		std::cout << "|div       |" << t.getAvgTime() << "|" << std::endl;

		//rcp
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < matsize; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a.data + i);
				const __m256 mb = _mm256_load_ps(b.data + i);

				__m256 temp = _mm256_setzero_ps();
				//rcpとmulをつかって
				//XXXXXXXX

				_mm256_store_ps(c.data + i, temp);
			}
			t.end();
		}
		std::cout << "|rcp       |" << t.getAvgTime() << "|" << std::endl;

		//sqrt
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < matsize; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a.data + i);
				const __m256 mb = _mm256_load_ps(b.data + i);

				__m256 temp = _mm256_setzero_ps();
				//sqrtを使って
				//XXXXXXXX

				_mm256_store_ps(c.data + i, temp);
			}
			t.end();
		}
		std::cout << "|sqrt      |" << t.getAvgTime() << "|" << std::endl;

		//rsqrt+rcp
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < matsize; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a.data + i);
				const __m256 mb = _mm256_load_ps(b.data + i);

				__m256 temp = _mm256_setzero_ps();
				//rsqrtとrcpを使って
				//XXXXXXXX
				temp = _mm256_rcp_ps(_mm256_rsqrt_ps(ma));

				_mm256_store_ps(c.data + i, temp);
			}
			t.end();
		}
		std::cout << "|rsqrt+rcp |" << t.getAvgTime() << "|" << std::endl;

		//rsqrt+mul
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < matsize; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a.data + i);
				const __m256 mb = _mm256_load_ps(b.data + i);

				__m256 temp = _mm256_setzero_ps();
				//rsqrtとmulを使って（ルートの逆数は乗算で戻る）
				//XXXXXXXX

				_mm256_store_ps(c.data + i, temp);
			}
			t.end();
		}
		std::cout << "|rsqrt+mul |" << t.getAvgTime() << "|" << std::endl;

		//abs(subとmaxを使って)
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < matsize; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a.data + i);
				const __m256 mb = _mm256_load_ps(b.data + i);

				__m256 temp = _mm256_setzero_ps();
				//subとmaxを使って
				//XXXXXXXX

				_mm256_store_ps(c.data + i, temp);
			}
			t.end();
		}
		std::cout << "|abs submax|" << t.getAvgTime() << "|" << std::endl;

		//abs(subとandを使って)
		for (int j = 0; j < loop; j++)
		{
			t.start();
			const int v32f_absmask[] = { 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff };
			// __attribute__ ((aligned(32))) const int v32f_absmask[] = { 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff }; // segmentation faultする場合はこっち
			const __m256 absmask = *(const __m256*)(&v32f_absmask[0]);
			for (int i = 0; i < matsize; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a.data + i);
				const __m256 mb = _mm256_load_ps(b.data + i);

				__m256 temp = _mm256_setzero_ps();
				//subとandを使って．andのマスクはabsmask
				//XXXXXXXX

				_mm256_store_ps(c.data + i, temp);
			}
			t.end();
		}
		std::cout << "|abs suband|" << t.getAvgTime() << "|" << std::endl;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX24
	//課題24
	//haddとdpで要素の総和を取るプログラムを作成し，それぞれの計算時間を比較せよ．なお，大きな差はない．
	//なお，最も速い実装は，shuffleでデータを入れ替えてaddする方法であるが特に課題指定はない．興味があれば挑戦してみるとよい．
	//また，正解を計算結果を比較するために，非ベクトル化コードが実装されている．
	//（この課題は，後にリダクションの最適化でもう一度登場する．）
	if (exercise == 24)
	{
		//空欄埋め問題
		const int default_loop = 100000;
		const int default_size = 128;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 24: loop = " << loop << ", size = " << size << std::endl << std::endl;

		Mat_32F a(size, size);
		mat_rand(a, 0.f, 0.0001f);

		const int matsize = size * size;

		CalcTime t;

		//ベクトル化無し（自動ベクトル化されている可能性もある）
		//VisualStudioだと最適化で最外loopを消されて実験できないが，g++は今のところ動いている
		//全て埋めてある
		double ans_double = 0.0;
		for (int i = 0; i < loop; i++)
		{
			t.start();
			double sum = 0.0;
			for (int i = 0; i < matsize; i++)
			{
				sum += (double)a.data[i];
			}
			t.end();
			ans_double = sum;
		}
		std::cout << "|method|time [ms]|total  |" << std::endl;
		std::cout << "|------|---------|-------|" << std::endl;
		std::cout << "|double|" << t.getAvgTime() << "|" << ans_double << "|" << std::endl;
		float ans = 0.f;
		for (int i = 0; i < loop; i++)
		{
			t.start();
			float sum = 0.f;
			for (int i = 0; i < matsize; i++)
			{
				sum += a.data[i];
			}
			t.end();
			ans = sum;
		}
		std::cout << "|float |" << t.getAvgTime() << "|" << ans << "|" << std::endl;

		//hadd
		float ans_hadd = 0.f;
		for (int i = 0; i < loop; i++)
		{
			t.start();
			float sum = 0.f;
			for (int i = 0; i < matsize; i += 8)
			{
				__m256 ma = _mm256_load_ps(a.data + i);
				//haddを使って
				//XXXXXXXX
				//XXXXXXXX
				//hint sum+= XXXXXXXX
			}
			t.end();
			ans_hadd = sum;
		}
		std::cout << "|hadd  |" << t.getAvgTime() << "|" << ans_hadd << "|" << std::endl;
		if (abs(ans - ans_hadd) > 1)std::cout << "total ans and_hadd " << ans << "," << ans_hadd << std::endl;

		//dp
		float ans_dp = 0.f;
		const __m256 one = _mm256_set1_ps(1.f);
		for (int i = 0; i < loop; i++)
		{
			t.start();
			float sum = 0.f;
			for (int i = 0; i < matsize; i += 8)
			{
				__m256 ma = _mm256_load_ps(a.data + i);
				//dpを使って
				//XXXXXXXX
				//hint: sum+=XXXXXXXX
			}
			t.end();
			ans_dp = sum;
		}
		std::cout << "|dp    |" << t.getAvgTime() << "|" << ans_dp << "|" << std::endl;
		if (abs(ans - ans_dp) > 1)std::cout << "total ans and_dp " << ans << "," << ans_hadd << std::endl;


		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX25
	//課題25
	//行列aにおいて要素の値があるしきい値以上の場合だけ3乗し，それ以外は何もしない処理をベクトル化実装せよ．
	//if(a[i]>=threshold) a[i]=a[i]*a[i]*a[i];
	if (exercise == 25)
	{
		//空欄埋め問題
		const int default_loop = 100000;
		const int default_size = 128;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 25: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const float threshold = 4.f;

		const int matsize = size * size;
		Mat_32F a(size, size);
		Mat_32F ans(size, size);
		Mat_32F b(size, size);
		mat_rand(a, 0, 100);

		CalcTime t;

		//C++：埋めてある
		for (int i = 0; i < loop; i++)
		{
			t.start();
			for (int i = 0; i < matsize; i++)
			{
				if (a.data[i] >= threshold) ans.data[i] = a.data[i] * a.data[i] * a.data[i];
				else ans.data[i] = a.data[i];
			}
			t.end();
		}
		std::cout << "|method|time [ms]|" << std::endl;
		std::cout << "|------|---------|" << std::endl;
		std::cout << "|scalar|" << t.getAvgTime() << "|" << std::endl;


		for (int i = 0; i < loop; i++)
		{
			t.start();
			const __m256 mth = _mm256_set1_ps(threshold);
			for (int i = 0; i < matsize; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a.data + i);

				__m256 temp = _mm256_setzero_ps();
				//cmp, mul, blendvを使って（blendを使わずにビット演算でもできる）
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				_mm256_store_ps(b.data + i, temp);
			}
			t.end();
		}
		std::cout << "|vector|" << t.getAvgTime() << "|" << std::endl;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		std::cout << "diff: " << mat_diff(ans, b) << std::endl;

		return 0;
	}
#endif
#ifdef EX26
	//課題26
	//上記のコードのように，SIMD命令を使う場合におけるループアンローリングを8，16，32，64と行い，計算時間を比較せよ．
	//ただし，ベクトル化していないコードのほうが速い可能性が高い．これは，これくらい単純なコードは，自動ベクトル化によってコードが最適化されるため．
	if (exercise == 26)
	{
		//空欄埋め問題
		const int default_loop = 100000;
		const int default_size = 128;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;

		std::cout << "exercise 26: loop = " << loop << ", size = " << size << std::endl << std::endl;

		const int matsize = size * size;
		Mat_32F a(size, size);
		Mat_32F b(size, size);
		Mat_32F c(size, size);
		Mat_32F ans(size, size);
		mat_rand(a, 0, 100);
		mat_rand(b, 0, 100);
		mat_zero(c);
		mat_zero(ans);

		CalcTime t;

		// unrolling 1: 埋めてある．自動ベクトル化の可能性大
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 1
			for (int i = 0; i < matsize; i++)
			{
				ans.data[i] = (a.data[i] - b.data[i]) * (a.data[i] - b.data[i]);
			}
			t.end();
		}
		std::cout << "|method   |time [ms]|" << std::endl;
		std::cout << "|---------|---------|" << std::endl;
		std::cout << "|no unroll|" << t.getAvgTime() << "|" << std::endl;

		// unrolling 1: 埋めてある．強制的にベクトル化するのを排除している．_mm_sub_ssのssがそれ．
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 1
			for (int i = 0; i < matsize; i++)
			{
				__m128 ma = _mm_loadu_ps(a.data + i);
				__m128 mb = _mm_loadu_ps(b.data + i);
				__m128 temp = _mm_sub_ss(ma, mb);
				_mm_storeu_ps(c.data + i, _mm_mul_ss(temp, temp));
			}
			t.end();
		}
		std::cout << "|unroll  1|" << t.getAvgTime() << "|" << std::endl;


		// unrolling 8: 埋めてある
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 8
			for (int i = 0; i < matsize; i += 8)
			{
				__m256 ma = _mm256_load_ps(a.data + i);
				__m256 mb = _mm256_load_ps(b.data + i);
				__m256 temp = _mm256_sub_ps(ma, mb);
				_mm256_store_ps(c.data + i, _mm256_mul_ps(temp, temp));
			}
			t.end();
		}
		std::cout << "|unroll  8|" << t.getAvgTime() << "|" << std::endl;


		// unrolling 16
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 16
			for (int i = 0; i < matsize; i += 16)
			{
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|unroll 16|" << t.getAvgTime() << "|" << std::endl;
		if (mat_diff(ans, c) > 1)std::cout << "invalid: diff=" << mat_diff(ans, c) << std::endl;

		// unrolling 32
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 32
			for (int i = 0; i < matsize; i += 32)
			{
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|unroll 32|" << t.getAvgTime() << "|" << std::endl;
		if (mat_diff(ans, c) > 1)std::cout << "invalid: diff=" << mat_diff(ans, c) << std::endl;

		// unrolling 64
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 64
			for (int i = 0; i < matsize; i += 64)
			{
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX

				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
				//XXXXXXXX
			}
			t.end();
		}
		std::cout << "|unroll 64|" << t.getAvgTime() << "|" << std::endl;
		if (mat_diff(ans, c) > 1)std::cout << "invalid: diff=" << mat_diff(ans, c) << std::endl;

		// unrolling 128
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 128
			for (int i = 0; i < matsize; i += 128)
			{
				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "|unroll128|" << t.getAvgTime() << "|" << std::endl;
		if (mat_diff(ans, c) > 1)std::cout << "invalid: diff=" << mat_diff(ans, c) << std::endl;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		return 0;
	}
#endif
#ifdef EX27
	//課題27
	//上の関数を用いた8x8のfloatの転置の動作を確認せよ．また，4x4のdouble型データの転置を作成せよ．
	//難しいと思った場合は，効率を無視して，set命令やstoreしてスカラで書くなどすれば書ける．
	if (exercise == 27)
	{
		std::cout << "exercise 27" << std::endl;
		{
			std::cout << "float 8x8" << std::endl;
			const int size = 64;
			float* a = (float*)_mm_malloc(sizeof(float) * size, 32);
			float* b = (float*)_mm_malloc(sizeof(float) * size, 32);

			for (int i = 0; i < size; i++)
			{
				a[i] = (float)i;
				b[i] = 0.f;
			}
			__m256 ma[8], mb[8];
			for (int i = 0; i < 8; i++)
			{
				ma[i] = _mm256_load_ps((float*)(&a[i * 8]));
				mb[i] = _mm256_setzero_ps();
			}

			for (int i = 0; i < 8; i++)
			{
				print_m256(ma[i]);
			}

			//転置
			_mm256_transpose_8x8_ps(mb, ma);

			std::cout << "transpose" << std::endl;
			for (int i = 0; i < 8; i++)
			{
				print_m256(mb[i]);
			}
			_mm_free(a);
			_mm_free(b);
		}
		{
			std::cout << "double 4x4" << std::endl;
			const int size = 16;
			double* a = (double*)_mm_malloc(sizeof(double) * size, 32);
			double* b = (double*)_mm_malloc(sizeof(double) * size, 32);

			for (int i = 0; i < size; i++)
			{
				a[i] = i;
				b[i] = 0;
			}
			__m256d ma[4], mb[4];
			for (int i = 0; i < 4; i++)
			{
				ma[i] = _mm256_load_pd((double*)(&a[i * 4]));
				mb[i] = _mm256_setzero_pd();
			}

			for (int i = 0; i < 4; i++)
			{
				print_m256d(ma[i]);
			}

			//転置
			//作成する（行数は適当）
			//XXXXXXXX
			//XXXXXXXX
			//XXXXXXXX
			//XXXXXXXX

			//転置結果の表示
			std::cout << "transpose" << std::endl;
			for (int i = 0; i < 4; i++)
			{
				print_m256d(mb[i]);
			}
			_mm_free(a);
			_mm_free(b);
		}
		return 0;
	}
#endif
#ifdef EX28
	//課題28
	//__m256i（int）型を_m256（float）型に変換せよ．
	//また，unsigned char型をfloat型に変換せよ．
	//更に，16個の`int`型を`short`型に変換する処理をSSEとAVXで実装せよ．
	if (exercise == 28)
	{
		std::cout << "exercise 28" << std::endl;

		//int ->float
		std::cout << "int ->float" << std::endl;
		__m256i m32i = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
		__m256 m32f = _mm256_setzero_ps();
		std::cout << "input data: " << std::endl;
		print_m256i_i32(m32i);
		std::cout << "before convert" << std::endl;
		print_m256(m32f);

		//cvtを使って
		//XXXXXXXX

		std::cout << "after convert" << std::endl;
		print_m256(m32f);

		//unsigne char->float
		unsigned char* a = (unsigned char*)_mm_malloc(sizeof(char) * 32, 32);
		float* b = (float*)_mm_malloc(sizeof(float) * 32, 32);
		for (int i = 0; i < 32; i++)
		{
			a[i] = i;
			b[i] = 0.f;
		}

		std::cout << std::endl << "unsigne char->float" << std::endl;
		std::cout << "input data: " << std::endl;
		for (int i = 0; i < 32; i++)std::cout << (float)a[i] << " ";
		std::cout << std::endl;

		std::cout << "before convert: b" << std::endl;
		for (int i = 0; i < 32; i++)std::cout << b[i] << " ";
		std::cout << std::endl;

		//aをfloatにキャストしてbへ書き込み
		for (int i = 0; i < 32; i += 8)
		{
			//XXXXXXXX
			//XXXXXXXX
			//XXXXXXXX
		}

		std::cout << "after convert: b" << std::endl;
		for (int i = 0; i < 32; i++)std::cout << b[i] << " ";
		std::cout << std::endl;

		//int->short
		std::cout << std::endl << "int->short" << std::endl;
		int* c = (int*)_mm_malloc(sizeof(int) * 16, 32);
		short* d = (short*)_mm_malloc(sizeof(short) * 16, 32);

		for (int i = 0; i < 16; i++)
		{
			c[i] = i;
			d[i] = 0;
		}
		std::cout << "input data: " << std::endl;
		for (int i = 0; i < 16; i++)std::cout << c[i] << " ";
		std::cout << std::endl;

		std::cout << "before convert: d" << std::endl;
		for (int i = 0; i < 16; i++)std::cout << d[i] << " ";
		std::cout << std::endl;

		//SSEでの実装 cを入力として，dに書き込み
		//ヒント：前半8個と後半8個に分けて2回処理する．packs_epi16で半分のサイズへ変換できる
		//ヒント：下記のAVXの作りかけもヒントになる．
		//XXXX 行数は任意
		//XXXX
		//XXXX

		//結果の表示
		std::cout << "after convert: d (SSE)" << std::endl;
		for (int i = 0; i < 16; i++)std::cout << d[i] << " ";
		std::cout << std::endl;

		//AVXでの実装
		__m256i mc2560 = _mm256_load_si256((__m256i*)c);
		__m256i mc2561 = _mm256_load_si256((__m256i*)(c + 8));
		__m256i temp256 = _mm256_packs_epi16(mc2560, mc2561);
		//permuteをしていないため結果がおかしい
		//XXXXXXXX temp256をpermute

		_mm256_store_si256((__m256i*)d, temp256);

		std::cout << "after convert: d (AVX)" << std::endl;
		for (int i = 0; i < 16; i++)std::cout << d[i] << " ";
		std::cout << std::endl;

		_mm_free(a);
		_mm_free(b);
		_mm_free(c);
		_mm_free(d);

		return 0;
	}
#endif
#ifdef EX29
	//課題29
	//上記のコードのスカラ実装，スカラ実装＋並列化，SIMDのみ，SIMD＋並列化を作成し，計算時間を比較せよ．
	//
	//課題25のように行列aにおいて要素の値があるしきい値以上の場合だけ3乗し，それ以外は何もしない処理をスカラ実装，スカラ実装＋並列化，SIMDのみ，SIMD＋並列化で作成し，計算時間を比較せよ．
	if (exercise == 29)
	{
		//空欄埋め問題
		const int default_loop = 10000;
		const int default_size = 256;

		const int loop = (arg_loop == 0) ? default_loop : arg_loop;
		const int size = (arg_size == 0) ? default_size : arg_size;
		if (size % 8 != 0)
		{
			std::cout << "input size must be 8 multiple." << std::endl;
		}
		std::cout << "exercise 29: loop = " << loop << ", size = " << size << std::endl << std::endl;

		Mat_32F a(size, size);
		Mat_32F b(size, size);
		Mat_32F ans(size, size);

		//init
		mat_rand(a, 0, 100);

		const float threshold = 20.f;

		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			//スカラー実装
			t.start();
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i++)
				{
					//XXXXXXXX ans.data[j * size + i] = ここだけ，ansに書き込み
				}
			}
			t.end();
		}
		std::cout << "|method    |time [ms]|" << std::endl;
		std::cout << "|----------|---------|" << std::endl;
		std::cout << "|scalar    |" << t.getAvgTime() << "|" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//スカラー，並列化実装
			t.start();
			//XXXXXXXX
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i++)
				{
					//XXXXXXXX　bに書き込み
				}
			}
			t.end();
		}
		std::cout << "|scalar+omp|" << t.getAvgTime() << "|" << std::endl;
		if (mat_diff(ans, b) > 1)std::cout << "invalid: diff = " << mat_diff(ans, b) << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//SIMD実装
			t.start();
			//XXXXXXXX（閾値のセット）
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i += 8)
				{
					//XXXXXXXX
					//XXXXXXXX
					//XXXXXXXX
					//XXXXXXXX
				}
			}
			t.end();
		}
		std::cout << "|SIMD      |" << t.getAvgTime() << "|" << std::endl;
		if (mat_diff(ans, b) > 1)std::cout << "invalid: diff = " << mat_diff(ans, b) << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//SIMD，並列化実装
			t.start();
			//XXXXXXXX
			for (int j = 0; j < size; ++j)
			{
				//XXXXXXXX（閾値のセット）	
				for (int i = 0; i < size; i += 8)
				{
					//XXXXXXXX
					//XXXXXXXX
					//XXXXXXXX
					//XXXXXXXX
				}
			}
			t.end();
		}
		std::cout << "|SIMD+omp  |" << t.getAvgTime() << "|" << std::endl;
		if (mat_diff(ans, b) > 1)std::cout << "invalid: diff = " << mat_diff(ans, b) << std::endl;

		std::cout << std::endl << "info:" << std::endl;
		std::cout << "default parameter: default_loop = " << default_loop << ", default_size = " << default_size << std::endl;
		//std::cout << "diff: " << mat_diff(ans, b) << std::endl;
		return 0;
	}

#endif

	std::cout << "no select" << std::endl;
	return 0;
}

#ifdef EX2
//課題2用
void GEMM(Mat_32F& a, Mat_32F& b, Mat_32F& c)
{
	const int size = a.cols;
#pragma omp parallel for
	for (int i = 0; i < size; ++i)
	{
		float* pc = c.data + i * size;
		float* pa = a.data + i * size;
		for (int k = 0; k < size; ++k)
		{
			float* pb = b.data + k * size;
			const float pak = pa[k];
			for (int j = 0; j < size; ++j)
			{
				pc[j] += (pak * pb[j]);
			}
		}
	}
}

#endif
#ifdef EX12
//課題12用
//3つの関数は同じ処理を行う．cos(radian)，sin(radian)をc，sとして受け取り，(a,b)の回転をROT_STAGES回繰り返す．
NOINLINE void rot_noinline(double a, double b, double& x, double& y, double c, double s)
{
	for (int j = 0; j < ROT_STAGES; j++)
	{
		const double t = a * c - b * s;
		b = a * s + b * c;
		a = t;
	}
	x = a;
	y = b;
}

void rot_function(double a, double b, double& x, double& y, double c, double s)
{
	for (int j = 0; j < ROT_STAGES; j++)
	{
		const double t = a * c - b * s;
		b = a * s + b * c;
		a = t;
	}
	x = a;
	y = b;
}

inline void rot_inline(double a, double b, double& x, double& y, double c, double s)
{
	for (int j = 0; j < ROT_STAGES; j++)
	{
		const double t = a * c - b * s;
		b = a * s + b * c;
		a = t;
	}
	x = a;
	y = b;
}
#endif
#ifdef EX27
//課題27用
inline void _mm256_transpose_8x8_ps(__m256* dst, const __m256* src)
{
	__m256  tmp[8], tmpp[8];

	for (int i = 0; i < 8; i += 2)
	{
		tmp[i + 0] = _mm256_unpacklo_ps(src[i], src[i + 1]);
		tmp[i + 1] = _mm256_unpackhi_ps(src[i], src[i + 1]);
	}
	for (int i = 0; i < 8; i += 4)
	{
		tmpp[i + 0] = _mm256_shuffle_ps(tmp[i + 0], tmp[i + 2], _MM_SHUFFLE(1, 0, 1, 0));
		tmpp[i + 1] = _mm256_shuffle_ps(tmp[i + 0], tmp[i + 2], _MM_SHUFFLE(3, 2, 3, 2));
	}
	for (int i = 0; i < 8; i += 4)
	{
		tmpp[i + 2] = _mm256_shuffle_ps(tmp[i + 1], tmp[i + 3], _MM_SHUFFLE(1, 0, 1, 0));
		tmpp[i + 3] = _mm256_shuffle_ps(tmp[i + 1], tmp[i + 3], _MM_SHUFFLE(3, 2, 3, 2));
	}
	for (int i = 0; i < 4; i++)
	{
		dst[i + 0] = _mm256_permute2f128_ps(tmpp[i], tmpp[i + 4], 0x20);
		dst[i + 4] = _mm256_permute2f128_ps(tmpp[i], tmpp[i + 4], 0x31);
	}
}
#endif
