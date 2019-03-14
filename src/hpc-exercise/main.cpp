#include "utils/mat.h"
#include "utils/mat_util.h"
#include "utils/simd_util.h"
#include <math.h>

//課題12
inline void rot(double a, double b, double &x, double &y, double radian)
{
	x = a * cos(radian);
	y = b * sin(radian);
}

//課題27
void _mm256_transpose_8x8_ps(__m256* dst, const __m256* src)
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

int main(const int argc, const char** argv)
{
	//課題1
	//行列積和演算AX+Bを計算するプログラムにおいて，行列積と和それぞれの実行時間をタイマーを挟むことで測定せよ．
	if (false)
	{
		std::cout << "課題1" << std::endl;
		const int loop = 10;
		const int row = 3;
		const int col = 3;
		Mat_32F a(row, col);
		mat_rand(a, 0, 100);
		Mat_32F x(row, col);
		mat_rand(x, 0, 100);
		Mat_32F b(row, col);
		mat_rand(b, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		CalcTime t;
		for (int i = 0; i < loop; i++)
		{
			t.start();	//時間計測開始
			ret = mat_add(mat_mul(a, x), b);
			t.end();	//時間計測終了
			std::cout << "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;

		a.show();
		x.show();
		b.show();
		ret.show();

		return 0;
	}

	//課題2
	//行列積を計算するプログラムにおいて，コンパイラオプションを変えて計算速度の計測し，その違いを観測せよ．
	if (false)
	{
		std::cout << "課題2" << std::endl;
		const int loop = 10;
		const int row = 3;
		const int col = 3;
		Mat_32F a(row, col);
		mat_rand(a, 0, 100);
		Mat_32F b(row, col);
		mat_rand(b, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		CalcTime t;
		for (int i = 0; i < loop; i++)
		{
			t.start();
			ret = mat_mul(a, b);
			t.end();
			std::cout << "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;

		a.show();
		b.show();
		ret.show();

		return 0;
	}

	//課題3
	//小さな行列に対して，各要素を
	//3x^4+3x^3+3
	//する計算するプログラムを作成し，乗算回数を削減する前と後で計算速度を比較せよ．
	if (false)
	{
		std::cout << "課題3" << std::endl;
		const int loop = 10;
		const int row = 3;
		const int col = 3;
		Mat_32F x(row, col);
		mat_rand(x, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		//before
		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			//3x^4+3x^3+3 そのまま
			t.start();
			for (int j = 0; j < x.rows; j++)
			{
				for (int i = 0; i < x.cols; i++)
				{
					//計算
					//XXXX
					//XXXX
				}
			}
			t.end();
			std::cout << "before: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "before: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//after
		for (int k = 0; k < loop; k++)
		{
			//3x^4+3x^3+3 ホーナー法つかった場合
			t.start();
			for (int j = 0; j < x.rows; j++)
			{
				for (int i = 0; i < x.cols; i++)
				{
					//計算
					//XXXX
				}
			}
			t.end();
			std::cout << "after: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "after: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題4
	//小さな行列に対して，各要素を下記の定数倍するプログラムを作成し，数式の展開前後で計算速度を比較せよ．
	//(2π+sqrt(5)+0.5^2)x
	if (false)
	{
		std::cout << "課題4" << std::endl;
		const int loop = 10;
		const int row = 3;
		const int col = 3;
		Mat_32F x(row, col);
		mat_rand(x, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		//before
		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			//毎回計算する場合
			t.start();
			for (int j = 0; j < x.rows; j++)
			{
				for (int i = 0; i < x.cols; i++)
				{
					//計算
					//XXXX
				}
			}
			t.end();
			std::cout << "before: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "before: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//after
		const float c = (2 * 3.14 + sqrt(5) + 0.5*0.5);
		for (int k = 0; k < loop; k++)
		{
			//先に計算する場合
			t.start();
			for (int j = 0; j < x.rows; j++)
			{
				for (int i = 0; i < x.cols; i++)
				{
					//計算
					//XXXX
				}
			}
			t.end();
			std::cout << "after: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "after: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題5
	//小さな行列に対して，各要素を3.141592で除算する計算するプログラムを作成し，除算を削減する前と後で計算速度を比較せよ．
	//大きな行列で行うと，効果が少ない可能性があるため注意すること．
	if (false)
	{
		std::cout << "課題5" << std::endl;
		const int loop = 10;
		const int row = 3;
		const int col = 3;
		Mat_32F x(row, col);
		mat_rand(x, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

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
					//計算
					//XXXX
				}
			}
			t.end();
			std::cout << "div: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "div: time (avg): " << t.getAvgTime() << " ms" << std::endl;

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
					//XXXX
				}
			}
			t.end();
			std::cout << "mul: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "mul: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題6
	//小さな４つの行列A,B,C,Dに対して，行列の各要素ごとに`(a/b)*(c/d)`を計算するプログラムを作成し，除算を削減する前と後で計算速度を比較せよ．
	if (false)
	{
		std::cout << "課題6" << std::endl;
		const int loop = 10;
		const int row = 3;
		const int col = 3;
		Mat_32F a(row, col);
		mat_rand(a, 0, 100);
		Mat_32F b(row, col);
		mat_rand(b, 0, 100);
		Mat_32F c(row, col);
		mat_rand(c, 0, 100);
		Mat_32F d(row, col);
		mat_rand(d, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);


		CalcTime t;
		//before
		for (int k = 0; k < loop; k++)
		{
			//普通に計算
			t.start();
			for (int j = 0; j < ret.rows; j++)
			{
				for (int i = 0; i < ret.cols; i++)
				{
					//計算
					//XXXX
				}
			}
			t.end();
			std::cout << "before: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "before: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//after
		for (int k = 0; k < loop; k++)
		{
			//除算を削減した場合
			t.start();
			for (int j = 0; j < ret.rows; j++)
			{
				for (int i = 0; i < ret.cols; i++)
				{
					//計算
					//XXXX
				}
			}
			t.end();
			std::cout << "after: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "after: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題7
	//累乗を３乗，4乗．．．ｎ乗としたときに速度がどうなるのか計測せよ．ただし，累乗をすると値が大きくなるため，浮動小数点の最大値を超える可能性がある．その時は入力データを0.1倍`（i*0.1）`などすること．
	if (false)
	{
		std::cout << "課題7" << std::endl;
		const int loop = 1000;
		const int n = 50;
		const float v = 2.5f;
		float ret = 0;

		for (int i = 0; i < n; i++)
		{
			CalcTime t;
			//v^i乗を計算
			for (int j = 0; j < loop; j++)
			{
				t.start();
				//pow，計算
				//XXXX
				t.end();
				//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
			}
			std::cout << i << " : time (avg): " << t.getAvgTime() << " ms" << std::endl;
			//std::cout << ret << std::endl;
		}
		return 0;
	}


	//課題8
	//２つの行列の和を`unsigned char, short, int, float, double`で計算しそれぞれ比較せよ．
	//なお，大きい行列サイズでないと，効果がでない場合がある．
	if (false)
	{
		std::cout << "課題8" << std::endl;
		const int loop = 1000;
		const int row = 64;
		const int col = 64;

		//unsigend char
		Mat_8U a_8u(row, col);
		mat_rand(a_8u, 0, 100);
		Mat_8U b_8u(row, col);
		mat_rand(b_8u, 0, 100);
		Mat_8U ret_8u(row, col);
		mat_zero(ret_8u);

		CalcTime t;
		for (int i = 0; i < loop; i++)
		{
			t.start();
			//unsigned char
			//XXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << " 8U : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//short
		Mat_16S a_16s(row, col);
		mat_rand(a_16s, 0, 100);
		Mat_16S b_16s(row, col);
		mat_rand(b_16s, 0, 100);
		Mat_16S ret_16s(row, col);
		mat_zero(ret_16s);

		for (int i = 0; i < loop; i++)
		{
			t.start();
			//short
			//XXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "16S : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//int
		Mat_32S a_32s(row, col);
		mat_rand(a_32s, 0, 100);
		Mat_32S b_32s(row, col);
		mat_rand(b_32s, 0, 100);
		Mat_32S ret_32s(row, col);
		mat_zero(ret_32s);

		for (int i = 0; i < loop; i++)
		{
			t.start();
			//int
			//XXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "32S : time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//float
		Mat_32F a_32f(row, col);
		mat_rand(a_32f, 0, 100);
		Mat_32F b_32f(row, col);
		mat_rand(b_32f, 0, 100);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);

		for (int i = 0; i < loop; i++)
		{
			t.start();
			//float
			//XXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "32F : time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//double
		Mat_64F a_64f(row, col);
		mat_rand(a_64f, 0, 100);
		Mat_64F b_64f(row, col);
		mat_rand(b_64f, 0, 100);
		Mat_64F ret_64f(row, col);
		mat_zero(ret_64f);

		for (int i = 0; i < loop; i++)
		{
			t.start();
			//float
			//XXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "64F : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}


	//課題9
	//intの行列を整数で2倍，浮動小数点で2.f倍,整数を１ビットだけビットシフトすることで2倍する場合の計算速度を比較せよ．
	//また，intの行列を整数で2で除算する場合，浮動小数点で2で除算する場合，浮動小数点の0.5で乗算する場合，１ビットだけビットシフトすることで1/2倍する場合の速度を比較せよ．
	//加えて，floatの行列で，2.0で除算する場合と0.5で乗算する場合を比較せよ．
	//なお，浮動小数点で乗算する場合は整数の場合よりも遅い． 
	//また，大きい行列サイズでないと，効果がでない場合がある．
	if (false)
	{
		std::cout << "課題9" << std::endl;
		const int loop = 1000;
		const int row = 50;
		const int col = 50;
		
		//int
		Mat_32S x_32s(row, col);
		mat_rand(x_32s, 0, 100);
		Mat_32S ret_32s(row, col);
		mat_zero(ret_32s);

		CalcTime t;
		//2x mul
		for (int k = 0; k < loop; k++)
		{
			//2倍 乗算
			t.start();
			for (int j = 0; j < ret_32s.rows; j++)
			{
				for (int i = 0; i < ret_32s.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "2x mul: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		
		//2.0x mul
		for (int k = 0; k < loop; k++)
		{
			//2.0倍 乗算(double)
			t.start();
			for (int j = 0; j < ret_32s.rows; j++)
			{
				for (int i = 0; i < ret_32s.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "2x mul (double): time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//2x bit shift
		for (int k = 0; k < loop; k++)
		{
			//2倍 ビットシフト
			t.start();
			for (int j = 0; j < ret_32s.rows; j++)
			{
				for (int i = 0; i < ret_32s.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "2x bit shift: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//1/2 div int
		for (int k = 0; k < loop; k++)
		{
			//1/2 除算
			t.start();
			for (int j = 0; j < ret_32s.rows; j++)
			{
				for (int i = 0; i < ret_32s.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "1/2 div: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		
		//1/2 div float
		for (int k = 0; k < loop; k++)
		{
			//1/2 除算
			t.start();
			for (int j = 0; j < ret_32s.rows; j++)
			{
				for (int i = 0; i < ret_32s.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "1/2 div (float): time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//1/2 -> mul 0.5
		for (int k = 0; k < loop; k++)
		{
			//1/2 0.5乗算で実現
			t.start();
			for (int j = 0; j < ret_32s.rows; j++)
			{
				for (int i = 0; i < ret_32s.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "1/2 mul: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//1/2->bit shift
		for (int k = 0; k < loop; k++)
		{
			//1/2 ビットシフト
			t.start();
			for (int j = 0; j < ret_32s.rows; j++)
			{
				for (int i = 0; i < ret_32s.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "1/2 bit shift: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//float
		Mat_32F x_32f(row, col);
		mat_rand(x_32f, 0, 100);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);
		
		//1/2 div
		for (int k = 0; k < loop; k++)
		{
			//1/2 除算
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "float: 1/2 div: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//1/2 -> mul 0.5
		for (int k = 0; k < loop; k++)
		{
			//1/2 0.5乗算
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "float: 1/2 mul: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題10
	//floatの行列をを3.141倍する場合と，intの行列を3.141倍を固定小数点で行う場合で計算し比較せよ．
	if (false)
	{
		std::cout << "課題10" << std::endl;
		const int loop = 1000;
		const int row = 3;
		const int col = 3;
		Mat_32S x_32s(row, col);
		mat_rand(x_32s, 0, 100);
		Mat_32S ret_32s(row, col);
		mat_zero(ret_32s);

		//int
		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			//固定小数点
			t.start();
			for (int j = 0; j < ret_32s.rows; j++)
			{
				for (int i = 0; i < ret_32s.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "fixed: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		Mat_32F x_32f(row, col);
		mat_rand(x_32f, 0, 100);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);

		//float
		for (int k = 0; k < loop; k++)
		{
			//浮動小数点
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "float: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題11
	//floatの行列への定数値の四則演算と，`sin, cos, exp, log, sqrt`関数の適用した場合と計算時間を比較せよ．
	//また，`sin, cos, exp, log, sqrt`計算はテーブル参照も作成した場合についても比較せよ．
	//なお，環境によっては，演算したほうが速い演算もある可能性がある．
	if (false)
	{
		std::cout << "課題11" << std::endl;
		const int loop = 100;
		const int row = 50;
		const int col = 50;
		Mat_8U x_8u(row, col);
		mat_rand(x_8u, 0, 255);

		Mat_32F x_32f(x_8u);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);

		//四則演算
		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			//四則演算
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "arithmetic: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//sin
		for (int k = 0; k < loop; k++)
		{
			//sin関数
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "sin: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//cos
		for (int k = 0; k < loop; k++)
		{
			//cos関数
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "cos: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//exp
		for (int k = 0; k < loop; k++)
		{
			//exp関数
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "exp: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//log
		for (int k = 0; k < loop; k++)
		{
			//log関数
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "log: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//sqrt
		for (int k = 0; k < loop; k++)
		{
			//sqrt関数
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "sqrt: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		std::cout << std::endl;
		//sin LUT
		float LUT[256];
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXX
		}
		for (int k = 0; k < loop; k++)
		{
			t.start();
			//sin LUT
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					ret_32f.data[j*ret_32f.cols + i] = LUT[(int)x_32f.data[j*ret_32f.cols + i]];
				}
			}
			t.end();
		}
		std::cout << "sin LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//cos LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//cos LUT
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "cos LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//exp LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//exp LUT
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "exp LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//log LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//log LUT
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "log LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//sqrt LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//sqrt LUT
			t.start();
			for (int j = 0; j < ret_32f.rows; j++)
			{
				for (int i = 0; i < ret_32f.cols; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "sqrt LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}


	//課題12
	//小さな行列A,Bの各要素を任意のradianだけ回転させて，x,yにして格納するプログラムを記述し，inline展開の有無で速度がどのように変わるか計測せよ．
	//また，関数をべた書きした場合とも比較せよ．
	//ただし，-O3のオプションを付けると強制的にinline展開される可能性がある．
	//inline void rot(double a, double b, double &x, double &y, double radian)
	//{
	//	x = a * cos(radian);
	//	y = b * sin(radian);
	//}
	if (false)
	{
		std::cout << "課題12" << std::endl;
		const int loop = 10000;
		const int row = 3;
		const int col = 3;
		Mat_64F a(row, col);
		mat_rand(a, 0, 100);
		Mat_64F b(row, col);
		mat_rand(b, 0, 100);

		Mat_64F x(row, col);
		mat_zero(x);
		Mat_64F y(row, col);
		mat_zero(y);

		const float radian = 2.2f;
		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < a.rows; j++)
			{
				for (int i = 0; i < a.cols; i++)
				{
					//rot
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;
		return 0;
	}

	//課題13
	//行列A，Bの各要素の乗算を行うときに，結果を行列Cに格納する場合と行列Aに上書きする場合との計算時間をせよ．
	if (false)
	{
		std::cout << "課題13" << std::endl;
		const int loop = 10000;
		const int row = 3;
		const int col = 3;
		//c = a x b
		{
			Mat_64F c(row, col);
			mat_zero(c);

			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				//C=A*B
				Mat_64F a(row, col);
				mat_rand(a, 0, 100);
				Mat_64F b(row, col);
				mat_rand(b, 0, 100);
				for (int j = 0; j < a.rows; j++)
				{
					for (int i = 0; i < a.cols; i++)
					{
						//XXXX
					}
				}
				t.end();
			}
			std::cout << "c=a*b: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		}

		//a = a x b
		{

			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				//A=AxB
				Mat_64F a(row, col);
				mat_rand(a, 0, 100);
				Mat_64F b(row, col);
				mat_rand(b, 0, 100);
				for (int j = 0; j < a.rows; j++)
				{
					for (int i = 0; i < a.cols; i++)
					{
						//XXXX
					}
				}
				t.end();
			}
			std::cout << "a=a*b: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		}

		return 0;
	}

	//課題14
	//上記の0で初期化するコードをループの順序を変えてどちらが速いか計測して検証せよ．
	//また，行列積のコードのループの順序を変えてどれが速いか計測して検証せよ．
	if (false)
	{
		std::cout << "課題14" << std::endl;
		const int loop = 1000;
		const int width = 64;
		const int height = 64;
		float x[width][height];

		CalcTime t;
		//col, row
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int i = 0; i < width; ++i)
			{
				for (int j = 0; j < height; ++j)
				{
					x[i][j] = 0.0f;
				}
			}
			t.end();
		}
		std::cout << "col-row: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//row, col
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int i = 0; i < width; ++i)
			{
				for (int j = 0; j < height; ++j)
				{
					x[i][j] = 0.0f;
				}
			}
			t.end();
		}
		std::cout << "row-col: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		const int size = 64;
		float a[size][size];
		float b[size][size];
		float c[size][size];
		//init
		for (int j = 0; j < size; j++)
		{
			for (int i = 0; i < size; i++)
			{
				a[j][i] = rand_32f(0, 100);
				b[j][i] = rand_32f(0, 100);
				c[j][i] = 0;
			}
		}

		//i, j, k
		for (int l = 0; l < loop; l++)
		{
			t.start();
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < size; ++j)
				{
					for (int k = 0; k < size; ++k)
					{
						c[i][j] = c[i][j] + a[i][k] * b[k][j];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "i-j-k: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//i, k, j
		for (int l = 0; l < loop; l++)
		{
			t.start();
			for (int i = 0; i < size; ++i)
			{
				for (int k = 0; k < size; ++k)
				{
					for (int j = 0; j < size; ++j)
					{
						//XXXX
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "i-k-j: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//j, i, k
		for (int l = 0; l < loop; l++)
		{
			t.start();
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; ++i)
				{
					for (int k = 0; k < size; ++k)
					{
						//XXXX
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "j-i-k: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//j, k, i
		for (int l = 0; l < loop; l++)
		{
			t.start();
			for (int j = 0; j < size; ++j)
			{
				for (int k = 0; k < size; ++k)
				{
					for (int i = 0; i < size; ++i)
					{
						//XXXX
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "j-k-i: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//k, i, j
		for (int l = 0; l < loop; l++)
		{
			t.start();
			for (int k = 0; k < size; ++k)
			{
				for (int i = 0; i < size; ++i)
				{
					for (int j = 0; j < size; ++j)
					{
						//XXXX
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "k-i-j: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//k, j, i
		for (int l = 0; l < loop; l++)
		{
			t.start();
			for (int k = 0; k < size; ++k)
			{
				for (int j = 0; j < size; ++j)
				{
					for (int i = 0; i < size; ++i)
					{
						//XXXX
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "k-j-i: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}


	//課題15
	//アンローリングの段数を2,4,8,16,32,...と変更することで，速度がどのように変わるか計測せよ．
	if (false)
	{
		std::cout << "課題15" << std::endl;
		const int loop = 100;
		const int size = 1024;
		float x[size], y[size];
		float a = 2.f;
		float b = 1.f;

		for (int i = 0; i < size; i++)
		{
			x[i] = rand_32f(0, 100);
		}

		CalcTime t;
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
		std::cout << "no unrolling: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 2
		for (int j = 0; j < loop; j++)
		{
			//unrolling 2
			t.start();
			for (int i = 0; i < size; i += 2)
			{
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "unrolling  2: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 4
		for (int j = 0; j < loop; j++)
		{
			//unrolling 4
			t.start();
			for (int i = 0; i < size; i += 4)
			{
				//XXXX
				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "unrolling  4: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 8
		for (int j = 0; j < loop; j++)
		{
			//unrolling 8
			t.start();
			for (int i = 0; i < size; i += 8)
			{
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
		std::cout << "unrolling  8: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 16
		for (int j = 0; j < loop; j++)
		{
			//unrolling 16
			t.start();
			for (int i = 0; i < size; i += 16)
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
			}
			t.end();
		}
		std::cout << "unrolling 16: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 32
		for (int j = 0; j < loop; j++)
		{
			//unrolling 32
			t.start();
			for (int i = 0; i < size; i += 32)
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
			}
			t.end();
		}
		std::cout << "unrolling 32: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 64
		for (int j = 0; j < loop; j++)
		{
			//unrolling 64
			t.start();
			for (int i = 0; i < size; i += 64)
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
		std::cout << "unrolling 64: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題16
	//上記のプログラムを実装し，ループピーリングの有無で速度がどのように変わるか計測せよ．
	if (false)
	{
		std::cout << "課題16" << std::endl;
		const int loop = 100;
		const int size = 1024;
		int x[size], y[size];

		for (int i = 0; i < size; i++)
		{
			x[i] = rand_32f(0, 100);
			y[i] = rand_32f(0, 100);
		}

		CalcTime t;
		//original
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; ++i)
			{
				if (i == 0)
				{
					y[i] = (x[i] + x[i + 1]) / 2;
				}
				else if (i == size - 1)
				{
					y[i] = (x[i - 1] + x[i]) / 2;
				}
				else
				{
					y[i] = (x[i - 1] + x[i] + x[i + 1]) / 3;
				}
			}
			t.end();
		}
		std::cout << "original: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//Loop peeling
		for (int j = 0; j < loop; j++)
		{
			t.start();
			{
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "loop peeling: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題17
	//上記のコードで，計算時間を比較せよ．
	//なお，必ずしも差がでるとは限らない（配列の大きさを変更すると差が出やすい）．
	if (false)
	{
		std::cout << "課題17" << std::endl;
		const int loop = 100;
		const int width = 768;
		const int height = 512;
		float a[height][width];

		CalcTime t;
		// before
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < height; j++)
			{
				for (int i = 0; i < width; i++)
				{
					a[j][i] = 0.f;
				}
			}
			t.end();
		}
		std::cout << "before: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//after
		const int size = width * height;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			float *p = &a[0][0];
			for (int i = 0; i < size; i++) {
				*p++ = 0.f;
			}
			t.end();
		}
		std::cout << "after: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題18
	//上記のコードを実行し，並列に動作していることを確認せよ．
	//また，並列化を有効にする場合としない場合の計算時間を比較せよ．
	if (false)
	{
		std::cout << "課題18" << std::endl;
		//XXXX
		for (int i = 0; i < 100; i++)
		{
			std::cout << i << std::endl; //並列化したい処理
		}
		return 0;
	}

	//課題19
	//総和を計算するコードで，reduction指定子を使用する場合としない場合で計算結果がどのようになるか比較せよ．
	if (false)
	{
		std::cout << "課題19" << std::endl;
		//並列化なし
		std::cout << "no parallelization" << std::endl;
		for (int j = 0; j < 10; j++)
		{
			int sum = 0;
			for (int i = 0; i < 100; i++)
			{
				sum += i;
			}
			std::cout << sum << std::endl;
		}

		//omp使用
		std::cout << "parallelization" << std::endl;
		for (int j = 0; j < 10; j++)
		{
			int sum = 0;
			//XXXX
			for (int i = 0; i < 100; i++)
			{
				sum += i;
			}
			std::cout << sum << std::endl;
		}

		//omp reduction使用
		std::cout << "parallelization with reduction" << std::endl;
		for (int j = 0; j < 10; j++)
		{
			int sum = 0;
			//XXXX
			for (int i = 0; i < 100; i++)
			{
				sum += i;
			}
			std::cout << sum << std::endl;
		}
		return 0;
	}

	//課題20
	//二つの行列の各要素の積を計算するコードで，スレッド数を変更して，計算時間がどのように推移するのかを確認せよ．
	if (false)
	{
		std::cout << "課題20" << std::endl;
		const int loop = 1000;
		const int size = 128;
		Mat_32F a(size, size);
		mat_rand(a, 0, 100);
		Mat_32F b(size, size);
		mat_rand(b, 0, 100);
		Mat_32F c(size, size);

		CalcTime t;
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			//omp num_threads(n)で並列化，nを入れる
			//XXXX
			for (int i = 0; i < size; ++i)
			{
				for (int k = 0; k < size; ++k)
				{
					for (int j = 0; j < size; ++j)
					{
						c.data[i*c.cols + j] = c.data[i*c.cols + j] + a.data[i*a.cols + k] * b.data[k*b.cols + j];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;
		return 0;
	}

	//課題21
	//四則演算のコードを書いてprintfデバッグで確認せよ．
	if (false)
	{
		std::cout << "課題21" << std::endl;
		const __m256 a = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
		const __m256 b = _mm256_set_ps(15, 14, 13, 12, 11, 10, 9, 8);

		__m256 c;
		//加算
		//XXXX
		std::cout << "add: ";
		print_m256(c);

		//減算
		//XXXX
		std::cout << "sub: ";
		print_m256(c);

		//乗算
		//XXXX
		std::cout << "mul: ";
		print_m256(c);

		//除算
		//XXXX
		std::cout << "div: ";
		print_m256(c);

		return 0;
	}


	//課題22
	//配列a,x,bに対して，`(((a*x+b)*x+b)*x+b)*x+b `の計算を配列ｃに格納するコードをmul/addで記述するものとFMAを使うもので記述し，FMAが速くなることを示せ．
	//なお，上記の関数は以下に等しい．
	//a=_mm256_fmadd_ps(a,b,c);
	//a=_mm256_fmadd_ps(a,b,c);
	//a=_mm256_fmadd_ps(a,b,c);
	//a=_mm256_fmadd_ps(a,b,c);
	//これは，単純にFMAが1度だとメモリで律速するこのコードでは計算速度の差が出にくいためである．差が小さければ，より演算を増やせば良い．
	if (false)
	{
		std::cout << "課題22" << std::endl;
		const int loop = 1000;
		const int size = 1024;
#ifdef __GNUC__
		float __attribute__((aligned(32))) a[size];
		float __attribute__((aligned(32))) x[size];
		float __attribute__((aligned(32))) b[size];
		float __attribute__((aligned(32))) c[size];
#elif _MSC_VER
		float __declspec(align(32)) a[size];
		float __declspec(align(32)) x[size];
		float __declspec(align(32)) b[size];
		float __declspec(align(32)) c[size];
#endif

		//init
		for (int i = 0; i < size; i++)
		{
			a[i] = rand_32f(0, 100);
			x[i] = rand_32f(0, 100);
			b[i] = rand_32f(0, 100);
			c[i] = 0;
		}

		CalcTime t;
		//mul, add
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mx = _mm256_load_ps(x + i);
				const __m256 mb = _mm256_load_ps(b + i);

				//mul,addを使って
				__m256 temp;
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				_mm256_store_ps(c + i, temp);
			}
			t.end();
		}
		std::cout << "mul-add: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//fma
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mx = _mm256_load_ps(x + i);
				const __m256 mb = _mm256_load_ps(b + i);

				//fmaを使って
				__m256 temp;
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				_mm256_store_ps(c + i, temp);
			}
			t.end();
		}
		std::cout << "fma: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}


	//課題23
	//divとrcp,sqrtとrsqrtの実行速度を比較せよ．
	if (false)
	{
		std::cout << "課題23" << std::endl;
		const int loop = 10000;
		const int size = 1024;
#ifdef __GNUC__
		float __attribute__((aligned(32))) a[size];
		float __attribute__((aligned(32))) b[size];
#elif _MSC_VER
		__declspec(align(32)) float a[size];
		__declspec(align(32)) float b[size];
#endif

		for (int i = 0; i < size; i++)
		{
			a[i] = rand_32f(0, 100);
			b[i] = rand_32f(0, 100);
		}

		CalcTime t;
		//div
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mb = _mm256_load_ps(b + i);

				__m256 temp;
				//divを使って
//XXXX
			}
			t.end();
		}
		std::cout << "div: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//rcp
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mb = _mm256_load_ps(b + i);

				__m256 temp;
				//rcpをつかって
//XXXX
			}
			t.end();
		}
		std::cout << "rcp: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//sqrt
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mb = _mm256_load_ps(b + i);

				__m256 temp;
				//sqrtを使って
//XXXX
			}
			t.end();
		}
		std::cout << " sqrt: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//rsqrt
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mb = _mm256_load_ps(b + i);

				__m256 temp;
				//rsqrtを使って
//XXXX
			}
			t.end();
		}
		std::cout << "rsqrt: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}


	//課題24
	//haddとdpで要素の総和を取るプログラムを作成せよ．
	//また，それぞれの計算時間を比較せよ ．
	//（この課題は，後にリダクションの最適化でもう一度登場します．）
	if (false)
	{
		std::cout << "課題24" << std::endl;
		const int loop = 100000;
#ifdef __GNUC__
		float __attribute__((aligned(32))) a[8];
		float __attribute__((aligned(32))) b[8];
#elif _MSC_VER
		__declspec(align(32)) float a[8];
		__declspec(align(32)) float b[8];
#endif

		for (int i = 0; i < 8; i++)
		{
			a[i] = rand_32f(0, 100);
			b[i] = rand_32f(0, 100);
		}

		CalcTime t;
		//hadd
		for (int i = 0; i < loop; i++)
		{
			t.start();
			__m256 ma = _mm256_load_ps(a);
			//haddを使って
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			t.end();
		}
		std::cout << "hadd: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		const __m256 one = _mm256_set1_ps(1);
		//dp
		for (int i = 0; i < loop; i++)
		{
			t.start();
			__m256 mb = _mm256_load_ps(b);
			//dpを使って
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			t.end();
		}
		std::cout << "dp: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題25
	//行列aにおいて要素の値があるしきい値以上の場合だけ3乗し，それ以外は何もしない処理をベクトル化実装せよ．
	//if(a[i]>=threshold) a[i]=a[i]*a[i]*a[i];
	if (false)
	{
		std::cout << "課題25" << std::endl;
		const int size = 8;
		const float th = 50;

#ifdef __GNUC__
		float __attribute__((aligned(32))) a[size];
#elif _MSC_VER
		__declspec(align(32)) float a[size];
#endif

		for (int i = 0; i < size; i++)
		{
			a[i] = rand_32f(0, 100);
		}

		const __m256 mth = _mm256_set1_ps(th);
		for (int i = 0; i < size; i += 8)
		{
			const __m256 ma = _mm256_load_ps(a + i);

			__m256 temp;
			//cmp, mul, blendvを使って
			//XXXX
			//XXXX
			//XXXX
			_mm256_store_ps(a + i, temp);
		}
		for (int i = 0; i < size; i++)
		{
			std::cout << a[i] << ", ";
			if (i != 0 && i%size == 0)
				std::cout << std::endl;
		}
		std::cout << std::endl;
		return 0;
	}

	//課題26
	//上記のコードのように，SIMD命令を使う場合におけるループアンローリングを8，16，32，64と行い，計算時間を比較せよ．
	if (false)
	{
		std::cout << "課題26" << std::endl;
		const int loop = 1000;
		const int size = 1024 * 3;
#ifdef __GNUC__
		__attribute__((aligned(32))) float a[size];
		__attribute__((aligned(32))) float b[size];
		__attribute__((aligned(32))) float c[size];
#elif _MSC_VER
		__declspec(align(32)) float a[size];
		__declspec(align(32)) float b[size];
		__declspec(align(32)) float c[size];
#endif

		//init
		for (int i = 0; i < size; i++)
		{
			a[i] = rand_32f(0, 100);
			b[i] = rand_32f(0, 100);
			c[i] = 0;
		}

		CalcTime t;
		// unrolling 8
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 8
			for (int i = 0; i < size; i += 8)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
				_mm256_store_ps(c + i, _mm256_add_ps(ma, mb));
			}
			t.end();
		}
		std::cout << "  8: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		// unrolling 16
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 16
			for (int i = 0; i < size; i += 16)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
				//XXXX

				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << " 16: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		// unrolling 32
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 32
			for (int i = 0; i < size; i += 32)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
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
		std::cout << " 32: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		// unrolling 64
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 64
			for (int i = 0; i < size; i += 64)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
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
		std::cout << " 64: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		// unrolling 128
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 128
			for (int i = 0; i < size; i += 128)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
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
		std::cout << "128: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題27
	//上の関数を用いて，データ構造の相互変換を確認せよ． ４ｘ４のdouble型データの転置を手書きで考えてみよ．
	//この問題は，自力で答えを埋めると時間がかかり過ぎるためプログラムの答えはすでに埋めてある．
	//そのため，この動作を追って転置ができていることを手書きで確認するだけで良い．
	if (false)
	{
		std::cout << "課題27" << std::endl;
		const int size = 64;
#ifdef __GNUC__
		__attribute__((aligned(32))) float a[size];
		__attribute__((aligned(32))) float b[size];
#elif _MSC_VER
		__declspec(align(32)) float a[size];
		__declspec(align(32)) float b[size];
#endif
		for (int i = 0; i < size; i++)
		{
			a[i] = i;
			b[i] = 0;
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
//XXXX

		for (int i = 0; i < 8; i++)
		{
			print_m256(mb[i]);
		}
		return 0;
	}

	//課題28
	//__m256i（int）型を_m256（float）型に変換せよ．
	if (false)
	{
		std::cout << "課題28" << std::endl;
		__m256i m32i = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
		__m256 m32f = _mm256_setzero_ps();

		print_m256(m32f);
		//cvtを使って
//XXXX
		print_m256(m32f);
	}


	//課題29
	//上記のコードを実行し，計算時間を計測せよ．
	//また，スカラ実装，スカラ実装＋並列化，SIMD実装のみを作成し，計算時間を比較せよ．
	if (false)
	{
		std::cout << "課題29" << std::endl;
		const int loop = 1000;
		const int size = 128;
		Mat_32F a(size, size);
		Mat_32F b(size, size);
		Mat_32F c(size, size);

		//init
		mat_rand(a, 0, 100);
		mat_rand(b, 0, 100);
		mat_zero(c);

		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			//スカラー実装
			t.start();
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "scalar    : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//スカラー，並列化実装
			t.start();
			//XXXX
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "scalar+omp: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		for (int k = 0; k < loop; k++)
		{
			//SIMD実装
			t.start();
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i += 8)
				{
					//XXXX
					//XXXX
					//XXXX
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "SIMD      : time (avg): " << t.getAvgTime() << " ms" << std::endl;


		for (int k = 0; k < loop; k++)
		{
			//SIMD，並列化実装
			t.start();
			//XXXX
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i += 8)
				{
					//XXXX
					//XXXX
					//XXXX
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "SIMD+omp  : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	std::cout << "no select" << std::endl;
	return 0;
}
