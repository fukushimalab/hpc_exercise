#include "utils/image_util.h"
#include "utils/simd_util.h"

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <algorithm>
#include <omp.h>

/////////////////////////////////
// box filter
/////////////////////////////////
//サンプルコード
//演習1
void boxFilter_scalar_AoS_32F(const Image_32F& src, Image_32F& dest, const int r);
//演習3
void boxFilter_simd_SoA_32F_pixelParallel_kernelUnrolling(const Image_32F& src, Image_32F& dest, const int r);
//演習4
void boxFilter_scalar_SoA_32F_pixelParallel_separable(const Image_32F& src, Image_32F& dest, const int r);

//各自実装する関数
//演習1
void boxFilter_scalar_SoA_32F(const Image_32F& src, Image_32F& dest, const int r);
//演習2
void boxFilter_scalar_SoA_32F_pixelParallelY(const Image_32F& src, Image_32F& dest, const int r);
void boxFilter_scalar_SoA_32F_pixelParallelX(const Image_32F& src, Image_32F& dest, const int r);
//演習3
void boxFilter_simd_SoA_32F_pixelParallel_pixelUnrolling(const Image_32F& src, Image_32F& dest, const int r);
//演習4
void boxFilter_simd_SoA_32F_pixelParallel_pixelUnrolling_separable(const Image_32F& src, Image_32F& dest, const int r);

//演習（オプション）
void boxFilter_scalar_SoA_8U_pixelUnrollingParallel(const Image_32F& src, Image_32F& dest, const int r);

int main(const int argc, const char** argv)
{
	if (argc == 1)
	{
		std::cout << "./boxfilter num_work num_iterations radius" << std::endl;
	}

	const int default_loop = 10;
	const int default_r = 3;
	const int default_work = 0;

	const int work = (argc < 2) ? default_work : atoi(argv[1]);
	const int loop = (argc < 3) ? default_loop : atoi(argv[2]);
	const int r = (argc < 4) ? default_r : atoi(argv[3]);

	std::cout << "work " << work << ": iteration = " << loop << " r = " << r << std::endl << std::endl;

	//演習1
	//1. カラー画像におけるボックスフィルタをスカラーで実装せよ．このとき，入力画像のデータ型をfloatとし，画像のデータ構造をSoAに変換する場合の実装を行え．これらの実装については，カラーループは展開せよ．なお，AoSのスカラ実装はサンプルコードで提供している．
	//2. カーネル半径rを1から10程度まで変更して，サンプルで提供されているAoS実装と1.の課題で実装したSoA実装の計算時間を計測せよ．また，splitとmergeにかかる計算時間も計測し，フィルタ処理本体とデータ構造変換を含めた場合の計算時間も求めよ．
	//3. AoS実装とSoA実装の計算時間を比較し，計算時間の違いをキャッシュ効率の観点などから考察せよ．また，データ構造変換による計算時間の短縮とデータ構造変換にかかる計算時間より，データ構造変換のトレードオフについても議論せよ．
	if (work == 1)
	{
		Image_8U src_8u;
		readPXM("img/lena.ppm", src_8u);
		Image_32F src_32f(src_8u);
		Image_32F dest_32f_AoS;
		Image_32F dest_32f_SoA;

		std::cout << "|time(avg) method|time [ms]  |" << std::endl;
		std::cout << "|----------------|-----------|" << std::endl;

		//AoS
		//実装済み
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_scalar_AoS_32F(src_32f, dest_32f_AoS, r);
				t.end();
			}
			std::cout << "|AoS filter      |" << t.getAvgTime() << " ms|" << std::endl;
		}
		//SoA
		//各自実装せよ
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_scalar_SoA_32F(src_32f, dest_32f_SoA, r);
				t.end();
			}
			std::cout << "|AoS->SoA->f->AoS|" << t.getAvgTime() << " ms|" << std::endl;
		}

		double psnr = calcPSNR(dest_32f_AoS, dest_32f_SoA);
		std::cout << std::endl << "PSNR : " << psnr << " dB" << std::endl;
		if (psnr < 50.0)std::cout << "invalid image" << std::endl;
		Image_8U dest_8u_AoS(dest_32f_AoS);
		writePXM("imgout/AoS.ppm", dest_8u_AoS);
		Image_8U dest_8u_SoA(dest_32f_SoA);
		writePXM("imgout/SoA.ppm", dest_8u_SoA);
		return 0;
	}

	//演習2
	//1. 演習1で実装したスカラーコードをOpenMPでスレッド並列化せよ．並列化するループは画素ループとし，yループとxループをそれぞれ並列化せよ．
	//2. 1.で実装した画素yループ並列化と画素xループ並列化の計算時間を測定し，その結果を比較して最も外側のループを並列化する場合が最も効率が良いことを確認せよ．
	//3. 画素ループ並列化のコードにおいて，並列化数を変更して計算時間を確認せよ．
	if (work == 2)
	{
		Image_8U src_8u;
		readPXM("img/lena.ppm", src_8u);
		Image_32F src_32f(src_8u);
		Image_32F dest_32f_pixel_parallel_y(src_8u);
		Image_32F dest_32f_pixel_parallel_x(src_8u);

		std::cout << "|time(avg) method|time [ms]  |" << std::endl;
		std::cout << "|----------------|-----------|" << std::endl;

		//single
		//課題１で実装済み
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_scalar_SoA_32F(src_32f, dest_32f_pixel_parallel_y, r);
				t.end();
			}
			std::cout << "|single          |" << t.getAvgTime() << " ms|" << std::endl;
		}

		//pixel parallel y
		//各自実装せよ
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_scalar_SoA_32F_pixelParallelY(src_32f, dest_32f_pixel_parallel_y, r);
				t.end();
			}
			std::cout << "|y loop paralel  |" << t.getAvgTime() << " ms|" << std::endl;
		}
		//pixel parallel x
		//各自実装せよ
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_scalar_SoA_32F_pixelParallelX(src_32f, dest_32f_pixel_parallel_x, r);
				t.end();
			}
			std::cout << "|x loop paralel  |" << t.getAvgTime() << " ms|" << std::endl;
		}

		double psnr = calcPSNR(dest_32f_pixel_parallel_y, dest_32f_pixel_parallel_x);
		std::cout << std::endl << "PSNR : " << psnr << " dB" << std::endl;
		if (psnr < 50.0)std::cout << "invalid image" << std::endl;
		Image_8U dest_8u_pixel_parallel_y(dest_32f_pixel_parallel_y);
		writePXM("imgout/pixel_parallel_y.ppm", dest_8u_pixel_parallel_y);
		Image_8U dest_8u_pixel_parallel_x(dest_32f_pixel_parallel_x);
		writePXM("imgout/pixel_parallel_x.ppm", dest_8u_pixel_parallel_x);
		return 0;
	}

	//演習3
	//1. SIMD命令を使用しボックスフィルタを実装せよ．
	//2. 1.で実装したコードに対応するスカラ実装，ベクトル化実装の計算時間をカーネルの半径rをr=1から1刻みづつr=16程度まで変更して測定せよ．また，サンプルコードとして提供されているカーネルループ展開についても，同様に計測せよ．
	//3. 2.の結果より，SIMD実装にしたことによる高速化率や，画素ループ展開とカーネルループ展開の効率の違いや，余り処理の影響などを考察せよ．
	if (work == 3)
	{
		Image_8U src_8u;
		readPXM("img/lena.ppm", src_8u);
		Image_32F src_32f(src_8u);
		Image_32F dest_32f_scalar(src_8u);
		Image_32F dest_32f_kernel_unrolling(src_8u);
		Image_32F dest_32f_pixel_unrolling(src_8u);

		std::cout << "|time(avg) method|time [ms]  |" << std::endl;
		std::cout << "|----------------|-----------|" << std::endl;

		//scalar
		//演習2で実装
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_scalar_SoA_32F_pixelParallelY(src_32f, dest_32f_scalar, r);
				t.end();
			}
			std::cout << "|pixel scalar    | " << t.getAvgTime() << " ms|" << std::endl;
		}
		//kernel unrolling
		//実装済み
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_simd_SoA_32F_pixelParallel_kernelUnrolling(src_32f, dest_32f_kernel_unrolling, r);
				t.end();
			}
			std::cout << "|kernel SIMD     | " << t.getAvgTime() << " ms|" << std::endl;
		}
		//pixel unrolling
		//各自実装せよ
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_simd_SoA_32F_pixelParallel_pixelUnrolling(src_32f, dest_32f_pixel_unrolling, r);
				t.end();
			}
			std::cout << "|pixel SIMD      | " << t.getAvgTime() << " ms|" << std::endl;
		}

		std::cout << std::endl;
		double psnr_k = calcPSNR(dest_32f_scalar, dest_32f_kernel_unrolling);
		double psnr_p = calcPSNR(dest_32f_scalar, dest_32f_pixel_unrolling);
		std::cout << "|PSNR: method    |PSNR  |" << std::endl;
		std::cout << "|----------------|------|" << std::endl;
		std::cout << "|PSNR kernel SIMD|" << psnr_k << " dB|" << std::endl;
		if (psnr_k < 50.0)std::cout << "invalid image: kernel unrolling" << std::endl;
		std::cout << "|PSNR pixel  SIMD|" << psnr_p << " dB|" << std::endl;
		if (psnr_p < 50.0)std::cout << "invalid image: pixel_unrolling" << std::endl;

		Image_8U dest_8u_kernel_unrolling(dest_32f_kernel_unrolling);
		writePXM("imgout/kernel_unrolling.ppm", dest_8u_kernel_unrolling);
		Image_8U dest_8u_pixel_unrolling(dest_32f_pixel_unrolling);
		writePXM("imgout/pixel_unrolling.ppm", dest_8u_pixel_unrolling);
		return 0;
	}

	//演習4
	//1. SIMD実装を縦，横のセパラブルフィルタに拡張せよ．なお，サンプルとして，スカラーの実装が示されている．
	//2. rを1から10程度まで変更しながら，スカラ実装，SIMD実装，スカラ・セパラブル実装，1.で実装したSIMD・セパラブル実装の計算時間を測定せよ．
	//3. 2.で計測した計算時間より，セパラブルフィルタによる高速化率が`O(r)`と`O(r^2)`との関係に従うか調べ，考察せよ．
	if (work == 4)
	{
		Image_8U src_8u;
		readPXM("img/lena.ppm", src_8u);
		Image_32F src_32f(src_8u);
		Image_32F dest_32f_scalar(src_8u);
		Image_32F dest_32f_simd(src_8u);
		Image_32F dest_32f_scalar_separable(src_8u);
		Image_32F dest_32f_simd_separable(src_8u);

		std::cout << "|time(avg) method|time [ms]|" << std::endl;
		std::cout << "|----------------|---------|" << std::endl;
		//scalar
		//演習1で実装
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_scalar_SoA_32F(src_32f, dest_32f_scalar, r);
				t.end();
			}
			std::cout << "|scalar single   |" << t.getAvgTime() << " ms|" << std::endl;
		}
		//演習2で実装
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_scalar_SoA_32F_pixelParallelY(src_32f, dest_32f_scalar, r);
				t.end();
			}
			std::cout << "|scalar parallel |" << t.getAvgTime() << " ms|" << std::endl;
		}
		//simd unrolling
		//演習3で実装
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_simd_SoA_32F_pixelParallel_pixelUnrolling(src_32f, dest_32f_simd, r);
				t.end();
			}
			std::cout << "|SIMD parallel   |" << t.getAvgTime() << " ms|" << std::endl;
		}
		//scalar separable
		//実装済み
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_scalar_SoA_32F_pixelParallel_separable(src_32f, dest_32f_scalar_separable, r);
				t.end();
			}
			std::cout << "|scalar separable|" << t.getAvgTime() << " ms|" << std::endl;
		}
		//SIMD separable
		//各自実装せよ
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_simd_SoA_32F_pixelParallel_pixelUnrolling_separable(src_32f, dest_32f_simd_separable, r);
				t.end();
			}
			std::cout << "|SIMD   separable|" << t.getAvgTime() << " ms|" << std::endl;
		}


		double psnr_simd = calcPSNR(dest_32f_scalar, dest_32f_simd);
		double psnr_sepscalar = calcPSNR(dest_32f_scalar, dest_32f_scalar_separable);
		double psnr_sepsimd = calcPSNR(dest_32f_scalar, dest_32f_simd_separable);
		std::cout << std::endl;
		std::cout << "|PSNR: method          |PSNR   |" << std::endl;
		std::cout << "|----------------------|-------|" << std::endl;
		std::cout << "|PSNR: SIMD            | " << psnr_simd << " dB|" << std::endl;
		if (psnr_simd < 50)std::cout << "invalid image: SIMD" << std::endl;
		std::cout << "|PSNR: scalar separable| " << psnr_sepscalar << " dB|" << std::endl;
		if (psnr_sepscalar < 50)std::cout << "invalid image: scalar separable" << std::endl;
		std::cout << "|PSNR: SIMD   separable| " << psnr_sepsimd << " dB|" << std::endl;
		if (psnr_sepsimd < 50)std::cout << "invalid image: SIMD separable" << std::endl;

		Image_8U dest_8u_scalar(dest_32f_scalar);
		writePXM("imgout/scalar.ppm", dest_8u_scalar);
		Image_8U dest_8u_simd(dest_32f_simd);
		writePXM("imgout/simd.ppm", dest_8u_simd);
		Image_8U dest_8u_scalar_separable(dest_32f_scalar_separable);
		writePXM("imgout/scalar_separable.ppm", dest_8u_scalar_separable);
		Image_8U dest_8u_simd_separable(dest_32f_simd_separable);
		writePXM("imgout/simd_separable.ppm", dest_8u_simd_separable);

		return 0;
	}

	std::cout << "no selected: box filter" << std::endl;
}


///////////////////////
// ボックスフィルタ
///////////////////////
void boxFilter_scalar_AoS_32F(const Image_32F& src, Image_32F& dest, const int r)
{
	Image_32F temp_32f;
	//copyMakeboder
	copyMakeBorder(src, temp_32f, r, r, r, r);

	//出力画像作成
	dest = Image_32F(src.rows, src.cols, src.channels);
	//重み合計
	const float wsum = (float)((2 * r + 1) * (2 * r + 1));

	if (src.channels == 3)
	{
		for (int y = 0; y < src.rows; y++)
		{
			for (int x = 0; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				for (int j = -r; j <= r; j++)
				{
					for (int i = -r; i <= r; i++)
					{
						sum[0] += temp_32f.data[((y + r + j) * temp_32f.cols + (x + r + i)) * 3 + 0];
						sum[1] += temp_32f.data[((y + r + j) * temp_32f.cols + (x + r + i)) * 3 + 1];
						sum[2] += temp_32f.data[((y + r + j) * temp_32f.cols + (x + r + i)) * 3 + 2];
					}
				}
				dest.data[3 * (y * dest.cols + x) + 0] = sum[0] / wsum;
				dest.data[3 * (y * dest.cols + x) + 1] = sum[1] / wsum;
				dest.data[3 * (y * dest.cols + x) + 2] = sum[2] / wsum;
			}
		}
	}
	else if (src.channels == 1)
	{
		//gray image
		std::cout << "not support gray image" << std::endl;
		return;
	}
}

////////////////////////////
//演習1
////////////////////////////
void boxFilter_scalar_SoA_32F(const Image_32F& src, Image_32F& dest, const int r)
{
	Image_32F temp_32f;
	//copyMakeboder
	copyMakeBorder(src, temp_32f, r, r, r, r);

	//重み合計
	const float wsum = (float)((2 * r + 1) * (2 * r + 1));

	if (src.channels == 3)
	{
		Image_32F splitImg[3];
		split(temp_32f, splitImg);//AoS->SoA

		//出力画像の初期化
		Image_32F destSplitImg[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

		for (int y = 0; y < src.rows; y++)
		{
			for (int x = 0; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				for (int j = -r; j <= r; j++)
				{
					for (int i = -r; i <= r; i++)
					{
						//sum処理
						//ここを実装
					}
				}
				//store処理
				//ここを実装
			}
		}
		merge(destSplitImg, 3, dest);//SoA->AoS
	}
	else if (src.channels == 1)
	{
		//gray image
		std::cout << "not support gray image" << std::endl;
		return;
	}
}


////////////////////////////
//演習2
////////////////////////////
void boxFilter_scalar_SoA_32F_pixelParallelY(const Image_32F& src, Image_32F& dest, const int r)
{
	//演習1のコードを使って実装
}

void boxFilter_scalar_SoA_32F_pixelParallelX(const Image_32F& src, Image_32F& dest, const int r)
{
	//演習1のコードを使って実装
}

////////////////////////////
//演習3
////////////////////////////
void boxFilter_simd_SoA_32F_pixelParallel_kernelUnrolling(const Image_32F& src, Image_32F& dest, const int r)
{
	Image_32F temp_32f;
	//copyMakeboder
	copyMakeBorder(src, temp_32f, r, r, r, r);

	//重み合計
	const float wsum = (float)((2 * r + 1) * (2 * r + 1));
	const int residualNum = (2 * r + 1) % 8;

	if (src.channels == 3)
	{
		Image_32F splitImg[3];
		split(temp_32f, splitImg);

		//出力画像の初期化
		Image_32F destSplitImg[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			for (int x = 0; x < src.cols; x++)
			{
				__m256 msum0 = _mm256_setzero_ps();
				__m256 msum1 = _mm256_setzero_ps();
				__m256 msum2 = _mm256_setzero_ps();
				for (int j = -r; j <= r; j++)
				{
					int i = -r;
					for (; i <= r - residualNum; i += 8)
					{
						__m256 temp0 = _mm256_loadu_ps(&splitImg[0].data[(y + r + j) * splitImg[0].cols + (x + r + i)]);
						__m256 temp1 = _mm256_loadu_ps(&splitImg[1].data[(y + r + j) * splitImg[1].cols + (x + r + i)]);
						__m256 temp2 = _mm256_loadu_ps(&splitImg[2].data[(y + r + j) * splitImg[2].cols + (x + r + i)]);

						msum0 = _mm256_add_ps(temp0, msum0);
						msum1 = _mm256_add_ps(temp1, msum1);
						msum2 = _mm256_add_ps(temp2, msum2);
					}
					//余り処理
					for (; i <= r; i++)
					{
						((float*)&msum0)[0] += splitImg[0].data[(y + r + j) * splitImg[0].cols + (x + r + i)];
						((float*)&msum1)[0] += splitImg[1].data[(y + r + j) * splitImg[1].cols + (x + r + i)];
						((float*)&msum2)[0] += splitImg[2].data[(y + r + j) * splitImg[2].cols + (x + r + i)];
					}
				}
				msum0 = _mm256_hadd_ps(msum0, msum0);
				msum0 = _mm256_hadd_ps(msum0, msum0);
				float sum0 = ((float*)&msum0)[0] + ((float*)&msum0)[4];
				msum1 = _mm256_hadd_ps(msum1, msum1);
				msum1 = _mm256_hadd_ps(msum1, msum1);
				float sum1 = ((float*)&msum1)[0] + ((float*)&msum1)[4];
				msum2 = _mm256_hadd_ps(msum2, msum2);
				msum2 = _mm256_hadd_ps(msum2, msum2);
				float sum2 = ((float*)&msum2)[0] + ((float*)&msum2)[4];
				destSplitImg[0].data[(y * destSplitImg[0].cols + x)] = sum0 / wsum;
				destSplitImg[1].data[(y * destSplitImg[1].cols + x)] = sum1 / wsum;
				destSplitImg[2].data[(y * destSplitImg[2].cols + x)] = sum2 / wsum;
			}
		}
		merge(destSplitImg, 3, dest);
	}
	else if (src.channels == 1)
	{
		//gray image
		std::cout << "not support gray image" << std::endl;
		return;
	}

}

void boxFilter_simd_SoA_32F_pixelParallel_pixelUnrolling(const Image_32F& src, Image_32F& dest, const int r)
{
	Image_32F temp_32f;
	//copyMakeboder
	copyMakeBorder(src, temp_32f, r, r, r, r);

	//重み合計
	const float wsum = (float)((2 * r + 1) * (2 * r + 1));
	const __m256 mwsum = _mm256_set1_ps(wsum);
	const int residualNum = temp_32f.cols % 8;

	if (src.channels == 3)
	{
		Image_32F splitImg[3];
		split(temp_32f, splitImg);

		//出力画像の初期化
		Image_32F destSplitImg[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			int x = 0;
			for (; x < src.cols - residualNum; x += 8)
			{
				__m256 sum0 = _mm256_setzero_ps();
				__m256 sum1 = _mm256_setzero_ps();
				__m256 sum2 = _mm256_setzero_ps();
				for (int j = -r; j <= r; j++)
				{
					for (int i = -r; i <= r; i++)
					{
						//sum処理
						//ここを実装

					}
				}
				//store処理
				//ここを実装

			}
			//余り処理
			for (; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				for (int j = -r; j <= r; j++)
				{
					for (int i = -r; i <= r; i++)
					{
						sum[0] += splitImg[0].data[(y + r + j) * splitImg[0].cols + (x + r + i)];
						sum[1] += splitImg[1].data[(y + r + j) * splitImg[1].cols + (x + r + i)];
						sum[2] += splitImg[2].data[(y + r + j) * splitImg[2].cols + (x + r + i)];
					}
				}
				destSplitImg[0].data[(y * destSplitImg[0].cols + x)] = sum[0] / wsum;
				destSplitImg[1].data[(y * destSplitImg[1].cols + x)] = sum[1] / wsum;
				destSplitImg[2].data[(y * destSplitImg[2].cols + x)] = sum[2] / wsum;
			}
		}
		merge(destSplitImg, 3, dest);
	}
	else if (src.channels == 1)
	{
		//gray image
		std::cout << "not support gray image" << std::endl;
		return;
	}

}

////////////////////////////
//演習4
////////////////////////////
void boxFilter_scalar_SoA_32F_pixelParallel_separable(const Image_32F& src, Image_32F& dest, const int r)
{
	//重み合計
	const float wsum = (float)((2 * r + 1) * (2 * r + 1));

	if (src.channels == 3)
	{
		Image_32F temp_32f;
		//copyMakeboder
		copyMakeBorder(src, temp_32f, r, r, 0, 0);

		Image_32F splitImg[3];
		Image_32F vfilterImg[3];
		split(temp_32f, splitImg);

		//出力画像の初期化
		Image_32F destSplitImg0[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

		//縦フィルタ
#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			for (int x = 0; x < src.cols; x++)
			{
				float sum[3] = { 0.f };
				for (int j = -r; j <= r; j++)
				{
					sum[0] += splitImg[0].data[(y + r + j) * splitImg[0].cols + x];
					sum[1] += splitImg[1].data[(y + r + j) * splitImg[1].cols + x];
					sum[2] += splitImg[2].data[(y + r + j) * splitImg[2].cols + x];
				}
				destSplitImg0[0].data[(y * destSplitImg0[0].cols + x)] = sum[0];
				destSplitImg0[1].data[(y * destSplitImg0[1].cols + x)] = sum[1];
				destSplitImg0[2].data[(y * destSplitImg0[2].cols + x)] = sum[2];
			}
		}
		copyMakeBorder(destSplitImg0[0], vfilterImg[0], 0, 0, r, r);
		copyMakeBorder(destSplitImg0[1], vfilterImg[1], 0, 0, r, r);
		copyMakeBorder(destSplitImg0[2], vfilterImg[2], 0, 0, r, r);
		Image_32F destSplitImg1[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

		//横フィルタ
#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			for (int x = 0; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				for (int i = -r; i <= r; i++)
				{
					sum[0] += vfilterImg[0].data[y * vfilterImg[0].cols + (x + r + i)];
					sum[1] += vfilterImg[1].data[y * vfilterImg[1].cols + (x + r + i)];
					sum[2] += vfilterImg[2].data[y * vfilterImg[2].cols + (x + r + i)];
				}
				destSplitImg1[0].data[(y * destSplitImg1[0].cols + x)] = sum[0] / wsum;
				destSplitImg1[1].data[(y * destSplitImg1[1].cols + x)] = sum[1] / wsum;
				destSplitImg1[2].data[(y * destSplitImg1[2].cols + x)] = sum[2] / wsum;
			}
		}
		merge(destSplitImg1, 3, dest);
	}
	else if (src.channels == 1)
	{
		//gray image
		std::cout << "not support gray image" << std::endl;
		return;
	}
}

void boxFilter_simd_SoA_32F_pixelParallel_pixelUnrolling_separable(const Image_32F& src, Image_32F& dest, const int r)
{
	//重み合計
	const float wsum = (float)((2 * r + 1) * (2 * r + 1));

	if (src.channels == 3)
	{
		Image_32F temp_32f;
		//copyMakeboder
		copyMakeBorder(src, temp_32f, r, r, 0, 0);
		const int residualNum = temp_32f.cols % 8;

		Image_32F splitImg[3];
		Image_32F vfilterImg[3];
		split(temp_32f, splitImg);

		//出力画像の初期化
		Image_32F destSplitImg0[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

		//縦フィルタ
		//これを実装
#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			int x = 0;
			for (; x < src.cols - residualNum; x += 8)
			{
				__m256 sum0 = _mm256_setzero_ps();
				__m256 sum1 = _mm256_setzero_ps();
				__m256 sum2 = _mm256_setzero_ps();
				for (int j = -r; j <= r; j++)
				{
					//sum処理
					//ここを実装

				}
				//store処理
				//ここを実装

			}
			//余り処理
			for (; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				for (int j = -r; j <= r; j++)
				{
					sum[0] += splitImg[0].data[(y + r + j) * splitImg[0].cols + x];
					sum[1] += splitImg[1].data[(y + r + j) * splitImg[1].cols + x];
					sum[2] += splitImg[2].data[(y + r + j) * splitImg[2].cols + x];
				}
				destSplitImg0[0].data[(y * destSplitImg0[0].cols + x)] = sum[0];
				destSplitImg0[1].data[(y * destSplitImg0[1].cols + x)] = sum[1];
				destSplitImg0[2].data[(y * destSplitImg0[2].cols + x)] = sum[2];
			}
		}
		copyMakeBorder(destSplitImg0[0], vfilterImg[0], 0, 0, r, r);
		copyMakeBorder(destSplitImg0[1], vfilterImg[1], 0, 0, r, r);
		copyMakeBorder(destSplitImg0[2], vfilterImg[2], 0, 0, r, r);
		Image_32F destSplitImg1[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

		//横フィルタ
		//これを実装
#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			int x = 0;
			for (; x < src.cols - residualNum; x += 8)
			{
				__m256 sum0 = _mm256_setzero_ps();
				__m256 sum1 = _mm256_setzero_ps();
				__m256 sum2 = _mm256_setzero_ps();
				for (int i = -r; i <= r; i++)
				{
					//sum処理
					//ここを実装

				}
				//store処理
				//ここを実装

			}
			for (; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				for (int i = -r; i <= r; i++)
				{
					sum[0] += vfilterImg[0].data[y * vfilterImg[0].cols + (x + r + i)];
					sum[1] += vfilterImg[1].data[y * vfilterImg[1].cols + (x + r + i)];
					sum[2] += vfilterImg[2].data[y * vfilterImg[2].cols + (x + r + i)];
				}
				destSplitImg1[0].data[(y * destSplitImg1[0].cols + x)] = sum[0] / wsum;
				destSplitImg1[1].data[(y * destSplitImg1[1].cols + x)] = sum[1] / wsum;
				destSplitImg1[2].data[(y * destSplitImg1[2].cols + x)] = sum[2] / wsum;
			}
		}
		merge(destSplitImg1, 3, dest);
	}
	else if (src.channels == 1)
	{
		//gray image
		std::cout << "not support gray image" << std::endl;
		return;
	}
}

////////////////////////////
//演習（オプション）
////////////////////////////
void boxFilter_scalar_SoA_8U_pixelUnrollingParallel(const Image_32F& src, Image_32F& dest, const int r)
{
}
