#include "utils/image.h"
#include "utils/image_util.h"
#include "utils/simd_util.h"
#include <time.h>
#include <math.h>

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
	//演習1
	//if (false)
	{
		const int loop = 10;

		const int r = 3;
		Image_8U src_8u;
		readPXM("img/lena.ppm", src_8u);
		Image_32F src_32f(src_8u);
		Image_32F dest_32f_AoS;
		Image_32F dest_32f_SoA;

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
			std::cout << "AoS: time (avg): " << t.getAvgTime() << " ms" << std::endl;
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
			std::cout << "SoA: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		}

		std::cout << "PSNR : " << calcPSNR(dest_32f_AoS, dest_32f_SoA) << " dB" << std::endl; 
		Image_8U dest_8u_AoS(dest_32f_AoS);
		writePXM("AoS.ppm", dest_8u_AoS);
		Image_8U dest_8u_SoA(dest_32f_SoA);
		writePXM("SoA.ppm", dest_8u_SoA);
		return 0;
	}

	//演習2
	if (false)
	{
		const int loop = 10;

		const int r = 3;
		Image_8U src_8u;
		readPXM("img/lena.ppm", src_8u);
		Image_32F src_32f(src_8u);
		Image_32F dest_32f_pixel_parallel_y;
		Image_32F dest_32f_pixel_parallel_x;

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
			std::cout << "y: time (avg): " << t.getAvgTime() << " ms" << std::endl;
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
			std::cout << "x: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		}


		std::cout << "PSNR : " << calcPSNR(dest_32f_pixel_parallel_y, dest_32f_pixel_parallel_x) << " dB" << std::endl; 
		Image_8U dest_8u_pixel_parallel_y(dest_32f_pixel_parallel_y);
		writePXM("pixel_parallel_y.ppm", dest_8u_pixel_parallel_y);
		Image_8U dest_8u_pixel_parallel_x(dest_32f_pixel_parallel_x);
		writePXM("pixel_parallel_x.ppm", dest_8u_pixel_parallel_x);
		return 0;
	}

	//演習3
	if (false)
	{
		const int loop = 10;

		const int r = 3;
		Image_8U src_8u;
		readPXM("img/lena.ppm", src_8u);
		Image_32F src_32f(src_8u);
		Image_32F dest_32f_scalar;
		Image_32F dest_32f_kernel_unrolling;
		Image_32F dest_32f_pixel_unrolling;
		
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
			std::cout << "pixel: time (avg): " << t.getAvgTime() << " ms" << std::endl;
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
			std::cout << "kernel: time (avg): " << t.getAvgTime() << " ms" << std::endl;
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
			std::cout << "pixel: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		}


		std::cout << "PSNR kernel unrolling: " << calcPSNR(dest_32f_scalar, dest_32f_kernel_unrolling) << " dB" << std::endl; 
		std::cout << "PSNR pixel unrolling: " << calcPSNR(dest_32f_scalar, dest_32f_pixel_unrolling) << " dB" << std::endl; 
		Image_8U dest_8u_kernel_unrolling(dest_32f_kernel_unrolling);
		writePXM("kernel_unrolling.ppm", dest_8u_kernel_unrolling);
		Image_8U dest_8u_pixel_unrolling(dest_32f_pixel_unrolling);
		writePXM("pixel_unrolling.ppm", dest_8u_pixel_unrolling);
		return 0;
	}

	//演習4
	if (false)
	{
		const int loop = 10;

		const int r = 3;
		Image_8U src_8u;
		readPXM("img/lena.ppm", src_8u);
		Image_32F src_32f(src_8u);
		Image_32F dest_32f_scalar;
		Image_32F dest_32f_simd;
		Image_32F dest_32f_scalar_separable;
		Image_32F dest_32f_simd_separable;

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
			std::cout << "scalar: time (avg): " << t.getAvgTime() << " ms" << std::endl;
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
			std::cout << "simd: time (avg): " << t.getAvgTime() << " ms" << std::endl;
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
			std::cout << "scalar separable: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		}
		//simd separable
		//各自実装せよ
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				boxFilter_simd_SoA_32F_pixelParallel_pixelUnrolling_separable(src_32f, dest_32f_simd_separable, r);
				t.end();
			}
			std::cout << "simd separable: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		}



		std::cout << "PSNR: simd: " << calcPSNR(dest_32f_scalar, dest_32f_simd) << " dB" << std::endl; 
		std::cout << "PSNR: scalar separable: " << calcPSNR(dest_32f_scalar, dest_32f_scalar_separable) << " dB" << std::endl; 
		std::cout << "PSNR: simd separable: " << calcPSNR(dest_32f_scalar, dest_32f_simd_separable) << " dB" << std::endl; 
		Image_8U dest_8u_scalar(dest_32f_scalar);
		writePXM("scalar.ppm", dest_8u_scalar);
		Image_8U dest_8u_simd(dest_32f_simd);
		writePXM("simd.ppm", dest_8u_simd);
		Image_8U dest_8u_scalar_separable(dest_32f_scalar_separable);
		writePXM("scalar_separable.ppm", dest_8u_scalar_separable);
		Image_8U dest_8u_simd_separable(dest_32f_simd_separable);
		writePXM("simd_separable.ppm", dest_8u_simd_separable);

		return 0;
	}
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
	const float wsum = (2*r+1)*(2*r+1);

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
						sum[0] += temp_32f.data[((y + r + j)*temp_32f.cols + (x + r + i))*3 + 0];
						sum[1] += temp_32f.data[((y + r + j)*temp_32f.cols + (x + r + i))*3 + 1];
						sum[2] += temp_32f.data[((y + r + j)*temp_32f.cols + (x + r + i))*3 + 2];
					}
				}
				dest.data[3 * (y*dest.cols + x) + 0] = sum[0] / wsum + 0.5f;
				dest.data[3 * (y*dest.cols + x) + 1] = sum[1] / wsum + 0.5f;
				dest.data[3 * (y*dest.cols + x) + 2] = sum[2] / wsum + 0.5f;
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
	const float wsum = (2*r+1)*(2*r+1);

	if (src.channels == 3)
	{
		Image_32F splitImg[3];
		split(temp_32f, splitImg);

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
	const float wsum = (2*r+1)*(2*r+1);
	const int residualNum = (2*r+1)%8;

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
					for (; i <= r-residualNum; i+=8)
					{
						__m256 temp0 = _mm256_loadu_ps(&splitImg[0].data[(y + r + j)*splitImg[0].cols + (x + r + i)]);
						__m256 temp1 = _mm256_loadu_ps(&splitImg[1].data[(y + r + j)*splitImg[1].cols + (x + r + i)]);
						__m256 temp2 = _mm256_loadu_ps(&splitImg[2].data[(y + r + j)*splitImg[2].cols + (x + r + i)]);

						msum0 = _mm256_add_ps(temp0, msum0);
						msum1 = _mm256_add_ps(temp1, msum1);
						msum2 = _mm256_add_ps(temp2, msum2);
					}
					//余り処理
					for(; i <= r; i++)
					{
						((float*)&msum0)[0] += splitImg[0].data[(y + r + j)*splitImg[0].cols + (x + r + i)];
						((float*)&msum1)[0] += splitImg[1].data[(y + r + j)*splitImg[1].cols + (x + r + i)];
						((float*)&msum2)[0] += splitImg[2].data[(y + r + j)*splitImg[2].cols + (x + r + i)];
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
				destSplitImg[0].data[(y*destSplitImg[0].cols + x)] = sum0 / wsum + 0.5f;
				destSplitImg[1].data[(y*destSplitImg[1].cols + x)] = sum1 / wsum + 0.5f;
				destSplitImg[2].data[(y*destSplitImg[2].cols + x)] = sum2 / wsum + 0.5f;
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
	const float wsum = (2*r+1)*(2*r+1);
	const __m256 mwsum = _mm256_set1_ps(wsum);
	const __m256 offset = _mm256_set1_ps(0.5f);
	const int residualNum = temp_32f.cols%8;

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
			for (; x < src.cols-residualNum; x+=8)
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
						sum[0] += splitImg[0].data[(y + r + j)*splitImg[0].cols + (x + r + i)];
						sum[1] += splitImg[1].data[(y + r + j)*splitImg[1].cols + (x + r + i)];
						sum[2] += splitImg[2].data[(y + r + j)*splitImg[2].cols + (x + r + i)];
					}
				}
				destSplitImg[0].data[(y*destSplitImg[0].cols + x)] = sum[0] / wsum + 0.5f;
				destSplitImg[1].data[(y*destSplitImg[1].cols + x)] = sum[1] / wsum + 0.5f;
				destSplitImg[2].data[(y*destSplitImg[2].cols + x)] = sum[2] / wsum + 0.5f;
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
	const float wsum = (2*r+1)*(2*r+1);

	if (src.channels == 3)
	{
		Image_32F temp_32f;
		//copyMakeboder
		copyMakeBorder(src, temp_32f, r, r, 0, 0);

		Image_32F splitImg[3];
		split(temp_32f, splitImg);

		//出力画像の初期化
		Image_32F destSplitImg0[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			for (int x = 0; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				for (int j = -r; j <= r; j++)
				{
					sum[0] += splitImg[0].data[(y + r + j)*splitImg[0].cols + x];
					sum[1] += splitImg[1].data[(y + r + j)*splitImg[1].cols + x];
					sum[2] += splitImg[2].data[(y + r + j)*splitImg[2].cols + x];
				}
				destSplitImg0[0].data[(y*destSplitImg0[0].cols + x)] = sum[0];
				destSplitImg0[1].data[(y*destSplitImg0[1].cols + x)] = sum[1];
				destSplitImg0[2].data[(y*destSplitImg0[2].cols + x)] = sum[2];
			}
		}
		copyMakeBorder(destSplitImg0[0], destSplitImg0[0], 0, 0, r, r);
		copyMakeBorder(destSplitImg0[1], destSplitImg0[1], 0, 0, r, r);
		copyMakeBorder(destSplitImg0[2], destSplitImg0[2], 0, 0, r, r);
		Image_32F destSplitImg1[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			for (int x = 0; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				for (int i = -r; i <= r; i++)
				{
					sum[0] += destSplitImg0[0].data[y*destSplitImg0[0].cols + (x + r + i)];
					sum[1] += destSplitImg0[1].data[y*destSplitImg0[1].cols + (x + r + i)];
					sum[2] += destSplitImg0[2].data[y*destSplitImg0[2].cols + (x + r + i)];
				}
				destSplitImg1[0].data[(y*destSplitImg1[0].cols + x)] = sum[0] / wsum + 0.5f;
				destSplitImg1[1].data[(y*destSplitImg1[1].cols + x)] = sum[1] / wsum + 0.5f;
				destSplitImg1[2].data[(y*destSplitImg1[2].cols + x)] = sum[2] / wsum + 0.5f;
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
	const float wsum = (2*r+1)*(2*r+1);

	if (src.channels == 3)
	{
		Image_32F temp_32f;
		//copyMakeboder
		copyMakeBorder(src, temp_32f, r, r, 0, 0);
		const int residualNum = temp_32f.cols%8;

		Image_32F splitImg[3];
		split(temp_32f, splitImg);

		//出力画像の初期化
		Image_32F destSplitImg0[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

		//縦フィルタ
		//これを実装
#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			int x = 0;
			for (; x < src.cols-residualNum; x+=8)
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
					sum[0] += splitImg[0].data[(y + r + j)*splitImg[0].cols + x];
					sum[1] += splitImg[1].data[(y + r + j)*splitImg[1].cols + x];
					sum[2] += splitImg[2].data[(y + r + j)*splitImg[2].cols + x];
				}
				destSplitImg0[0].data[(y*destSplitImg0[0].cols + x)] = sum[0];
				destSplitImg0[1].data[(y*destSplitImg0[1].cols + x)] = sum[1];
				destSplitImg0[2].data[(y*destSplitImg0[2].cols + x)] = sum[2];
			}
		}
		copyMakeBorder(destSplitImg0[0], destSplitImg0[0], 0, 0, r, r);
		copyMakeBorder(destSplitImg0[1], destSplitImg0[1], 0, 0, r, r);
		copyMakeBorder(destSplitImg0[2], destSplitImg0[2], 0, 0, r, r);
		Image_32F destSplitImg1[3] = { Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1), Image_32F(src.rows, src.cols, 1) };

		//横フィルタ
		//実装済み
#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			for (int x = 0; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				for (int i = -r; i <= r; i++)
				{
					sum[0] += destSplitImg0[0].data[y*destSplitImg0[0].cols + (x + r + i)];
					sum[1] += destSplitImg0[1].data[y*destSplitImg0[1].cols + (x + r + i)];
					sum[2] += destSplitImg0[2].data[y*destSplitImg0[2].cols + (x + r + i)];
				}
				destSplitImg1[0].data[(y*destSplitImg1[0].cols + x)] = sum[0] / wsum + 0.5f;
				destSplitImg1[1].data[(y*destSplitImg1[1].cols + x)] = sum[1] / wsum + 0.5f;
				destSplitImg1[2].data[(y*destSplitImg1[2].cols + x)] = sum[2] / wsum + 0.5f;
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


