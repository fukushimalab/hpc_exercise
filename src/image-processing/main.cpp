#include "utils/image_util.h"
#include "utils/simd_util.h"

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <algorithm>

/////////////////////////////////
// チュートリアル
/////////////////////////////////
void GammaCorrection(const Image_8U& src, Image_8U& dest, const float gamma);
void GammaCorrectionFast(const Image_8U& src, Image_8U& dest, const float gamma);
void MeanVar(const Image_8U& src, float& mean, float& var);
void MeanVarFast(const Image_8U& src, float& mean, float& var);
void GaussianFilter(const Image_8U& src, Image_8U& dest, const int r, const float sigma);
void FIRFilter(const Image_8U& src, Image_8U& dest, const int r, float parameter);

/////////////////////////////////
// 総合課題（個別）
/////////////////////////////////
void BilateralFilter(const Image_8U& src, Image_8U& dest, const int r, const float sigma_r, const float sigma_s);
void NonLocalMeansFilter(const Image_8U& src, Image_8U& dest, const int template_r, const int search_r, const float h);


int main(const int argc, const char** argv)
{
	//////////////////////////////////////////////////////////////////////
	// チュートリアル
	//////////////////////////////////////////////////////////////////////

	if (false)
	{
		std::cout << "tutorial 1" << std::endl;
		Image_8U img;
		readPXM("img/lena.ppm", img);//カラー画像読み込み
		Image_8U gray;
		cvtColorGray(img, gray);//グレイに変換
		writePXM("img/lena_gray.pgm", gray);//グレイ画像保存

		return 0;
	}

	if (false)
	{
		std::cout << "tutorial 2" << std::endl;
		Image_8U img;
		readPXM("img/lena.ppm", img);//カラー画像読み込み
		Image_8U gray;
		cvtColorGray(img, gray);//グレイに変換

		for (int j = 0; j < img.rows; j++)
		{
			for (int i = 0; i < img.cols; i++)
			{
				if (i < 50 || j < 50 || i >= img.cols - 50 || j >= img.rows - 50)
				{
					gray.data[img.cols * j + i] = 0;
				}
			}
		}
		writePXM("img/lena_gray_pad.pgm", gray);//グレイ画像保存

		return 0;
	}

	if (false)
	{
		std::cout << "tutorial 3" << std::endl;
		Image_8U img;//入力画像
		readPXM("img/lena.ppm", img);//カラー画像読み込み
		Image_8U dst(img);//出力画像

		for (int j = 0; j < img.rows; j++)
		{
			for (int i = 0; i < img.cols; i++)
			{
				//カラーなのでデータが3倍．緑はRGBなので1番目を2倍
				dst.data[3 * (img.cols * j + i) + 1] = std::min(255, 2 * (img.data[3 * (img.cols * j + i) + 1]));
			}
		}
		writePXM("img/lena_green.ppm", dst);//カラー画像保存

		//splitとmergeを使ってもよい．こっちは青を2倍
		Image_8U splitImg[3];
		split(img, splitImg);
		for (int j = 0; j < img.rows; j++)
		{
			for (int i = 0; i < img.cols; i++)
			{
				//3倍がいらない．青を2倍するので2個目のアレイを使う．
				splitImg[2].data[img.cols * j + i] = std::min(255, 2 * (splitImg[2].data[img.cols * j + i]));
			}
		}
		merge(splitImg, 3, dst);
		writePXM("img/lena_blue.ppm", dst);//カラー画像保存

		return 0;
	}

	if (false)
	{
		std::cout << "tutorial 4" << std::endl;
		Image_8U img;//入力画像
		readPXM("img/lena.ppm", img);//カラー画像読み込み

		Image_8U dst(img);//出力画像
		copyMakeBorder(img, dst, 50, 50, 50, 50);

		//ループはimgのサイズで．
		for (int j = 0; j < img.rows; j++)
		{
			for (int i = 0; i < img.cols; i++)
			{
				//カラーなのでデータが3倍．緑はRGBなので1番目を2倍
				//ステップはdataのサイズで．また，位置は先頭から(50,50)ずれていることに注意．
				//カラー画像RGBなので3つ分アンロール
				dst.data[3 * (dst.cols * (j + 50) + i + 50) + 0] = 0;
				dst.data[3 * (dst.cols * (j + 50) + i + 50) + 1] = 0;
				dst.data[3 * (dst.cols * (j + 50) + i + 50) + 2] = 0;
			}
		}
		writePXM("img/lena_border.ppm", dst);//カラー画像保存

		return 0;
	}


	//////////////////////////////////////////////////////////////////////
	// 課題
	//////////////////////////////////////////////////////////////////////

	if (false)
	{
		std::cout << "Gaussian filter" << std::endl;

		const int loop = 10;

		const int r = 5;
		const float sigma = 20.f;
		Image_8U src, dest, reference;
		readPXM("img/lena.ppm", src);

		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			GaussianFilter(src, reference, r, sigma);
			t.end();
		}
		std::cout << "base: " << t.getAvgTime() << " ms" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			t.start();
			//最適化後．コードを一般化したのでガウシアンフィルタではなくてFIRフィルタと命名
			FIRFilter(src, dest, r, sigma);
			t.end();
		}
		std::cout << "opt.: " << t.getAvgTime() << " ms" << std::endl;
		//畳み込み等の浮動小数点演算の順序に依存するので，PSNRが無限大にすることは難しい．
		//ただし，PSNRは58以上は，unsigned charにキャストした場合の四捨五入の量子化の範囲内なので同一の画像と見なして良い．
		//また，人間の目で50dB以上の画像を区別するのはほぼ不可能．
		std::cout << "PSNR: " << calcPSNR(reference, dest) << "dB" << std::endl;

		writePXM("img/gauss.ppm", dest);
		return 0;
	}

	// gamma correction
	if (false)
	{
		std::cout << "gamma correction" << std::endl;
		const int loop = 10;

		const float gamma = 2.f;
		Image_8U src, dest, reference;
		readPXM("img/lena.ppm", src);

		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			GammaCorrection(src, reference, gamma);
			t.end();
		}
		std::cout << "base: " << t.getAvgTime() << " ms" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			t.start();
			//この関数を最適する
			GammaCorrectionFast(src, dest, gamma);
			t.end();
		}
		std::cout << "opt.: " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "PSNR: " << calcPSNR(reference, dest) << std::endl;

		writePXM("img/gamma.ppm", dest);
		return 0;
	}

	///////////////////////
	// mean-var
	///////////////////////
	if (false)
	{
		std::cout << "mean-var" << std::endl;
		const int loop = 10;

		Image_8U src;
		readPXM("img/lena.ppm", src);

		float mean, var;
		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			MeanVar(src, mean, var);
			t.end();
		}

		std::cout << "base: " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "mean: " << mean << std::endl;
		std::cout << "var : " << var << std::endl;

		for (int k = 0; k < loop; k++)
		{
			t.start();
			MeanVarFast(src, mean, var);
			t.end();
		}
		std::cout << "opt.: " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "mean: " << mean << std::endl;
		std::cout << "var : " << var << std::endl;
		return 0;
	}

	//////////////////////////////////////////////////////////////////////
	// 総合課題（個別）
	//////////////////////////////////////////////////////////////////////
	///////////////////////
	// bilateral filter
	///////////////////////
	if (false)
	{
		const int loop = 10;

		const float sigma_s = 1.f;
		const int r = (int)(3.f * sigma_s);
		const float sigma_r = 16.0f;
		Image_8U src, dest;
		readPXM("img/lena.ppm", src);

		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			BilateralFilter(src, dest, r, sigma_r, sigma_s);
			t.end();
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "PSNR : " << calcPSNR(src, dest) << " dB" << std::endl; //高速化した関数のとサンプルコードの出力の精度を確認すること（現状は，入力画像を入れている．）
		writePXM("bf.ppm", dest);
		return 0;
	}
	///////////////////////
	// non-local means filter
	///////////////////////
	if (false)
	{
		const int loop = 10;

		const float h = 20;
		const int search_r = 3;
		const int template_r = 1;
		Image_8U src, gray, dest;
		readPXM("img/lena.ppm", src);
		cvtColorGray(src, gray);

		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			NonLocalMeansFilter(gray, dest, template_r, search_r, h);
			t.end();
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "PSNR : " << calcPSNR(gray, dest) << " dB" << std::endl; //高速化した関数のとサンプルコードの出力の精度を確認すること（現状は，入力画像を入れている．）
		writePXM("nlmf.pgm", dest);
		return 0;
	}

	std::cout << "no selected image processing" << std::endl;
}

//////////////////////////////////////////////////////////////////////
// チュートリアル
//////////////////////////////////////////////////////////////////////
///////////////////////
// ガンマ変換
///////////////////////
void GammaCorrection(const Image_8U& src, Image_8U& dest, const float gamma)
{
	dest = Image_8U(src.rows, src.cols, src.channels);
	const int cn = src.channels;
	for (int y = 0; y < src.rows; y++)
	{
		for (int x = 0; x < src.cols * cn; x++)
		{
			dest.data[cn * (y * src.cols) + x] = (unsigned char)(pow((float)src.data[cn * (y * src.cols) + x] / 255.f, 1.f / gamma) * 255.0f);
		}
	}
}

void GammaCorrectionFast(const Image_8U& src, Image_8U& dest, const float gamma)
{
	dest = Image_8U(src.rows, src.cols, src.channels);
	const int cn = src.channels;
	for (int y = 0; y < src.rows; y++)
	{
		for (int x = 0; x < src.cols * cn; x++)
		{
			dest.data[cn * (y * src.cols) + x] = (unsigned char)(pow((float)src.data[cn * (y * src.cols) + x] / 255.f, 1.f / gamma) * 255.0f);
		}
	}
}

///////////////////////
// 平均・分散計算
///////////////////////
void MeanVar(const Image_8U& src, float& mean, float& var)
{
	mean = 0;
	var = 0;
	const int cn = src.channels;

	for (int y = 0; y < src.rows; y++)
	{
		for (int x = 0; x < src.cols * cn; x++)
		{
			mean += src.data[cn * (y * src.cols) + x];
			var += src.data[cn * (y * src.cols) + x] * src.data[cn * (y * src.cols) + x];
		}
	}
	mean /= (float)(src.rows * src.cols * cn);
	var = var / (float)(src.rows * src.cols * cn) - mean * mean;
}

void MeanVarFast(const Image_8U& src, float& mean, float& var)
{
	mean = 0;
	var = 0;
	const int cn = src.channels;

	for (int y = 0; y < src.rows; y++)
	{
		for (int x = 0; x < src.cols * cn; x++)
		{
			mean += src.data[cn * (y * src.cols) + x];
			var += src.data[cn * (y * src.cols) + x] * src.data[cn * (y * src.cols) + x];
		}
	}
	mean /= (float)(src.rows * src.cols * cn);
	var = var / (float)(src.rows * src.cols * cn) - mean * mean;
}

///////////////////////
// ガウシアンフィルタ
///////////////////////
void GaussianFilter(const Image_8U& src, Image_8U& dest, const int r, const float sigma)
{
	Image_8U temp_8u;
	copyMakeBorder(src, temp_8u, r, r, r, r);

	Image_32F temp_32f(temp_8u);
	Image_32F splitImg[3];
	split(temp_32f, splitImg);

	Image_32F dest_32f(src.rows, src.cols, src.channels);
	if (src.channels == 3)
	{
		for (int y = 0; y < src.rows; y++)
		{
			for (int x = 0; x < src.cols; x++)
			{
				float sum[3] = { 0 };
				float wsum = 0;
				for (int j = -r; j <= r; j++)
				{
					for (int i = -r; i <= r; i++)
					{
						const int distance = j * j + i * i;
						const float w = exp(-distance / (2.f * sigma * sigma));

						for (int c = 0; c < 3; c++)
						{
							sum[c] += splitImg[c].data[(y + r + j) * splitImg[c].cols + (x + r + i)] * w;
						}
						wsum += w;
					}
				}
				for (int c = 0; c < 3; c++)
				{
					dest_32f.data[3 * (y * dest_32f.cols + x) + c] = sum[c] / wsum;
				}
			}
		}
	}
	else if (src.channels == 1)
	{
		//gray image
	}
	dest = Image_8U(dest_32f);
}

///////////////////////
// FIR filter (sample)
///////////////////////
inline __m256 _mm256_load_epu8cvtps(const __m128i* P)
{
	return _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadu_si128((__m128i*)P)));
}

void inline _mm256_stream_ps_color(void* dst, const __m256 rsrc, const __m256 gsrc, const __m256 bsrc)
{
	static const int smask1 = _MM_SHUFFLE(1, 2, 3, 0);
	static const int smask2 = _MM_SHUFFLE(2, 3, 0, 1);
	static const int smask3 = _MM_SHUFFLE(3, 0, 1, 2);
	static const int bmask1 = 0x44;
	static const int bmask2 = 0x22;
	static const int pmask1 = 0x20;
	static const int pmask2 = 0x30;
	static const int pmask3 = 0x31;
	const __m256 aa = _mm256_shuffle_ps(rsrc, rsrc, smask1);
	const __m256 bb = _mm256_shuffle_ps(gsrc, gsrc, smask2);
	const __m256 cc = _mm256_shuffle_ps(bsrc, bsrc, smask3);
	__m256 bval = _mm256_blend_ps(_mm256_blend_ps(aa, cc, bmask1), bb, bmask2);
	__m256 gval = _mm256_blend_ps(_mm256_blend_ps(cc, bb, bmask1), aa, bmask2);
	__m256 rval = _mm256_blend_ps(_mm256_blend_ps(bb, aa, bmask1), cc, bmask2);
	_mm256_stream_ps((float*)dst + 0, _mm256_permute2f128_ps(bval, rval, pmask1));
	_mm256_stream_ps((float*)dst + 8, _mm256_permute2f128_ps(gval, bval, pmask2));
	_mm256_stream_ps((float*)dst + 16, _mm256_permute2f128_ps(rval, gval, pmask3));
}

float* createBoxKernel(const int radius)
{
	float* kernel = (float*)_mm_malloc(sizeof(float) * (2 * radius + 1) * (2 * radius + 1), 32);

	int idx = 0;
	const float div = 1.f / ((2 * radius + 1) * (2 * radius + 1));
	for (int i = -radius; i <= radius; i++)
	{
		for (int j = -radius; j <= radius; j++)
		{
			kernel[idx++] = div;
		}
	}

	return kernel;
}

float* createGaussianKernel(const int radius, float sigma)
{
	float* kernel = (float*)_mm_malloc(sizeof(float) * (2 * radius + 1) * (2 * radius + 1), 32);

	int idx = 0;
	float sum = 0.f;
	for (int i = -radius; i <= radius; i++)
	{
		for (int j = -radius; j <= radius; j++)
		{
			const float v = exp((i * i + j * j) / (-2.f * sigma * sigma));
			kernel[idx++] = v;
			sum += v;
		}
	}

	idx = 0;
	sum = 1.f / sum;
	for (int i = -radius; i <= radius; i++)
	{
		for (int j = -radius; j <= radius; j++)
		{
			kernel[idx++] *= sum;
		}
	}

	return kernel;
}

void FIRFilter(const Image_8U& src, Image_8U& dest, const int r, float parameter)
{
	//generate Gaussian kernel
	float* kernel = createGaussianKernel(r, parameter);
	//generate Box kernel
	//float* kernel = createBoxKernel(r);

	Image_8U srcBorder;
	copyMakeBorder(src, srcBorder, r, r, r, r);
	Image_32F temp_32f(srcBorder);

	Image_32F dest_32f(src.rows, src.cols, src.channels);
	const int step = srcBorder.cols;
	if (src.channels == 3)
	{
		Image_8U splitImg[3];

		split(srcBorder, splitImg);

#pragma omp parallel for
		for (int y = 0; y < src.rows; y++)
		{
			unsigned char* rptr = splitImg[0].data + (y + r) * step + r;
			unsigned char* gptr = splitImg[1].data + (y + r) * step + r;
			unsigned char* bptr = splitImg[2].data + (y + r) * step + r;
			float* dptr = dest_32f.data + 3 * (y * dest_32f.cols);

			for (int x = 0; x < src.cols; x += 8)
			{
				__m256 msum_r = _mm256_setzero_ps();
				__m256 msum_g = _mm256_setzero_ps();
				__m256 msum_b = _mm256_setzero_ps();
				int kidx = 0;
				for (int i = -r; i <= r; i++)
				{
					unsigned char* rrptr = rptr + i * step;
					unsigned char* grptr = gptr + i * step;
					unsigned char* brptr = bptr + i * step;
					for (int j = -r; j <= r; j++)
					{
						__m256 mref_r = _mm256_load_epu8cvtps((__m128i*)(rrptr + x + j));
						__m256 mref_g = _mm256_load_epu8cvtps((__m128i*)(grptr + x + j));
						__m256 mref_b = _mm256_load_epu8cvtps((__m128i*)(brptr + x + j));
						__m256 mw = _mm256_set1_ps(kernel[kidx++]);
						msum_r = _mm256_fmadd_ps(mw, mref_r, msum_r);
						msum_g = _mm256_fmadd_ps(mw, mref_g, msum_g);
						msum_b = _mm256_fmadd_ps(mw, mref_b, msum_b);
					}
				}
				_mm256_stream_ps_color(dptr + 3 * x, msum_r, msum_g, msum_b);
			}
		}
	}
	else if (src.channels == 1)
	{
		//gray image
		std::cout << "not support gray image" << std::endl;
		return;
	}
	dest = Image_8U(dest_32f);

	_mm_free(kernel);
}

//////////////////////////////////////////////////////////////////////
// 総合課題（個別）
//////////////////////////////////////////////////////////////////////
///////////////////////
// バイラテラルフィルタ（カラー画像のみ対応）
///////////////////////
void BilateralFilter(const Image_8U& src, Image_8U& dest, const int r, const float sigma_r, const float sigma_s)
{
	Image_8U temp_8u;
	copyMakeBorder(src, temp_8u, r, r, r, r);
	Image_32F temp_32f(temp_8u);

	Image_32F dest_32f(src.rows, src.cols, src.channels);
	if (src.channels == 3)
	{
		Image_32F splitImg[3];
		split(temp_32f, splitImg);
		for (int x = 0; x < src.cols; x++)
		{
			for (int y = 0; y < src.rows; y++)
			{
				float sum[3] = { 0 };
				float wsum = 0;
				const float tgt_r = splitImg[0].data[(y + r) * splitImg[0].cols + (x + r)];
				const float tgt_g = splitImg[1].data[(y + r) * splitImg[1].cols + (x + r)];
				const float tgt_b = splitImg[2].data[(y + r) * splitImg[2].cols + (x + r)];
				for (int i = -r; i <= r; i++)
				{
					for (int j = -r; j <= r; j++)
					{
						const int space_distance = j * j + i * i;
						const float ws = exp(-space_distance / (2.f * sigma_s * sigma_s));

						const float ref_r = splitImg[0].data[(y + r + j) * splitImg[0].cols + (x + r + i)];
						const float ref_g = splitImg[1].data[(y + r + j) * splitImg[1].cols + (x + r + i)];
						const float ref_b = splitImg[2].data[(y + r + j) * splitImg[2].cols + (x + r + i)];
						const float range_distance = (ref_r - tgt_r) * (ref_r - tgt_r) + (ref_g - tgt_g) * (ref_g - tgt_g) + (ref_b - tgt_b) * (ref_b - tgt_b);
						const float wr = exp(-range_distance / (2.f * sigma_r * sigma_r));
						const float w = ws * wr;

						for (int c = 0; c < 3; c++)
						{
							sum[c] += splitImg[c].data[(y + r + j) * splitImg[c].cols + (x + r + i)] * w;
						}
						wsum += w;
					}
				}
				for (int c = 0; c < 3; c++)
				{
					dest_32f.data[3 * (y * dest_32f.cols + x) + c] = sum[c] / wsum + 0.5f;
				}
			}
		}
	}
	else if (src.channels == 1)
	{
		//gray image
		std::cout << "not support gray image" << std::endl;
		return;
	}
	dest = Image_8U(dest_32f);
}


///////////////////////
// ノンローカルミーンフィルタ（グレー画像のみ対応）
///////////////////////
void NonLocalMeansFilter(const Image_8U& src, Image_8U& dest, const int template_r, const int search_r, const float h)
{
	Image_8U temp_8u;
	copyMakeBorder(src, temp_8u, search_r + template_r, search_r + template_r, search_r + template_r, search_r + template_r);
	Image_32F temp_32f(temp_8u);

	const int r = search_r;
	const int tr = template_r;
	const int pad = r + tr;
	Image_32F dest_32f(src.rows, src.cols, src.channels);
	if (src.channels == 1)
	{
		for (int x = 0; x < src.cols; x++)
		{
			for (int y = 0; y < src.rows; y++)
			{
				float sum = 0;
				float wsum = 0;
				for (int i = -r; i <= r; i++)
				{
					for (int j = -r; j <= r; j++)
					{
						float template_distance = 0;
						for (int k = -tr; k <= tr; k++)
						{
							for (int l = -tr; l <= tr; l++)
							{
								const float tgt = temp_32f.data[(y + pad + l) * temp_32f.cols + (x + pad + k)];
								const float ref = temp_32f.data[(y + pad + j + l) * temp_32f.cols + (x + pad + i + k)];
								template_distance += (ref - tgt) * (ref - tgt);
							}
						}
						const float w = exp(-template_distance / (h * h));

						sum += temp_32f.data[(y + pad + j) * temp_32f.cols + (x + pad + i)] * w;
						wsum += w;
					}
				}
				dest_32f.data[y * dest_32f.cols + x] = sum / wsum + 0.5f;
			}
		}
	}
	else if (src.channels == 3)
	{
		//color image
		std::cout << "not support color image" << std::endl;
		return;
	}
	dest = Image_8U(dest_32f);
}
