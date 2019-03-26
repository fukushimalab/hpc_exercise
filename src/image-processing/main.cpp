#include "utils/image.h"
#include "utils/image_util.h"
#include "utils/simd_util.h"
#include <time.h>
#include <math.h>

/////////////////////////////////
// チュートリアル
/////////////////////////////////
void GammaCorrection(const Image_8U& src, Image_8U& dest, const float gamma);
void MeanVar(const Image_8U& src, float& mean, float& var);
void GaussianFilter(const Image_8U& src, Image_8U& dest, const int r, const float sigma);

/////////////////////////////////
// 総合課題（個別）
/////////////////////////////////
void BilateralFilter(const Image_8U& src, Image_8U& dest, const int r, const float sigma_r, const float sigma_s);
void NonLocalMeansFilter(const Image_8U& src, Image_8U& dest, const int template_r, const int search_r, const float h);


int main(const int argc, const char** argv)
{
	///////////////////////
	//image read
	///////////////////////
	//Image_8U img;
	//readPXM("img/lena.ppm", img);
	///////////////////////
	//copyMakeBorder
	///////////////////////
	//Image_8U cpyImg;
	//const int top = 100;
	//const int bottom = 100;
	//const int left = 100;
	//const int right = 100;
	//copyMakeBorder(img, cpyImg, top, bottom, left, right);
	//writePXM("lena_copyMakeBorder.ppm", cpyImg);
	///////////////////////
	//color2gray
	///////////////////////
	//Image_8U gray;
	//cvtColorGray(img, gray);
	//writePXM("lena_gray.pgm", gray);
	///////////////////////
	//split and merge
	///////////////////////
	//Image_8U splitImg[3];
	//split(img, splitImg);
	//writePXM("lena_r.pgm", splitImg[0]);
	//writePXM("lena_g.pgm", splitImg[1]);
	//writePXM("lena_b.pgm", splitImg[2]);
	//Image_8U mergeImg;
	//merge(splitImg, 3, mergeImg);
	//writePXM("lena_merge.ppm", mergeImg);


	//////////////////////////////////////////////////////////////////////
	// チュートリアル
	//////////////////////////////////////////////////////////////////////
	///////////////////////
	// gamma correction
	///////////////////////
	if (false)
	{
		const int loop = 10;

		const float gamma = 2;
		Image_8U src, dest;
		readPXM("img/lena.ppm", src);

		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			GammaCorrection(src, dest, gamma);
			t.end();
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;
		writePXM("gamma.ppm", dest);
		return 0;
	}

	///////////////////////
	// mean-var
	///////////////////////
	if (false)
	{
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
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "mean: " << mean << std::endl;
		std::cout << "var : " << var << std::endl;
		return 0;
	}

	///////////////////////
	// Gaussian filter
	///////////////////////
	if (false)
	{
		const int loop = 10;

		const float sigma = 2;
		const int r = sigma * 3;
		Image_8U src, dest;
		readPXM("img/lena.ppm", src);


		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			GaussianFilter(src, dest, r, sigma);
			t.end();
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;
		writePXM("gauss.ppm", dest);
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
		const int r = 3 * sigma_s;
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
		for (int x = 0; x < src.cols*cn; x++)
		{
			dest.data[cn*(y*src.cols) + x] = pow((float)src.data[cn*(y*src.cols) + x] / 255.f, 1.f / gamma)*255.0f;
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
		for (int x = 0; x < src.cols*cn; x++)
		{
			mean += src.data[cn*(y*src.cols) + x];
			var += src.data[cn*(y*src.cols) + x] * src.data[cn*(y*src.cols) + x];
		}
	}
	mean /= (float)(src.rows*src.cols*cn);
	var = var / (float)(src.rows*src.cols*cn) - mean * mean;
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
						const float w = exp(-distance / (2.f*sigma*sigma));

						for (int c = 0; c < 3; c++)
						{
							sum[c] += splitImg[c].data[(y + r + j)*splitImg[c].cols + (x + r + i)] * w;
						}
						wsum += w;
					}
				}
				for (int c = 0; c < 3; c++)
				{
					dest_32f.data[3 * (y*dest_32f.cols + x) + c] = sum[c] / wsum;
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
				const float tgt_r = splitImg[0].data[(y + r)*splitImg[0].cols + (x + r)];
				const float tgt_g = splitImg[1].data[(y + r)*splitImg[1].cols + (x + r)];
				const float tgt_b = splitImg[2].data[(y + r)*splitImg[2].cols + (x + r)];
				for (int i = -r; i <= r; i++)
				{
					for (int j = -r; j <= r; j++)
					{
						const int space_distance = j * j + i * i;
						const float ws = exp(-space_distance / (2.f*sigma_s*sigma_s));

						const float ref_r = splitImg[0].data[(y + r + j)*splitImg[0].cols + (x + r + i)];
						const float ref_g = splitImg[1].data[(y + r + j)*splitImg[1].cols + (x + r + i)];
						const float ref_b = splitImg[2].data[(y + r + j)*splitImg[2].cols + (x + r + i)];
						const float range_distance = (ref_r - tgt_r)*(ref_r - tgt_r) + (ref_g - tgt_g)*(ref_g - tgt_g) + (ref_b - tgt_b)*(ref_b - tgt_b);
						const float wr = exp(-range_distance / (2.f*sigma_r*sigma_r));
						const float w = ws * wr;

						for (int c = 0; c < 3; c++)
						{
							sum[c] += splitImg[c].data[(y + r + j)*splitImg[c].cols + (x + r + i)] * w;
						}
						wsum += w;
					}
				}
				for (int c = 0; c < 3; c++)
				{
					dest_32f.data[3 * (y*dest_32f.cols + x) + c] = sum[c] / wsum + 0.5f;
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
								const float tgt = temp_32f.data[(y + pad + l)*temp_32f.cols + (x + pad + k)];
								const float ref = temp_32f.data[(y + pad + j + l)*temp_32f.cols + (x + pad + i + k)];
								template_distance += (ref - tgt)*(ref - tgt);
							}
						}
						const float w = exp(-template_distance / (h*h));

						sum += temp_32f.data[(y + pad + j)*temp_32f.cols + (x + pad + i)] * w;
						wsum += w;
					}
				}
				dest_32f.data[y*dest_32f.cols + x] = sum / wsum + 0.5f;
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

