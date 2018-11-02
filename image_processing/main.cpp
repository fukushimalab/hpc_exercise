#include "utils/image.h"
#include "utils/image_util.h"
#include "utils/simd_util.h"
#include <time.h>
#include <math.h>

void GammaCorrection(const Image_8U& src, Image_8U& dest, const float gamma);
void MeanVar(const Image_8U& src, float& mean, float& var);
void GaussianFilter(const Image_8U& src, Image_8U& dest, const int r, const float sigma);

int main(const int argc, const char** argv)
{
	//image read
	//Image_8U img;
	//readPXM("img/lena.ppm", img);

	//copyMakeBorder
	//Image_8U cpyImg;
	//const int top = 100;
	//const int bottom = 100;
	//const int left = 100;
	//const int right = 100;
	//copyMakeBorder(img, cpyImg, top, bottom, left, right);
	//writePXM("lena1.ppm", cpyImg);

	//color2gray
	//Image_8U gray;
	//cvtColorGray(img, gray);

	//split and merge
	//writePXM("lena1.pgm", gray);
	//Image_8U splitImg[3];
	//split(img, splitImg);
	//writePXM("lena1.pgm", splitImg[0]);
	//writePXM("lena2.pgm", splitImg[1]);
	//writePXM("lena3.pgm", splitImg[2]);
	//Image_8U mergeImg;
	//merge(splitImg, 3, mergeImg);
	//writePXM("lena1.ppm", mergeImg);


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
}


void GammaCorrection(const Image_8U& src, Image_8U& dest, const float gamma)
{
	dest = Image_8U(src.rows, src.cols, src.channels);
	const int cn = src.channels;
	for (int y = 0; y < src.rows; y++)
	{
		for (int x = 0; x < src.cols*cn; x++)
		{
			dest.data[cn*(y*src.cols)+x] = pow((float)src.data[cn*(y*src.cols)+x]/255.f, 1.f/gamma)*255.0f;
		}
	}
}


void MeanVar(const Image_8U& src, float& mean, float& var)
{
	mean = 0;
	var = 0;
	const int cn = src.channels;

	for (int y = 0; y < src.rows; y++)
	{
		for (int x = 0; x < src.cols*cn; x++)
		{
			mean += src.data[cn*(y*src.cols)+x];
			var += src.data[cn*(y*src.cols)+x]*src.data[cn*(y*src.cols)+x];
		}
	}
	mean /=(float)(src.rows*src.cols*cn);
	var = var/(float)(src.rows*src.cols*cn) - mean*mean;
}

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
	else if(src.channels == 1)
	{
			//gray image
	}
	dest = Image_8U(dest_32f);
}

