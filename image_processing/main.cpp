#include "utils/image.h"
#include "utils/image_util.h"
#include "utils/simd_util.h"
#include <time.h>
#include <math.h>

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
	//gaussian filter
	///////////////////////
	//if (false)
	{
		const int loop = 10;

		const float sigma = 2;
		const int r = sigma * 3;
		Image_8U src, temp_8u;
		readPXM("img/lena.ppm", src);
		copyMakeBorder(src, temp_8u, r, r, r, r);
		Image_32F temp_32f(temp_8u);
		Image_32F splitImg[3];
		split(temp_32f, splitImg);
		Image_32F dest_32f(src.rows, src.cols, src.channels);

		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			t.start();
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
			t.end();
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;
		Image_8U dest_8u(dest_32f);
		writePXM("gauss.ppm", dest_8u);
		return 0;
	}
}

