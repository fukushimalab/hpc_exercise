#include <stdlib.h>
#include "image_util.h"
#include <iostream>
#include <algorithm>
#ifdef __GNUC__
#include <x86intrin.h>
#elif _MSC_VER
#include <intrin.h>
#pragma warning(disable: 4996)
#endif


//Image util functions
////////////////////////////
//read PXM
////////////////////////////
int readPXM(const char* name, Image_8U& src)
{
	char buffer[256];
	FILE* fp = fopen(name, "rb");
	if (fp == NULL)
	{
		fprintf(stderr, "file open error %s\n", name);
		return -1;
	}

	fscanf(fp, "%s", buffer);
	if (buffer[1] == '3' || buffer[1] == '6')src.channels = 3;//color
	if (buffer[1] == '2' || buffer[1] == '5')src.channels = 1;//binary

	fscanf(fp, "%s", buffer);
	if (buffer[0] == '#')
	{
		fgets(buffer, 256, fp);//skip comment
		fscanf(fp, "%d %d", &src.cols, &src.rows);
	}
	else
	{
		src.cols = atoi(buffer);
		fscanf(fp, "%d", &src.rows);
	}
	fscanf(fp, "%s\n", buffer);//skip intensity

	if (src.data != nullptr)
	{
		_mm_free(src.data);
	}
	const int size = src.rows*src.cols*src.channels * sizeof(unsigned char);
	src.data = (unsigned char*)_mm_malloc(size, 32);
	fread(src.data, sizeof(unsigned char), size, fp);
	fclose(fp);
	return 0;
}

////////////////////////////
//write PXM
////////////////////////////
void writePXM(const char* name, const Image_8U& src)
{
	FILE* fp = fopen(name, "wb");
	if (fp == NULL)
	{
		fprintf(stderr, "file open error %s\n", name);
		return;
	}

	if (src.channels == 1)
	{
		fprintf(fp, "%s\n", "P5");
	}
	else if (src.channels == 3)
	{
		fprintf(fp, "%s\n", "P6");
	}

	fprintf(fp, "# Created by Programming OUYOU\n");
	fprintf(fp, "%d %d\n", src.cols, src.rows);
	fprintf(fp, "255\n");

	const int size = src.rows*src.cols*src.channels;
	fwrite(src.data, sizeof(unsigned char), size, fp);
	fprintf(fp, "\n");
	fclose(fp);
}


////////////////////////////
//copy make border
////////////////////////////
void copyMakeBorder(const Image_8U src, Image_8U& dest, const int top, const int bottom, const int left, const int right)
{
	{
		dest.rows = src.rows + top + bottom;
		dest.cols = src.cols + left + right;
		dest.channels = src.channels;
		if (dest.data != nullptr)
		{
			_mm_free(dest.data);
		}
		const int size = dest.rows*dest.cols*dest.channels * sizeof(unsigned char);
		dest.data = (unsigned char*)_mm_malloc(size, 32);
	}

	for (int j = -top; j < src.rows + bottom; j++)
	{
		for (int i = -left; i < src.cols + right; i++)
		{
			for (int k = 0; k < src.channels; k++)
			{
				const int y = j < 0 ? std::max(j, 0) : std::min(j, src.rows - 1);
				const int x = i < 0 ? std::max(i, 0) : std::min(i, src.cols - 1);
				dest.data[dest.channels*((j + top) * dest.cols + (i + left)) + k] = src.data[src.channels*(y*src.cols + x) + k];
			}
		}
	}
}

void copyMakeBorder(const Image_32F src, Image_32F& dest, const int top, const int bottom, const int left, const int right)
{
	{
		dest.rows = src.rows + top + bottom;
		dest.cols = src.cols + left + right;
		dest.channels = src.channels;
		if (dest.data != nullptr)
		{
			_mm_free(dest.data);
		}
		const int size = dest.rows*dest.cols*dest.channels * sizeof(float);
		dest.data = (float*)_mm_malloc(size, 32);
	}

	for (int j = -top; j < src.rows + bottom; j++)
	{
		for (int i = -left; i < src.cols + right; i++)
		{
			for (int k = 0; k < src.channels; k++)
			{
				const int y = j < 0 ? std::max(j, 0) : std::min(j, src.rows - 1);
				const int x = i < 0 ? std::max(i, 0) : std::min(i, src.cols - 1);
				dest.data[dest.channels*((j + top) * dest.cols + (i + left)) + k] = src.data[src.channels*(y*src.cols + x) + k];
			}
		}
	}
}


////////////////////////////
//color (RGB) to gray
////////////////////////////
void cvtColorGray(const Image_8U src, Image_8U& dest)
{
	if (src.rows != dest.rows || src.cols != dest.cols || dest.channels != 1)
	{
		dest.rows = src.rows;
		dest.cols = src.cols;
		dest.channels = 1;
		if (dest.data != nullptr)
		{
			_mm_free(dest.data);
		}
		const int size = dest.rows*dest.cols*dest.channels * sizeof(unsigned char);
		dest.data = (unsigned char*)_mm_malloc(size, 32);
	}

	const int size = src.cols*src.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = (unsigned char)(((src.data[3 * i + 0] + src.data[3 * i + 1] + src.data[3 * i + 2])*0.3333f) + 0.5f);
	}
}

void cvtColorGray(const Image_32F src, Image_32F& dest)
{
	if (src.rows != dest.rows || src.cols != dest.cols || dest.channels != 1)
	{
		dest.rows = src.rows;
		dest.cols = src.cols;
		dest.channels = 1;
		if (dest.data != nullptr)
		{
			_mm_free(dest.data);
		}
		const int size = dest.rows*dest.cols*dest.channels * sizeof(float);
		dest.data = (float*)_mm_malloc(size, 32);
	}

	const int size = src.cols*src.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = (float)(((src.data[3 * i + 0] + src.data[3 * i + 1] + src.data[3 * i + 2])*0.3333f) + 0.5f);
	}
}

////////////////////////////
//split image
////////////////////////////
void split(const Image_8U src, Image_8U* dest)
{
	const int size = src.rows*src.cols*src.channels * sizeof(unsigned char);
	for (int i = 0; i < src.channels; i++)
	{
		dest[i].rows = src.rows;
		dest[i].cols = src.cols;
		dest[i].channels = 1;
		dest[i].type = sizeof(unsigned char);
		dest[i].data = (unsigned char*)_mm_malloc(size, 32);
	}

	for (int j = 0; j < src.rows; j++)
	{
		for (int i = 0; i < src.cols; i++)
		{
			for (int k = 0; k < src.channels; k++)
			{
				dest[k].data[j*dest[k].cols + i] = src.data[src.channels*(j*src.cols + i) + k];
			}
		}
	}
}

void split(const Image_32F src, Image_32F* dest)
{
	const int size = src.rows*src.cols*src.channels * sizeof(float);
	for (int i = 0; i < src.channels; i++)
	{
		dest[i].rows = src.rows;
		dest[i].cols = src.cols;
		dest[i].channels = 1;
		dest[i].type = sizeof(float);
		dest[i].data = (float*)_mm_malloc(size, 32);
	}

	for (int j = 0; j < src.rows; j++)
	{
		for (int i = 0; i < src.cols; i++)
		{
			for (int k = 0; k < src.channels; k++)
			{
				dest[k].data[j*dest[k].cols + i] = src.data[src.channels*(j*src.cols + i) + k];
			}
		}
	}
}


////////////////////////////
//merge image
////////////////////////////
void merge(const Image_8U* src, const int channel, Image_8U& dest)
{
	dest.rows = src[0].rows;
	dest.cols = src[0].cols;
	dest.type = sizeof(unsigned char);
	dest.channels = channel;
	const int size = dest.rows * dest.cols * dest.channels * dest.type;
	dest.data = (unsigned char*)_mm_malloc(size, 32);

	for (int j = 0; j < dest.rows; j++)
	{
		for (int i = 0; i < dest.cols; i++)
		{
			for (int k = 0; k < dest.channels; k++)
			{
				dest.data[dest.channels*(j*dest.cols + i) + k] = src[k].data[j*src[k].cols + i];
			}
		}
	}
}

void merge(const Image_32F* src, const int channel, Image_32F& dest)
{
	dest.rows = src[0].rows;
	dest.cols = src[0].cols;
	dest.type = sizeof(float);
	dest.channels = channel;
	const int size = dest.rows * dest.cols * dest.channels * dest.type;
	dest.data = (float*)_mm_malloc(size, 32);

	for (int j = 0; j < dest.rows; j++)
	{
		for (int i = 0; i < dest.cols; i++)
		{
			for (int k = 0; k < dest.channels; k++)
			{
				dest.data[dest.channels*(j*dest.cols + i) + k] = src[k].data[j*src[k].cols + i];
			}
		}
	}
}


////////////////////////////
//init (zero)
////////////////////////////
void image_zero(Image_8U& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 0;
			}
		}
	}
}

void image_zero(Image_16S& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 0;
			}
		}
	}
}

void image_zero(Image_32S& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 0;
			}
		}
	}
}

void image_zero(Image_32F& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 0;
			}
		}
	}
}

void image_zero(Image_64F& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 0;
			}
		}
	}
}


////////////////////////////
//init (one)
////////////////////////////
void image_one(Image_8U& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 1;
			}
		}
	}
}

void image_one(Image_16S& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 1;
			}
		}
	}
}

void image_one(Image_32S& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 1;
			}
		}
	}
}

void image_one(Image_32F& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 1.0f;
			}
		}
	}
}

void image_one(Image_64F& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = 1.0;
			}
		}
	}
}


////////////////////////////
//init (rand)
////////////////////////////
void image_rand(Image_8U& m, const unsigned char rand_min, const unsigned char rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = rand_min + (rand() * (rand_max - rand_min + 1) / (1 + RAND_MAX));
			}
		}
	}
}

void image_rand(Image_16S& m, const short rand_min, const short rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = rand_min + (rand() * (rand_max - rand_min + 1) / (1 + RAND_MAX));
			}
		}
	}
}

void image_rand(Image_32S& m, const int rand_min, const int rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = rand_min + (rand() * (rand_max - rand_min + 1) / (1 + RAND_MAX));
			}
		}
	}
}

void image_rand(Image_32F& m, const float rand_min, const float rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = rand_min + (rand() * (rand_max - rand_min + 1.0f) / (1.0f + RAND_MAX));
			}
		}
	}
}

void image_rand(Image_64F& m, const double rand_min, const double rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			for (int k = 0; k < m.channels; k++)
			{
				m.data[m.channels*(j*m.cols + i) + k] = rand_min + (rand() * (rand_max - rand_min + 1.0) / (1.0 + RAND_MAX));
			}
		}
	}
}


////////////////////////////
//image info
////////////////////////////
void image_info(const Image_8U m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " x " << m.channels << " ]  8U" << std::endl;
}

void image_info(const Image_16S m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " x " << m.channels << " ] 16S" << std::endl;
}

void image_info(const Image_32S m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " x " << m.channels << " ] 32S" << std::endl;
}

void image_info(const Image_32F m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " x " << m.channels << " ] 32F" << std::endl;
}

void image_info(const Image_64F m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " x " << m.channels << " ] 64F" << std::endl;
}

//rand
unsigned char rand_8u(const unsigned char rand_min, const unsigned char rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min + 1) / (1 + RAND_MAX));
}

short rand_16s(const short rand_min, const short rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min + 1) / (1 + RAND_MAX));
}

int rand_32s(const int rand_min, const int rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min + 1) / (1 + RAND_MAX));
}

float rand_32f(const float rand_min, const float rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min + 1.0f) / (1.0f + RAND_MAX));
}

double rand_64f(const double rand_min, const double rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min + 1.0) / (1.0 + RAND_MAX));
}


//timer
void CalcTime::start()
{
	clock_gettime(CLOCK_REALTIME, &s);
	return;
}

void CalcTime::end()
{
	clock_gettime(CLOCK_REALTIME, &e);
	que.push_back((double)(e.tv_sec-s.tv_sec)*1e3 + (double)(e.tv_nsec-s.tv_nsec)*1e-6); //ms
	return;
}

void CalcTime::clear()
{
	que.clear();
	return;
}

double CalcTime::getAvgTime(const bool dropFirstMeasure, const bool isClear)
{
	double count = 0;
	double time = 0;
	if((dropFirstMeasure && que.size() <= 1) || (!dropFirstMeasure && que.size() == 0))
	{
		return -1;
	}
	std::vector<double>::iterator it = que.begin();
	if(dropFirstMeasure)
	{
		it++;
	}
	for (; it != que.end(); ++it)
	{
		time += *it;
		count++;
	}
	if(isClear)
	{
		que.clear();
	}
	return time/count;
}

double CalcTime::getLastTime()
{
	if(que.size() == 0)
	{
		return -1;
	}
	return que.back();
}

// calcPSNR
double calcPSNR(const Image_8U& src1, const Image_8U& src2)
{
        if (src1.channels != src2.channels)
        {
                std::cout << "error: different number of chennels" << std::endl;
                return -1;
        }
        if (src1.rows != src2.rows || src1.cols != src2.cols)
        {
                std::cout << "error: different image size" << std::endl;
                return -1;
        }

        Image_64F temp1(src1);
        Image_64F temp2(src2);

        const int height = src1.rows;
        const int width = src1.cols;
        const int ch = src1.channels;

        double sse = 0;
        for (int i = 0; i < height*width*ch; i++)
        {
                double t = src1.data[i] - src2.data[i];
                sse += t * t;
        }

        if (sse <= 1e-10)
        {
                return INFINITY;
        }
        else
        {
                const double  mse = sse / (double)(ch * height * width);
                return 10.0 * log10((255.0 * 255.0) / mse);
        }
}

double calcPSNR(const Image_32F& src1, const Image_32F& src2)
{
        if (src1.channels != src2.channels)
        {
                std::cout << "error: different number of chennels" << std::endl;
                return -1;
        }
        if (src1.rows != src2.rows || src1.cols != src2.cols)
        {
                std::cout << "error: different image size" << std::endl;
                return -1;
        }

        Image_64F temp1(src1);
        Image_64F temp2(src2);

        const int height = src1.rows;
        const int width = src1.cols;
        const int ch = src1.channels;

        double sse = 0;
        for (int i = 0; i < height*width*ch; i++)
        {
                double t = src1.data[i] - src2.data[i];
                sse += t * t;
        }

        if (sse <= 1e-10)
        {
                return INFINITY;
        }
        else
        {
                const double  mse = sse / (double)(ch * height * width);
                return 10.0 * log10((255.0 * 255.0) / mse);
        }
}
