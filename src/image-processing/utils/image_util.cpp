#include "image_util.h"
#include "simd_util.h"
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cmath>

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
	const int size = src.rows * src.cols * src.channels * sizeof(unsigned char);
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

	const int size = src.rows * src.cols * src.channels;
	fwrite(src.data, sizeof(unsigned char), size, fp);
	fprintf(fp, "\n");
	fclose(fp);
}


////////////////////////////
//copy make border
////////////////////////////
void copyMakeBorder(const Image_8U& src, Image_8U& dest, const int top, const int bottom, const int left, const int right)
{
	dest.rows = src.rows + top + bottom;
	dest.cols = src.cols + left + right;
	dest.channels = src.channels;
	if (dest.data != nullptr)
	{
		_mm_free(dest.data);
	}
	const int size = dest.rows * dest.cols * dest.channels * sizeof(unsigned char);
	dest.data = (unsigned char*)_mm_malloc(size, 32);

#pragma omp parallel for
	for (int j = -top; j < src.rows + bottom; j++)
	{
		for (int i = -left; i < src.cols + right; i++)
		{
			const int y = j < 0 ? std::max(j, 0) : std::min(j, src.rows - 1);
			const int x = i < 0 ? std::max(i, 0) : std::min(i, src.cols - 1);
			for (int k = 0; k < src.channels; k++)
			{
				dest.data[dest.channels * ((j + top) * dest.cols + (i + left)) + k] = src.data[src.channels * (y * src.cols + x) + k];
			}
		}
	}
}

void copyMakeBorder(const Image_32F& src, Image_32F& dest, const int top, const int bottom, const int left, const int right)
{
	dest.rows = src.rows + top + bottom;
	dest.cols = src.cols + left + right;
	dest.channels = src.channels;
	if (dest.data != nullptr)
	{
		_mm_free(dest.data);
	}
	const int size = dest.rows * dest.cols * dest.channels * sizeof(float);
	dest.data = (float*)_mm_malloc(size, 32);

	if (src.channels == 3)
	{
		const int LEFT = get_simd_ceil(left, 8);
		const int RIGHT = get_simd_ceil(right, 8);
		const int END = dest.cols - RIGHT;

		const int e = (src.cols - 1) * 3;
		const int sw = src.cols * 3;
		const int dw = dest.cols * 3;
		//src top line
		{
			float* s = src.data;
			float* d = dest.data + top * dw;

			for (int i = 0; i < LEFT; i += 8)
				_mm256_storeu_ps_color(d + 3 * i, _mm256_set1_ps(s[0]), _mm256_set1_ps(s[1]), _mm256_set1_ps(s[2]));

			for (int i = END; i < END + RIGHT; i += 8)
				_mm256_storeu_ps_color(d + 3 * i, _mm256_set1_ps(s[e]), _mm256_set1_ps(s[e + 1]), _mm256_set1_ps(s[e + 2]));

			memcpy(d + 3 * left, s, sizeof(float) * sw);
		}
		//border upper
		for (int j = 0; j < top; j++)
		{
			float* s = dest.data + top * dw;
			float* d = dest.data + j * dw;
			memcpy(d, s, sizeof(float) * dw);
		}
#pragma omp parallel for
		for (int j = top + 1; j < dest.rows - bottom; j++)
		{
			float* s = src.data + (j - top) * sw;
			float* d = dest.data + j * dw;

			for (int i = 0; i < LEFT; i += 8)
				_mm256_storeu_ps_color(d + 3 * i, _mm256_set1_ps(s[0]), _mm256_set1_ps(s[1]), _mm256_set1_ps(s[2]));

			for (int i = END; i < END + RIGHT; i += 8)
				_mm256_storeu_ps_color(d + 3 * i, _mm256_set1_ps(s[e]), _mm256_set1_ps(s[e + 1]), _mm256_set1_ps(s[e + 2]));

			memcpy(d + 3 * left, s, sizeof(float) * sw);
		}

		//border lower
		for (int j = dest.rows - bottom; j < dest.rows; j++)
		{
			float* s = dest.data + (dest.rows - bottom - 1) * dw;
			float* d = dest.data + j * dw;
			memcpy(d, s, sizeof(float) * dw);
		}
	}
	else
	{
#pragma omp parallel for
		for (int j = -top; j < src.rows + bottom; j++)
		{
			for (int i = -left; i < src.cols + right; i++)
			{
				const int y = j < 0 ? std::max(j, 0) : std::min(j, src.rows - 1);
				const int x = i < 0 ? std::max(i, 0) : std::min(i, src.cols - 1);
				for (int k = 0; k < src.channels; k++)
				{
					dest.data[dest.channels * ((j + top) * dest.cols + (i + left)) + k] = src.data[src.channels * (y * src.cols + x) + k];
				}
			}
		}
	}
}


////////////////////////////
//color (RGB) to gray
////////////////////////////
void cvtColorGray(const Image_8U& src, Image_8U& dest)
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
		const int size = dest.rows * dest.cols * dest.channels * sizeof(unsigned char);
		dest.data = (unsigned char*)_mm_malloc(size, 32);
	}

	const int size = src.cols * src.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = (unsigned char)(((src.data[3 * i + 0] + src.data[3 * i + 1] + src.data[3 * i + 2]) * 0.3333f) + 0.5f);
	}
}

void cvtColorGray(const Image_32F& src, Image_32F& dest)
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
		const int size = dest.rows * dest.cols * dest.channels * sizeof(float);
		dest.data = (float*)_mm_malloc(size, 32);
	}

	const int size = src.cols * src.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = (float)(((src.data[3 * i + 0] + src.data[3 * i + 1] + src.data[3 * i + 2]) * 0.3333f) + 0.5f);
	}
}

////////////////////////////
//split image
////////////////////////////
void split(const Image_8U& src, Image_8U* dest)
{
	const int size = src.rows * src.cols * src.channels * sizeof(unsigned char);
	for (int i = 0; i < src.channels; i++)
	{
		dest[i].rows = src.rows;
		dest[i].cols = src.cols;
		dest[i].channels = 1;
		dest[i].type = sizeof(unsigned char);
		dest[i].data = (unsigned char*)_mm_malloc(size, 32);
	}

	const int step = src.cols;
#pragma omp parallel for
	for (int j = 0; j < src.rows; j++)
	{
		for (int i = 0; i < src.cols; i++)
		{
			for (int k = 0; k < src.channels; k++)
			{
				dest[k].data[j * step + i] = src.data[src.channels * (j * step + i) + k];
			}
		}
	}
}

void split(const Image_32F& src, Image_32F* dest)
{
	const int size = src.rows * src.cols * src.channels * sizeof(float);
	for (int i = 0; i < src.channels; i++)
	{
		dest[i].rows = src.rows;
		dest[i].cols = src.cols;
		dest[i].channels = 1;
		dest[i].type = sizeof(float);
		dest[i].data = (float*)_mm_malloc(size, 32);
	}

	const int step = src.cols;
	if (src.channels == 3)
	{
#pragma omp parallel for
		for (int j = 0; j < src.rows; j++)
		{
			float* s = &src.data[src.channels * (j * step)];
			float* d0 = &dest[0].data[j * step];
			float* d1 = &dest[1].data[j * step];
			float* d2 = &dest[2].data[j * step];
			const int simdend = (src.cols / 8) * 8;
			for (int i = 0; i < src.cols; i += 8)
			{
				__m256 b, g, r;
				_mm256_load_cvtps_bgr2planar_ps(s + 3 * i, b, g, r);
				_mm256_storeu_ps(d0 + i, b);
				_mm256_storeu_ps(d1 + i, g);
				_mm256_storeu_ps(d2 + i, r);
			}
			for (int i = simdend; i < src.cols; i++)
			{
				d0[i] = s[3 * i + 0];
				d1[i] = s[3 * i + 1];
				d2[i] = s[3 * i + 2];
			}
		}
	}
	else
	{
#pragma omp parallel for
		for (int j = 0; j < src.rows; j++)
		{
			float* s = &src.data[src.channels * (j * step)];
			float** d = new float* [src.channels];
			for (int k = 0; k < src.channels; k++)
			{
				d[k] = &dest[k].data[j * step];
			}
			for (int i = 0; i < src.cols; i++)
			{
				for (int k = 0; k < src.channels; k++)
				{
					d[k][i] = s[src.channels * i + k];
				}
			}
			delete[]d;
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

#pragma omp parallel for
	for (int j = 0; j < dest.rows; j++)
	{
		for (int i = 0; i < dest.cols; i++)
		{
			for (int k = 0; k < dest.channels; k++)
			{
				dest.data[dest.channels * (j * dest.cols + i) + k] = src[k].data[j * src[k].cols + i];
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

	if (dest.channels == 3)
	{
#pragma omp parallel for
		for (int j = 0; j < dest.rows; j++)
		{
			float* d = &dest.data[dest.channels * (j * dest.cols)];
			float* s0 = &src[0].data[j * src[0].cols];
			float* s1 = &src[1].data[j * src[1].cols];
			float* s2 = &src[2].data[j * src[2].cols];

			const int simdend = (dest.cols / 8) * 8;
			for (int i = 0; i < simdend; i += 8)
			{
				__m256 b, g, r;
				b = _mm256_load_ps(s0 + i);
				g = _mm256_load_ps(s1 + i);
				r = _mm256_load_ps(s2 + i);
				_mm256_stream_ps_color(d + 3 * i, b, g, r);
			}
			for (int i = simdend; i < dest.cols; i++)
			{
				d[3 * i + 0] = s0[i];
				d[3 * i + 1] = s1[i];
				d[3 * i + 2] = s2[i];
			}
		}
	}
	else
	{
#pragma omp parallel for
		for (int j = 0; j < dest.rows; j++)
		{
			float* d = &dest.data[dest.channels * (j * dest.cols)];
			float** s = new float* [dest.channels];
			for (int k = 0; k < dest.channels; k++)
			{
				s[k] = &src[k].data[j * src[k].cols];
			}

			for (int i = 0; i < dest.cols; i++)
			{
				for (int k = 0; k < dest.channels; k++)
				{
					d[dest.channels * i + k] = s[k][i];
				}
			}
			delete[]s;
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
				m.data[m.channels * (j * m.cols + i) + k] = 0;
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
				m.data[m.channels * (j * m.cols + i) + k] = 0;
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
				m.data[m.channels * (j * m.cols + i) + k] = 0;
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
				m.data[m.channels * (j * m.cols + i) + k] = 0;
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
				m.data[m.channels * (j * m.cols + i) + k] = 0;
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
				m.data[m.channels * (j * m.cols + i) + k] = 1;
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
				m.data[m.channels * (j * m.cols + i) + k] = 1;
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
				m.data[m.channels * (j * m.cols + i) + k] = 1;
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
				m.data[m.channels * (j * m.cols + i) + k] = 1.0f;
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
				m.data[m.channels * (j * m.cols + i) + k] = 1.0;
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
				m.data[m.channels * (j * m.cols + i) + k] = rand_min + (rand() * static_cast<float>(rand_max - rand_min) / static_cast<float>(RAND_MAX));
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
				m.data[m.channels * (j * m.cols + i) + k] = rand_min + (rand() * static_cast<float>(rand_max - rand_min) / static_cast<float>(RAND_MAX));
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
				m.data[m.channels * (j * m.cols + i) + k] = rand_min + (rand() * static_cast<float>(rand_max - rand_min) / static_cast<float>(RAND_MAX));
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
				m.data[m.channels * (j * m.cols + i) + k] = rand_min + (rand() * (rand_max - rand_min) / static_cast<float>(RAND_MAX));
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
				m.data[m.channels * (j * m.cols + i) + k] = rand_min + (rand() * (rand_max - rand_min) / static_cast<double>(RAND_MAX));
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
	return rand_min + (rand() * (rand_max - rand_min) / (RAND_MAX));
}
short rand_16s(const short rand_min, const short rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min) / (RAND_MAX));
}
int rand_32s(const int rand_min, const int rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min) / (RAND_MAX));
}
float rand_32f(const float rand_min, const float rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min) / static_cast<float>(RAND_MAX));
}
double rand_64f(const double rand_min, const double rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min) / static_cast<double>(RAND_MAX));
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
	for (int i = 0; i < height * width * ch; i++)
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
	for (int i = 0; i < height * width * ch; i++)
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

//timer
#ifdef USE_TIME_CHRONO
void CalcTime::start()
{
	s = std::chrono::system_clock::now();
	return;
}

void CalcTime::end()
{
	e = std::chrono::system_clock::now();
	double time = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() / 1000.0);
	que.push_back(time); //ms
	return;
}

#elif defined(__APPLE__)
#include <mach/mach_time.h>
#define ORWL_NANO (+1.0E-9)
#define ORWL_GIGA UINT64_C(1000000000)

static double orwl_timebase = 0.0;
static uint64_t orwl_timestart = 0;

struct timespec orwl_gettime(void)
{
	if (!orwl_timestart)
	{
		mach_timebase_info_data_t tb = { 0 };
		mach_timebase_info(&tb);
		orwl_timebase = tb.numer;
		orwl_timebase /= tb.denom;
		orwl_timestart = mach_absolute_time();
	}
	struct timespec t;
	double diff = (mach_absolute_time() - orwl_timestart) * orwl_timebase;
	t.tv_sec = diff * ORWL_NANO;
	t.tv_nsec = diff - (t.tv_sec * ORWL_GIGA);
	return t;
}

void CalcTime::start()
{
	s = orwl_gettime();
	return;
}

void CalcTime::end()
{
	e = orwl_gettime();
	que.push_back((double)(e.tv_sec - s.tv_sec) * 1e3 + (double)(e.tv_nsec - s.tv_nsec) * 1e-6); //ms
	return;
}
#elif defined(_MSC_VER)
CalcTime::CalcTime()
{
	QueryPerformanceFrequency(&frequency);
}

void CalcTime::start()
{
	QueryPerformanceCounter(&s);
	return;
}

void CalcTime::end()
{
	QueryPerformanceCounter(&e);
	LONGLONG span = e.QuadPart - s.QuadPart;
	double sec = (double)span / (double)frequency.QuadPart;

	que.push_back(sec * 1000.0); //msec
	return;
}
#elif defined(__GNUC__)
void CalcTime::start()
{
	clock_gettime(CLOCK_REALTIME, &s);
	return;
}

void CalcTime::end()
{
	clock_gettime(CLOCK_REALTIME, &e);
	que.push_back((double)(e.tv_sec - s.tv_sec) * 1e3 + (double)(e.tv_nsec - s.tv_nsec) * 1e-6); //ms
	return;
}

#endif

void CalcTime::clear()
{
	que.clear();
	return;
}

double CalcTime::getAvgTime(const bool dropFirstMeasure, const bool isClear)
{
	double count = 0;
	double time = 0;
	if ((dropFirstMeasure && que.size() <= 1) || (!dropFirstMeasure && que.size() == 0))
	{
		return -1;
	}
	std::vector<double>::iterator it = que.begin();
	if (dropFirstMeasure)
	{
		it++;
	}
	for (; it != que.end(); ++it)
	{
		time += *it;
		count++;
	}
	if (isClear)
	{
		que.clear();
	}
	return time / count;
}

double CalcTime::getLastTime()
{
	if (que.size() == 0)
	{
		return -1;
	}
	return que.back();
}