#include "mat_util.h"
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <omp.h>

#ifdef __GNUC__
#include <x86intrin.h>
#elif _MSC_VER
#include <intrin.h>
#endif

//Mat util functions
////////////////////////////
//init (zero)
////////////////////////////
void mat_zero(Mat_8U& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, sizeof(unsigned char) * size, 0);
}

void mat_zero(Mat_16S& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, sizeof(short) * size, 0);
}

void mat_zero(Mat_32S& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, sizeof(int) * size, 0);
}

void mat_zero(Mat_32F& m)
{
	const int size = m.cols * m.rows;
	const int simdsize = (size / 8);
	float* ptr = m.data;
	for (int i = 0; i < simdsize; i++)
	{
		_mm256_store_ps(ptr, _mm256_setzero_ps());
		ptr += 8;
	}
	for (int i = simdsize * 8; i < size; i++)
	{
		m.data[i] = 0.f;
	}
}

void mat_zero(Mat_64F& m)
{
	const int size = m.cols * m.rows;
	const int simdsize = (size / 4);
	double* ptr = m.data;
	for (int i = 0; i < simdsize; i++)
	{
		_mm256_store_pd(ptr, _mm256_setzero_pd());
		ptr += 4;
	}
	for (int i = simdsize * 4; i < size; i++)
	{
		m.data[i] = 0.0;
	}
}


////////////////////////////
//init (one)
////////////////////////////
void mat_one(Mat_8U& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, sizeof(unsigned char) * size, 1);
}

void mat_one(Mat_16S& m)
{
	const int size = m.cols * m.rows;
	const int simdsize = (size / 16);
	short* ptr = m.data;
	for (int i = 0; i < simdsize; i++)
	{
		_mm256_store_si256((__m256i*)ptr, _mm256_set1_epi16(1));
		ptr += 16;
	}
	for (int i = simdsize * 16; i < size; i++)
	{
		m.data[i] = 1;
	}
}

void mat_one(Mat_32S& m)
{
	const int size = m.cols * m.rows;
	const int simdsize = (size / 8);
	int* ptr = m.data;
	for (int i = 0; i < simdsize; i++)
	{
		_mm256_store_si256((__m256i*)ptr, _mm256_set1_epi32(1));
		ptr += 8;
	}
	for (int i = simdsize * 8; i < size; i++)
	{
		m.data[i] = 1;
	}
}

void mat_one(Mat_32F& m)
{
	const int size = m.cols * m.rows;
	const int simdsize = (size / 8);
	float* ptr = m.data;
	for (int i = 0; i < simdsize; i++)
	{
		_mm256_store_ps(ptr, _mm256_set1_ps(1.f));
		ptr += 8;
	}
	for (int i = simdsize * 8; i < size; i++)
	{
		m.data[i] = 1.f;
	}
}

void mat_one(Mat_64F& m)
{
	const int size = m.cols * m.rows;
	const int simdsize = (size / 4);
	double* ptr = m.data;
	for (int i = 0; i < simdsize; i++)
	{
		_mm256_store_pd(ptr, _mm256_set1_pd(1.0));
		ptr += 4;
	}
	for (int i = simdsize * 4; i < size; i++)
	{
		m.data[i] = 1.0;
	}
}


////////////////////////////
//init (rand)
////////////////////////////
void mat_rand(Mat_8U& m, const unsigned char rand_min, const unsigned char rand_max)
{
	const int size = m.rows * m.cols;

	unsigned char* ptr = m.data;
	const float v = (float)(rand_max - rand_min) / (RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (unsigned char)(rand() * v);
	}
}

void mat_rand(Mat_16S& m, const short rand_min, const short rand_max)
{
	const int size = m.rows * m.cols;

	short* ptr = m.data;
	const float v = (float)(rand_max - rand_min) / (RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (short)(rand() * v);
	}
}

void mat_rand(Mat_32S& m, const int rand_min, const int rand_max)
{
	const int size = m.rows * m.cols;

	int* ptr = m.data;
	const float v = (float)(rand_max - rand_min) / (RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (int)(rand() * v);
	}
}

void mat_rand(Mat_32F& m, const float rand_min, const float rand_max)
{
	const int size = m.rows * m.cols;

	float* ptr = m.data;
	const float v = (float)(rand_max - rand_min) / (RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (rand() * v);
	}
}

void mat_rand(Mat_64F& m, const double rand_min, const double rand_max)
{
	const int size = m.rows * m.cols;

	double* ptr = m.data;
	const double v = (double)(rand_max - rand_min) / (RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (rand() * v);
	}
}

////////////////////////////
//val add 
////////////////////////////
Mat_8U mat_add(const Mat_8U& m, const unsigned char v)
{
	Mat_8U dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] + v;
	}
	return dest;
}

Mat_16S mat_add(const Mat_16S& m, const short v)
{
	Mat_16S dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] + v;
	}
	return dest;
}

Mat_32S mat_add(const Mat_32S& m, const int v)
{
	Mat_32S dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] + v;
	}
	return dest;
}

Mat_32F mat_add(const Mat_32F& m, const float v)
{
	Mat_32F dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] + v;
	}
	return dest;
}

Mat_64F mat_add(const Mat_64F& m, const double v)
{
	Mat_64F dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] + v;
	}
	return dest;
}


////////////////////////////
//mat add
////////////////////////////
Mat_8U mat_add(const Mat_8U& m1, const Mat_8U& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_8U dest(m1.rows, m1.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
	return dest;
}

Mat_16S mat_add(const Mat_16S& m1, const Mat_16S& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_16S dest(m1.rows, m1.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
	return dest;
}

Mat_32S mat_add(const Mat_32S& m1, const Mat_32S& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_32S dest(m1.rows, m1.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
	return dest;
}

Mat_32F mat_add(const Mat_32F& m1, const Mat_32F& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_32F dest(m1.rows, m1.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
	return dest;
}

Mat_64F mat_add(const Mat_64F& m1, const Mat_64F& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_64F dest(m1.rows, m1.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
	return dest;
}

////////////////////////////
//val mul
////////////////////////////
Mat_8U mat_mul(const Mat_8U& m, const unsigned char v)
{
	Mat_8U dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] * v;
	}
	return dest;
}

Mat_16S mat_mul(const Mat_16S& m, const short v)
{
	Mat_16S dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] * v;
	}
	return dest;
}

Mat_32S mat_mul(const Mat_32S& m, const int v)
{
	Mat_32S dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] * v;
	}
	return dest;
}

Mat_32F mat_mul(const Mat_32F& m, const float v)
{
	Mat_32F dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] * v;
	}
	return dest;
}

Mat_64F mat_mul(const Mat_64F& m, const double v)
{
	Mat_64F dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] * v;
	}
	return dest;
}


////////////////////////////
//mat mul
////////////////////////////
Mat_8U mat_mul(const Mat_8U& m1, const Mat_8U& m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_8U dest(m1.rows, m1.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] * m2.data[i];
	}
	return dest;
}

Mat_16S mat_mul(const Mat_16S& m1, const Mat_16S& m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_16S dest(m1.rows, m1.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] * m2.data[i];
	}
	return dest;
}

Mat_32S mat_mul(const Mat_32S& m1, const Mat_32S& m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_32S dest(m1.rows, m1.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] * m2.data[i];
	}
	return dest;
}

Mat_32F mat_mul(const Mat_32F& m1, const Mat_32F& m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_32F dest(m1.rows, m2.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] * m2.data[i];
	}
	return dest;
}

Mat_64F mat_mul(const Mat_64F& m1, const Mat_64F& m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_64F dest(m1.rows, m2.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] * m2.data[i];
	}
	return dest;
}


////////////////////////////
//val div
////////////////////////////
Mat_8U mat_div(const Mat_8U& m, const unsigned char v)
{
	Mat_8U dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] / v;
	}
	return dest;
}

Mat_16S mat_div(const Mat_16S& m, const short v)
{
	Mat_16S dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] / v;
	}
	return dest;
}

Mat_32S mat_div(const Mat_32S& m, const int v)
{
	Mat_32S dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] / v;
	}
	return dest;
}

Mat_32F mat_div(const Mat_32F& m, const float v)
{
	Mat_32F dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] / v;
	}
	return dest;
}

Mat_64F mat_div(const Mat_64F& m, const double v)
{
	Mat_64F dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] / v;
	}
	return dest;
}


////////////////////////////
//show
////////////////////////////
void mat_show(const Mat_8U& m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 8U" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << (int)m.data[j * m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

void mat_show(const Mat_16S& m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 16S" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << m.data[j * m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

void mat_show(const Mat_32S& m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 32S" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << m.data[j * m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

void mat_show(const Mat_32F& m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 32F" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << m.data[j * m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

void mat_show(const Mat_64F& m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 64F" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << m.data[j * m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
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
	return rand_min + (rand() * (rand_max - rand_min) / (RAND_MAX));
}

double rand_64f(const double rand_min, const double rand_max)
{
	return rand_min + (rand() * (rand_max - rand_min) / (RAND_MAX));
}

int mat_diff(Mat_8U& src1, Mat_8U& src2)
{
	if (src1.rows != src2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	int ret = 0;
	const int size = src1.cols * src1.rows;

	for (int i = 0; i < size; i++)
	{
		ret += ((int)(src1.data[i] - src2.data[i]) * (int)(src1.data[i] - src2.data[i]));
	}
	return ret;
}

int mat_diff(Mat_16S& src1, Mat_16S& src2)
{
	if (src1.rows != src2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	int ret = 0;
	const int size = src1.cols * src1.rows;

	for (int i = 0; i < size; i++)
	{
		ret += ((int)(src1.data[i] - src2.data[i]) * (int)(src1.data[i] - src2.data[i]));
	}
	return ret;
}
int mat_diff(Mat_32S& src1, Mat_32S& src2)
{
	if (src1.rows != src2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	int ret = 0;
	const int size = src1.cols * src1.rows;

	for (int i = 0; i < size; i++)
	{
		ret += ((src1.data[i] - src2.data[i]) * (src1.data[i] - src2.data[i]));
	}
	return ret;
}

double mat_diff(Mat_32F& src1, Mat_32F& src2)
{
	if (src1.rows != src2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	double ret = 0.0;
	const int size = src1.cols * src1.rows;

	for (int i = 0; i < size; i++)
	{
		ret += (double)((src1.data[i] - src2.data[i]) * (src1.data[i] - src2.data[i]));
	}
	return ret;
}

double mat_diff(Mat_64F& src1, Mat_64F& src2)
{
	if (src1.rows != src2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	double ret = 0.0;
	const int size = src1.cols * src1.rows;

	for (int i = 0; i < size; i++)
	{
		ret += ((src1.data[i] - src2.data[i]) * (src1.data[i] - src2.data[i]));
	}
	return ret;
}

//timer
#ifdef __GNUC__
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
#elif _MSC_VER

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
#endif

void loofline_test_omp(const int size, const int iteration, const int num_thread)
{
	const int thread_max = (num_thread == -1) ? omp_get_max_threads() : num_thread;
	omp_set_num_threads(thread_max);

	//FLOPSŒvŽZ
	//x+1: (1*x+1)
	//x*x+x+1: x*(x+1)+1
	//x*x*x+x*x+x+1: x*(x*(x+1))+1

	std::cout << "loofline test omp" << std::endl;
	const int loop = iteration;

	CalcTime t;

	float* x = (float*)_mm_malloc(sizeof(float) * size, 32);
	float* y = (float*)_mm_malloc(sizeof(float) * size * thread_max, 32);

	float* ptr = x;
	const float rand_max = 1.f;
	const float rand_min = 0.f;
	const float v = (float)(rand_max - rand_min) / (RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (rand() * v);
	}

	int simdsize8 = size / 8;
	int simdsize16 = size / 16;

	printf("size %d, iteration\n", iteration);
	printf("order, GFLOPS, FLOPS/BYTE\n");

	int n = 0;
	{
		//---------------------------------------------
		n = 1;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			int v = omp_get_thread_num();
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 2;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}

		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 3;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 4;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 5;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 6;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 7;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 8;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 9;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 10;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 11;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 12;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 13;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 14;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 15;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 16;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 17;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 18;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 19;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18
				ret = _mm256_fmadd_ps(mx, ret, mones);//19
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//19

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 20;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18
				ret = _mm256_fmadd_ps(mx, ret, mones);//19
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//19
				ret = _mm256_fmadd_ps(mx, ret, mones);//20
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//20

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}

		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	}

	_mm_free(x);
	_mm_free(y);
}

void loofline_test_omp2(const int size, const int iteration, const int num_thread)
{
	const int thread_max = (num_thread == -1) ? omp_get_max_threads() : num_thread;
	omp_set_num_threads(thread_max);

	//FLOPSŒvŽZ
	//x+1: (1*x+1)
	//x*x+x+1: x*(x+1)+1
	//x*x*x+x*x+x+1: x*(x*(x+1))+1

	std::cout << "loofline test" << std::endl;
	const int loop = iteration;
	int simdsize8 = size / 8;
	int simdsize16 = size / 16;

	printf("size %d, iteration\n", iteration);
	printf("order, GFLOPS, FLOPS/BYTE\n");

#pragma omp parallel for
	for (int th = 0; th < thread_max; th++)
	{
		CalcTime t;
		float* x = (float*)_mm_malloc(sizeof(float) * size, 32);
		float* y = (float*)_mm_malloc(sizeof(float) * size, 32);

		float* ptr = x;
		const float rand_max = 1.f;
		const float rand_min = 0.f;
		const float v = (float)(rand_max - rand_min) / (RAND_MAX);
		for (int i = 0; i < size; i++)
		{
			*ptr++ = rand_min + (rand() * v);
		}

		int n = 0;
		//---------------------------------------------
		n = 1;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 2;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 3;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 4;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 5;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 6;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 7;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 8;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 9;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 10;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 11;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 12;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 13;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 14;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 15;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 16;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 17;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 18;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 19;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18
				ret = _mm256_fmadd_ps(mx, ret, mones);//19
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//19

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 20;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18
				ret = _mm256_fmadd_ps(mx, ret, mones);//19
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//19
				ret = _mm256_fmadd_ps(mx, ret, mones);//20
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//20

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		if (omp_get_thread_num() == 0)printf("%02d, %f, %f\n", n, thread_max * n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
		_mm_free(x);
		_mm_free(y);
	}
}

void loofline_test(const int size, const int iteration)
{
	//FLOPSŒvŽZ
	//x+1: (1*x+1)
	//x*x+x+1: x*(x+1)+1
	//x*x*x+x*x+x+1: x*(x*(x+1))+1

	std::cout << "loofline test" << std::endl;
	const int loop = iteration;

	CalcTime t;

	float* x = (float*)_mm_malloc(sizeof(float) * size, 32);
	float* y = (float*)_mm_malloc(sizeof(float) * size, 32);

	float* ptr = x;
	const float rand_max = 1.f;
	const float rand_min = 0.f;
	const float v = (float)(rand_max - rand_min) / (RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (rand() * v);
	}

	int simdsize8 = size / 8;
	int simdsize16 = size / 16;

	printf("size %d, iteration\n", iteration);
	printf("order, GFLOPS, FLOPS/BYTE\n");

	int n = 0;

	bool isUnroll2 = true;
	if (isUnroll2)
	{
		//---------------------------------------------
		n = 1;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 2;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 3;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 4;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 5;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 6;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 7;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 8;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 9;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 10;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 11;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 12;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 13;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 14;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 15;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 16;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 17;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 18;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 19;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18
				ret = _mm256_fmadd_ps(mx, ret, mones);//19
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//19

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 20;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize16; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18
				ret = _mm256_fmadd_ps(mx, ret, mones);//19
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//19
				ret = _mm256_fmadd_ps(mx, ret, mones);//20
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//20

				_mm256_store_ps(py, ret);
				_mm256_store_ps(py + 8, ret2);

				px += 16;
				py += 16;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	}
	else
	{
		//---------------------------------------------
		n = 1;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 2;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 3;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 4;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 5;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 6;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 7;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 8;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 9;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 10;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 11;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 12;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 13;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 14;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 15;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 16;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 17;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 18;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 19;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret = _mm256_fmadd_ps(mx, ret, mones);//19

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

		//---------------------------------------------
		n = 20;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x;
			float* py = y;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize8; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret = _mm256_fmadd_ps(mx, ret, mones);//11
				ret = _mm256_fmadd_ps(mx, ret, mones);//12
				ret = _mm256_fmadd_ps(mx, ret, mones);//13
				ret = _mm256_fmadd_ps(mx, ret, mones);//14
				ret = _mm256_fmadd_ps(mx, ret, mones);//15
				ret = _mm256_fmadd_ps(mx, ret, mones);//16
				ret = _mm256_fmadd_ps(mx, ret, mones);//17
				ret = _mm256_fmadd_ps(mx, ret, mones);//18
				ret = _mm256_fmadd_ps(mx, ret, mones);//19
				ret = _mm256_fmadd_ps(mx, ret, mones);//20

				_mm256_store_ps(py, ret);

				px += 8;
				py += 8;
			}
			t.end();
		}
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	}

	_mm_free(x);
	_mm_free(y);
}