#include "mat_util.h"
#include <stdlib.h>
#include <iostream>

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
