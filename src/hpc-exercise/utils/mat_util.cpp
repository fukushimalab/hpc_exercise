#include "mat_util.h"
#include "simd_util.h"
#include <iostream>
#include <cstring>

#include <immintrin.h>
#include <pmmintrin.h>
#include <bitset>

//Mat util functions

#pragma region init(zero)
void mat_zero(Mat_8S& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, 0, sizeof(char) * size);
}
void mat_zero(Mat_8U& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, 0, sizeof(unsigned char) * size);
}
void mat_zero(Mat_16S& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, 0, sizeof(short) * size);
}
void mat_zero(Mat_32S& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, 0, sizeof(int) * size);
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
#pragma endregion

#pragma region init(one)
void mat_one(Mat_8S& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, 1, sizeof(char) * size);
}
void mat_one(Mat_8U& m)
{
	const int size = m.cols * m.rows;
	memset(m.data, 1, sizeof(unsigned char) * size);
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
#pragma endregion

#pragma region init(rand)
void mat_rand(Mat_8S& m, const char rand_min, const char rand_max)
{
	const int size = m.rows * m.cols;

	char* ptr = m.data;
	const float v = (float)(rand_max - rand_min) / static_cast<float>(RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (char)(rand() * v);
	}
}

void mat_rand(Mat_8U& m, const unsigned char rand_min, const unsigned char rand_max)
{
	const int size = m.rows * m.cols;

	unsigned char* ptr = m.data;
	const float v = (float)(rand_max - rand_min) / static_cast<float>(RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (unsigned char)(rand() * v);
	}
}
void mat_rand(Mat_16S& m, const short rand_min, const short rand_max)
{
	const int size = m.rows * m.cols;

	short* ptr = m.data;
	const float v = (float)(rand_max - rand_min) / static_cast<float>(RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (short)(rand() * v);
	}
}
void mat_rand(Mat_32S& m, const int rand_min, const int rand_max)
{
	const int size = m.rows * m.cols;

	int* ptr = m.data;
	const float v = (float)(rand_max - rand_min) / static_cast<float>(RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (int)(rand() * v);
	}
}
void mat_rand(Mat_32F& m, const float rand_min, const float rand_max)
{
	const int size = m.rows * m.cols;

	float* ptr = m.data;
	const float v = (float)(rand_max - rand_min) / static_cast<float>(RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (rand() * v);
	}
}
void mat_rand(Mat_64F& m, const double rand_min, const double rand_max)
{
	const int size = m.rows * m.cols;

	double* ptr = m.data;
	const double v = (double)(rand_max - rand_min) / static_cast<double>(RAND_MAX);
	for (int i = 0; i < size; i++)
	{
		*ptr++ = rand_min + (rand() * v);
	}
}
#pragma endregion


#pragma region valAdd
Mat_8S mat_add(const Mat_8S& m, const char v)
{
	Mat_8S dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] + v;
	}
	return dest;
}
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
#pragma endregion

#pragma region matAdd
void mat_add_scalar(const Mat_8S& m1, const Mat_8S& m2, Mat_8S& dest)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
}
void mat_add(const Mat_8S& m1, const Mat_8S& m2, Mat_8S& dest)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	const int size = m1.cols * m1.rows;
	const int simdsize = get_simd_floor(size, 32);
	for (int i = 0; i < simdsize; i += 32)
	{
		__m256i ms1 = _mm256_load_si256((__m256i*)(m1.data + i));
		__m256i ms2 = _mm256_load_si256((__m256i*)(m2.data + i));
		_mm256_store_si256((__m256i*)(dest.data + i), _mm256_add_epi8(ms1, ms2));
	}
	for (int i = simdsize; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
}
Mat_8S mat_add(const Mat_8S& m1, const Mat_8S& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_8S dest(m1.rows, m1.cols);
	mat_add(m1, m2, dest);

	return dest;
}
void mat_add_scalar(const Mat_8U& m1, const Mat_8U& m2, Mat_8U& dest)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
}
void mat_add(const Mat_8U& m1, const Mat_8U& m2, Mat_8U& dest)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	const int size = m1.cols * m1.rows;
	const int simdsize = get_simd_floor(size, 32);
	for (int i = 0; i < simdsize; i += 32)
	{
		__m256i ms1 = _mm256_load_si256((__m256i*)(m1.data + i));
		__m256i ms2 = _mm256_load_si256((__m256i*)(m2.data + i));
		_mm256_store_si256((__m256i*)(dest.data + i), _mm256_adds_epu8(ms1, ms2));
	}
	for (int i = simdsize; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
}
Mat_8U mat_add(const Mat_8U& m1, const Mat_8U& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_8U ret(m1.rows, m1.cols);
	mat_add(m1, m2, ret);

	return ret;
}
void mat_add(const Mat_16S& m1, const Mat_16S& m2, Mat_16S& dest)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
}
Mat_16S mat_add(const Mat_16S& m1, const Mat_16S& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_16S dest(m1.rows, m1.cols);
	mat_add(m1, m2, dest);

	return dest;
}
void mat_add(const Mat_32S& m1, const Mat_32S& m2, Mat_32S& dest)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
}
Mat_32S mat_add(const Mat_32S& m1, const Mat_32S& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_32S dest(m1.rows, m1.cols);
	mat_add(m1, m2, dest);

	return dest;
}
void mat_add(const Mat_32F& m1, const Mat_32F& m2, Mat_32F& dest)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
}
Mat_32F mat_add(const Mat_32F& m1, const Mat_32F& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_32F dest(m1.rows, m1.cols);
	mat_add(m1, m2, dest);

	return dest;
}
void mat_add(const Mat_64F& m1, const Mat_64F& m2, Mat_64F& dest)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] + m2.data[i];
	}
}
Mat_64F mat_add(const Mat_64F& m1, const Mat_64F& m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_64F dest(m1.rows, m1.cols);
	mat_add(m1, m2, dest);

	return dest;
}
#pragma endregion

#pragma region valMul
Mat_8S mat_mul(const Mat_8S& m, const char v)
{
	Mat_8S dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] * v;
	}
	return dest;
}

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
#pragma endregion

#pragma region matMul
Mat_8S mat_mul(const Mat_8S& m1, const Mat_8S& m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_8S dest(m1.rows, m1.cols);
	const int size = m1.cols * m1.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m1.data[i] * m2.data[i];
	}
	return dest;
}
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
	const int simdsize = get_simd_floor(size, 16);
	for (int i = 0; i < simdsize; i += 16)
	{
		__m256i ms1 = _mm256_load_si256((__m256i*)(m1.data + i));
		__m256i ms2 = _mm256_load_si256((__m256i*)(m2.data + i));
		_mm256_store_si256((__m256i*)(dest.data + i), _mm256_mullo_epi16(ms1, ms2));
	}
	for (int i = simdsize; i < size; i++)
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
	const int simdsize = get_simd_floor(size, 8);
	for (int i = 0; i < simdsize; i += 8)
	{
		__m256i ms1 = _mm256_load_si256((__m256i*)(m1.data + i));
		__m256i ms2 = _mm256_load_si256((__m256i*)(m2.data + i));
		_mm256_store_si256((__m256i*)(dest.data + i), _mm256_mullo_epi32(ms1, ms2));
	}
	for (int i = simdsize; i < size; i++)
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
	const int simdsize = get_simd_floor(size, 8);
	for (int i = 0; i < simdsize; i += 8)
	{
		__m256 ms1 = _mm256_load_ps(m1.data + i);
		__m256 ms2 = _mm256_load_ps(m2.data + i);
		_mm256_store_ps((dest.data + i), _mm256_mul_ps(ms1, ms2));
	}
	for (int i = simdsize; i < size; i++)
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
	const int simdsize = get_simd_floor(size, 4);
	for (int i = 0; i < simdsize; i += 4)
	{
		__m256d ms1 = _mm256_load_pd(m1.data + i);
		__m256d ms2 = _mm256_load_pd(m2.data + i);
		_mm256_store_pd((dest.data + i), _mm256_mul_pd(ms1, ms2));
	}
	for (int i = simdsize; i < size; i++)
	{
		dest.data[i] = m1.data[i] * m2.data[i];
	}
	return dest;
}
#pragma endregion

#pragma region valDiv
Mat_8S mat_div(const Mat_8S& m, const unsigned char v)
{
	Mat_8S dest(m.rows, m.cols);
	const int size = m.cols * m.rows;
	for (int i = 0; i < size; i++)
	{
		dest.data[i] = m.data[i] / v;
	}
	return dest;
}
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
#pragma endregion


#pragma region show
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
#pragma endregion

#pragma region rand
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
#pragma endregion

#pragma region matDiff
int mat_diff(Mat_8U& src1, Mat_8U& src2)
{
	if (src1.rows * src1.cols != src2.cols * src2.rows)
	{
		std::cout << "invalid mat size (mat diff)" << std::endl;
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
	if (src1.rows * src1.cols != src2.cols * src2.rows)
	{
		std::cout << "invalid mat size (mat diff)" << std::endl;
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
	if (src1.rows * src1.cols != src2.cols * src2.rows)
	{
		std::cout << "invalid mat size (mat diff)" << std::endl;
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
	if (src1.rows * src1.cols != src2.cols * src2.rows)
	{
		std::cout << "invalid mat size (mat diff)" << std::endl;
		exit(-1);
	}

	double ret = 0.0;
	const int size = src1.cols * src1.rows;

	for (int i = 0; i < size; i++)
	{
		ret += ((double)(src1.data[i] - src2.data[i]) * (double)(src1.data[i] - src2.data[i]));
	}
	return ret;
}
double mat_diff(Mat_64F& src1, Mat_64F& src2)
{
	if (src1.rows * src1.cols != src2.cols * src2.rows)
	{
		std::cout << "invalid mat size (mat diff)" << std::endl;
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
#pragma endregion


#pragma region timer
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
#pragma endregion

bool show_mxcsr(const bool showState, const bool showMask, const bool isClaer)
{
	bool ret = false;

	if (isClaer)
	{
		_MM_SET_EXCEPTION_STATE(0);
	}

	const unsigned int mxcsr = _mm_getcsr();
	std::cout << mxcsr << std::endl;
	std::bitset<16> bit(mxcsr);
	for (int i = 0; i < bit.size(); i++)
	{
		std::cout << bit[i];
	}
	std::cout << std::endl;

	if (mxcsr & 0b10)
	{
		ret = true;
	}

	if (showState)
	{
		if (mxcsr & 0b1)
		{
			std::cout << "Invalid Operation happens" << std::endl;
		}
		if (mxcsr & 0b10)
		{
			std::cout << "Denormal happens" << std::endl;
		}
		if (mxcsr & 0b100)
		{
			std::cout << "Divede By Zero happens" << std::endl;
		}
		if (mxcsr & 0b1000)
		{
			std::cout << "Overflow happens" << std::endl;
		}
		if (mxcsr & 0b10000)
		{
			std::cout << "Underflow happens" << std::endl;
		}
		if (mxcsr & 0b100000)
		{
			std::cout << "Precision happens" << std::endl;
		}
	}

	if (showMask)
	{
		if (mxcsr & 0b10000000)
		{
			std::cout << "Invalid Operation Mask" << std::endl;
		}
		if (mxcsr & 0b100000000)
		{
			std::cout << "Denormal Mask" << std::endl;
		}
		if (mxcsr & 0b1000000000)
		{
			std::cout << "Divide By Zero Mask" << std::endl;
		}
		if (mxcsr & 0b10000000000)
		{
			std::cout << "Overflow Mask" << std::endl;
		}
		if (mxcsr & 0b100000000000)
		{
			std::cout << "Underflow Mask" << std::endl;
		}
		if (mxcsr & 0b1000000000000)
		{
			std::cout << "Precision Mask" << std::endl;
		}

		{
			const int round = (mxcsr & 0b110000000000000) >> 13;
			if (round == 0b00)
			{
				std::cout << "Round To Nearest" << std::endl;
			}
			else if (round == 0b01)
			{
				std::cout << "Round To Negative" << std::endl;
			}
			else if (round == 0b10)
			{
				std::cout << "Round To Positive" << std::endl;
			}
			else if (round == 0b11)
			{
				std::cout << "Round To Zero" << std::endl;
			}
		}

		if (mxcsr & 0b1000000000000000)
		{
			std::cout << "enable Flush To Zero" << std::endl;
		}
		if (mxcsr & 0b1000000)
		{
			std::cout << "enable Denormals Are Zero" << std::endl;
		}
	}

	return ret;
}