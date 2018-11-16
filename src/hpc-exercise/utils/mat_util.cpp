#include <stdlib.h>
#include "mat_util.h"
#include <iostream>

//Mat util functions
////////////////////////////
//init (zero)
////////////////////////////
void mat_zero(Mat_8U& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 0;
		}
	}
}

void mat_zero(Mat_16S& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 0;
		}
	}
}

void mat_zero(Mat_32S& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 0;
		}
	}
}

void mat_zero(Mat_32F& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 0;
		}
	}
}

void mat_zero(Mat_64F& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 0;
		}
	}
}


////////////////////////////
//init (one)
////////////////////////////
void mat_one(Mat_8U& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 1;
		}
	}
}

void mat_one(Mat_16S& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 1;
		}
	}
}

void mat_one(Mat_32S& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 1;
		}
	}
}

void mat_one(Mat_32F& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 1.f;
		}
	}
}

void mat_one(Mat_64F& m)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = 1.0;
		}
	}
}


////////////////////////////
//init (rand)
////////////////////////////
void mat_rand(Mat_8U& m, const unsigned char rand_min, const unsigned char rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = rand_min + (rand() * (rand_max - rand_min + 1.0f) / (1.0f + RAND_MAX));
		}
	}
}

void mat_rand(Mat_16S& m, const short rand_min, const short rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = rand_min + (rand() * (rand_max - rand_min + 1.0f) / (1.0f + RAND_MAX));
		}
	}
}

void mat_rand(Mat_32S& m, const int rand_min, const int rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = rand_min + (rand() * (rand_max - rand_min + 1.0f) / (1.0f + RAND_MAX));
		}
	}
}

void mat_rand(Mat_32F& m, const float rand_min, const float rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = rand_min + (rand() * (rand_max - rand_min + 1.0f) / (1.0f + RAND_MAX));
		}
	}
}

void mat_rand(Mat_64F& m, const double rand_min, const double rand_max)
{
	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			m.data[j*m.cols + i] = rand_min + (rand() * (rand_max - rand_min + 1.0) / (1.0 + RAND_MAX));
		}
	}
}

////////////////////////////
//val add 
////////////////////////////
Mat_8U mat_add(const Mat_8U m, const unsigned char v)
{
	Mat_8U ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] + v;
		}
	}
	return ret;
}

Mat_16S mat_add(const Mat_16S m, const short v)
{
	Mat_16S ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] + v;
		}
	}
	return ret;
}

Mat_32S mat_add(const Mat_32S m, const int v)
{
	Mat_32S ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] + v;
		}
	}
	return ret;
}

Mat_32F mat_add(const Mat_32F m, const float v)
{
	Mat_32F ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] + v;
		}
	}
	return ret;
}

Mat_64F mat_add(const Mat_64F m, const double v)
{
	Mat_64F ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] + v;
		}
	}
	return ret;
}


////////////////////////////
//mat add
////////////////////////////
Mat_8U mat_add(const Mat_8U m1, const Mat_8U m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_8U ret(m1.rows, m1.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = m1.data[j*m1.cols + i] + m2.data[j*m2.cols + i];
		}
	}
	return  ret;
}

Mat_16S mat_add(const Mat_16S m1, const Mat_16S m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_16S ret(m1.rows, m1.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = m1.data[j*m1.cols + i] + m2.data[j*m2.cols + i];
		}
	}
	return  ret;
}

Mat_32S mat_add(const Mat_32S m1, const Mat_32S m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_32S ret(m1.rows, m1.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = m1.data[j*m1.cols + i] + m2.data[j*m2.cols + i];
		}
	}
	return  ret;
}

Mat_32F mat_add(const Mat_32F m1, const Mat_32F m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_32F ret(m1.rows, m1.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = m1.data[j*m1.cols + i] + m2.data[j*m2.cols + i];
		}
	}
	return  ret;
}

Mat_64F mat_add(const Mat_64F m1, const Mat_64F m2)
{
	if (m1.rows != m2.rows || m1.cols != m2.cols)
	{
		std::cout << "invalid mat size (mat add)" << std::endl;
		exit(-1);
	}

	Mat_64F ret(m1.rows, m1.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = m1.data[j*m1.cols + i] + m2.data[j*m2.cols + i];
		}
	}
	return  ret;
}


////////////////////////////
//val mul
////////////////////////////
Mat_8U mat_mul(const Mat_8U m, const unsigned char v)
{
	Mat_8U ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] * v;
		}
	}
	return ret;
}

Mat_16S mat_mul(const Mat_16S m, const short v)
{
	Mat_16S ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] * v;
		}
	}
	return ret;
}

Mat_32S mat_mul(const Mat_32S m, const int v)
{
	Mat_32S ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] * v;
		}
	}
	return ret;
}

Mat_32F mat_mul(const Mat_32F m, const float v)
{
	Mat_32F ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] * v;
		}
	}
	return ret;
}

Mat_64F mat_mul(const Mat_64F m, const double v)
{
	Mat_64F ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] * v;
		}
	}
	return ret;
}


////////////////////////////
//mat mul
////////////////////////////
Mat_8U mat_mul(const Mat_8U m1, const Mat_8U m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_8U ret(m1.rows, m2.cols);
	mat_zero(ret);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			for (int k = 0; k < m1.cols; k++)
			{
				ret.data[j*ret.cols + i] += m1.data[j*m1.cols + k] * m2.data[k*m2.cols + i];
			}
		}
	}
	return ret;
}

Mat_16S mat_mul(const Mat_16S m1, const Mat_16S m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_16S ret(m1.rows, m2.cols);
	mat_zero(ret);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			for (int k = 0; k < m1.cols; k++)
			{
				ret.data[j*ret.cols + i] += m1.data[j*m1.cols + k] * m2.data[k*m2.cols + i];
			}
		}
	}
	return ret;
}

Mat_32S mat_mul(const Mat_32S m1, const Mat_32S m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_32S ret(m1.rows, m2.cols);
	mat_zero(ret);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			for (int k = 0; k < m1.cols; k++)
			{
				ret.data[j*ret.cols + i] += m1.data[j*m1.cols + k] * m2.data[k*m2.cols + i];
			}
		}
	}
	return ret;
}

Mat_32F mat_mul(const Mat_32F m1, const Mat_32F m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_32F ret(m1.rows, m2.cols);
	mat_zero(ret);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			for (int k = 0; k < m1.cols; k++)
			{
				ret.data[j*ret.cols + i] += m1.data[j*m1.cols + k] * m2.data[k*m2.cols + i];
			}
		}
	}
	return ret;
}

Mat_64F mat_mul(const Mat_64F m1, const Mat_64F m2)
{
	if (m1.rows != m2.cols)
	{
		std::cout << "invalid mat size (mat mul)" << std::endl;
		exit(-1);
	}

	Mat_64F ret(m1.rows, m2.cols);
	mat_zero(ret);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			for (int k = 0; k < m1.cols; k++)
			{
				ret.data[j*ret.cols + i] += m1.data[j*m1.cols + k] * m2.data[k*m2.cols + i];
			}
		}
	}
	return ret;
}


////////////////////////////
//val div
////////////////////////////
Mat_8U mat_div(const Mat_8U m, const unsigned char v)
{
	Mat_8U ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] / v;
		}
	}
	return ret;
}

Mat_16S mat_div(const Mat_16S m, const short v)
{
	Mat_16S ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] / v;
		}
	}
	return ret;
}

Mat_32S mat_div(const Mat_32S m, const int v)
{
	Mat_32S ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] / v;
		}
	}
	return ret;
}

Mat_32F mat_div(const Mat_32F m, const float v)
{
	Mat_32F ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] / v;
		}
	}
	return ret;
}

Mat_64F mat_div(const Mat_64F m, const double v)
{
	Mat_64F ret(m.rows, m.cols);
	for (int j = 0; j < ret.rows; j++)
	{
		for (int i = 0; i < ret.cols; i++)
		{
			ret.data[j*ret.cols + i] = ret.data[j*ret.cols + i] / v;
		}
	}
	return ret;
}


////////////////////////////
//show
////////////////////////////
void mat_show(const Mat_8U m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 8U" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << (int)m.data[j*m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

void mat_show(const Mat_16S m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 16S" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << m.data[j*m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

void mat_show(const Mat_32S m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 32S" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << m.data[j*m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

void mat_show(const Mat_32F m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 32F" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << m.data[j*m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

void mat_show(const Mat_64F m)
{
	std::cout << "[ " << m.rows << " x " << m.cols << " ] 64F" << std::endl;

	for (int j = 0; j < m.rows; j++)
	{
		for (int i = 0; i < m.cols; i++)
		{
			std::cout << m.data[j*m.cols + i] << ", ";
		}
		std::cout << std::endl;
	}
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

