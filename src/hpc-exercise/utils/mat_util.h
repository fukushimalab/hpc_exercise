#pragma once
#include "mat.h"
#include <vector>
#include <time.h>

//Mat util functions
//init (zero)
void mat_zero(Mat_8U& m);
void mat_zero(Mat_16S& m);
void mat_zero(Mat_32S& m);
void mat_zero(Mat_32F& m);
void mat_zero(Mat_64F& m);

//init (one)
void mat_one(Mat_8U& m);
void mat_one(Mat_16S& m);
void mat_one(Mat_32S& m);
void mat_one(Mat_32F& m);
void mat_one(Mat_64F& m);

//init (rand)
void mat_rand(Mat_8U& m, const unsigned char rand_min, const unsigned char rand_max);
void mat_rand(Mat_16S& m, const short rand_min, const short rand_max);
void mat_rand(Mat_32S& m, const int rand_min, const int rand_max);
void mat_rand(Mat_32F& m, const float rand_min, const float rand_max);
void mat_rand(Mat_64F& m, const double rand_min, const double rand_max);

//val add
Mat_8U mat_add(const Mat_8U m1, const unsigned char v);
Mat_16S mat_add(const Mat_16S m1, const short v);
Mat_32S mat_add(const Mat_32S m1, const int v);
Mat_32F mat_add(const Mat_32F m1, const float v);
Mat_64F mat_add(const Mat_64F m1, const double v);

//mat add
Mat_8U mat_add(const Mat_8U m1, const Mat_8U m2);
Mat_16S mat_add(const Mat_16S m1, const Mat_16S m2);
Mat_32S mat_add(const Mat_32S m1, const Mat_32S m2);
Mat_32F mat_add(const Mat_32F m1, const Mat_32F m2);
Mat_64F mat_add(const Mat_64F m1, const Mat_64F m2);

//val mul
Mat_8U mat_mul(const Mat_8U m1, const unsigned char v);
Mat_16S mat_mul(const Mat_16S m1, const short v);
Mat_32S mat_mul(const Mat_32S m1, const int v);
Mat_32F mat_mul(const Mat_32F m1, const float v);
Mat_64F mat_mul(const Mat_64F m1, const double v);

//mat mul
Mat_8U  mat_mul(const Mat_8U m1, const Mat_8U m2);
Mat_16S mat_mul(const Mat_16S m1, const Mat_16S m2);
Mat_32S mat_mul(const Mat_32S m1, const Mat_32S m2);
Mat_32F mat_mul(const Mat_32F m1, const Mat_32F m2);
Mat_64F mat_mul(const Mat_64F m1, const Mat_64F m2);

//val div
Mat_8U mat_div(const Mat_8U m1, const unsigned char v);
Mat_16S mat_div(const Mat_16S m1, const short v);
Mat_32S mat_div(const Mat_32S m1, const int v);
Mat_32F mat_div(const Mat_32F m1, const float v);
Mat_64F mat_div(const Mat_64F m1, const double v);

//show
void mat_show(const Mat_8U m);
void mat_show(const Mat_16S m);
void mat_show(const Mat_32S m);
void mat_show(const Mat_32F m);
void mat_show(const Mat_64F m);

//rand
unsigned char rand_8u(const unsigned char rand_min, const unsigned char rand_max);
short rand_16s(const short rand_min, const short rand_max);
int rand_32s(const int rand_min, const int rand_max);
float rand_32f(const float rand_min, const float rand_max);
double rand_64f(const double rand_min, const double rand_max);

//timer
struct CalcTime
{
	std::vector<double> que;
	timespec s;
	timespec e;
	void start();
	void end();
	void clear();

	double getAvgTime(const bool dropFirstMeasure = true, const bool isClear = true);
	double getLastTime();
};
