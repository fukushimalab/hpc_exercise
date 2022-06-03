#pragma once
//for using not accurate timer (1ms) if your system is not supported. Usually using more accurate timer(ns).
//#define USE_TIME_CHRONO

#include "mat.h"
#include <vector>
#include <ctime>

#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif


//Mat util functions
//init (zero)
void mat_zero(Mat_8S& m);
void mat_zero(Mat_8U& m);
void mat_zero(Mat_16S& m);
void mat_zero(Mat_32S& m);
void mat_zero(Mat_32F& m);
void mat_zero(Mat_64F& m);

//init (one)
void mat_one(Mat_8S& m);
void mat_one(Mat_8U& m);
void mat_one(Mat_16S& m);
void mat_one(Mat_32S& m);
void mat_one(Mat_32F& m);
void mat_one(Mat_64F& m);

//init (rand)
void mat_rand(Mat_8S& m, const char rand_min, const char rand_max);
void mat_rand(Mat_8U& m, const unsigned char rand_min, const unsigned char rand_max);
void mat_rand(Mat_16S& m, const short rand_min, const short rand_max);
void mat_rand(Mat_32S& m, const int rand_min, const int rand_max);
void mat_rand(Mat_32F& m, const float rand_min, const float rand_max);
void mat_rand(Mat_64F& m, const double rand_min, const double rand_max);

//val add
Mat_8S mat_add(const Mat_8S& m1, const char v);
Mat_8U mat_add(const Mat_8U& m1, const unsigned char v);
Mat_16S mat_add(const Mat_16S& m1, const short v);
Mat_32S mat_add(const Mat_32S& m1, const int v);
Mat_32F mat_add(const Mat_32F& m1, const float v);
Mat_64F mat_add(const Mat_64F& m1, const double v);

//mat add
Mat_8S mat_add(const Mat_8S& m1, const Mat_8S& m2);
Mat_8U mat_add(const Mat_8U& m1, const Mat_8U& m2);
Mat_16S mat_add(const Mat_16S& m1, const Mat_16S& m2);
Mat_32S mat_add(const Mat_32S& m1, const Mat_32S& m2);
Mat_32F mat_add(const Mat_32F& m1, const Mat_32F& m2);
Mat_64F mat_add(const Mat_64F& m1, const Mat_64F& m2);

void mat_add(const Mat_8S& m1, const Mat_8S& m2, Mat_8S& dest);
void mat_add(const Mat_8U& m1, const Mat_8U& m2, Mat_8U& dest);
void mat_add(const Mat_16S& m1, const Mat_16S& m2, Mat_16S& dest);
void mat_add(const Mat_32S& m1, const Mat_32S& m2, Mat_32S& dest);
void mat_add(const Mat_32F& m1, const Mat_32F& m2, Mat_32F& dest);
void mat_add(const Mat_64F& m1, const Mat_64F& m2, Mat_64F& dest);

void mat_add_scalar(const Mat_8S& m1, const Mat_8S& m2, Mat_8S& dest);//for auto vectorization test
void mat_add_scalar(const Mat_8U& m1, const Mat_8U& m2, Mat_8U& dest);//for auto vectorization test

//val mul
Mat_8S mat_mul(const Mat_8S& m1, const char v);
Mat_8U mat_mul(const Mat_8U& m1, const unsigned char v);
Mat_16S mat_mul(const Mat_16S& m1, const short v);
Mat_32S mat_mul(const Mat_32S& m1, const int v);
Mat_32F mat_mul(const Mat_32F& m1, const float v);
Mat_64F mat_mul(const Mat_64F& m1, const double v);

//mat mul
Mat_8S  mat_mul(const Mat_8S& m1, const Mat_8S& m2);
Mat_8U  mat_mul(const Mat_8U& m1, const Mat_8U& m2);
Mat_16S mat_mul(const Mat_16S& m1, const Mat_16S& m2);
Mat_32S mat_mul(const Mat_32S& m1, const Mat_32S& m2);
Mat_32F mat_mul(const Mat_32F& m1, const Mat_32F& m2);
Mat_64F mat_mul(const Mat_64F& m1, const Mat_64F& m2);

//val div
Mat_8S mat_div(const Mat_8S& m1, const char v);
Mat_8U mat_div(const Mat_8U& m1, const unsigned char v);
Mat_16S mat_div(const Mat_16S& m1, const short v);
Mat_32S mat_div(const Mat_32S& m1, const int v);
Mat_32F mat_div(const Mat_32F& m1, const float v);
Mat_64F mat_div(const Mat_64F& m1, const double v);

//show
void mat_show(const Mat_8S& m);
void mat_show(const Mat_8U& m);
void mat_show(const Mat_16S& m);
void mat_show(const Mat_32S& m);
void mat_show(const Mat_32F& m);
void mat_show(const Mat_64F& m);

//rand
unsigned char rand_8u(const unsigned char rand_min, const unsigned char rand_max);
short rand_16s(const short rand_min, const short rand_max);
int rand_32s(const int rand_min, const int rand_max);
float rand_32f(const float rand_min, const float rand_max);
double rand_64f(const double rand_min, const double rand_max);

//diff
int mat_diff(Mat_8U& src1, Mat_8U& src2);
int mat_diff(Mat_16S& src1, Mat_16S& src2);
int mat_diff(Mat_32S& src1, Mat_32S& src2);
double mat_diff(Mat_32F& src1, Mat_32F& src2);
double mat_diff(Mat_64F& src1, Mat_64F& src2);

//timer
#ifdef USE_TIME_CHRONO
#include <chrono>
struct CalcTime
{
	std::vector<double> que;
	std::chrono::system_clock::time_point s;
	std::chrono::system_clock::time_point e;
	void start();
	void end();
	void clear();

	double getAvgTime(const bool dropFirstMeasure = true, const bool isClear = true);
	double getLastTime();
};
#elif defined(__GNUC__)
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
#elif defined(_MSC_VER)
struct CalcTime
{
	std::vector<double> que;
	LARGE_INTEGER s;
	LARGE_INTEGER e;

	LARGE_INTEGER frequency;
	void start();
	void end();
	void clear();

	double getAvgTime(const bool dropFirstMeasure = true, const bool isClear = true);
	double getLastTime();
	CalcTime();
};
#endif

bool show_mxcsr(const bool showState = true, const bool showMask = false, const bool isClaer = false);