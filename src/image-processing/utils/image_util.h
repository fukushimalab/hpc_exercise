#pragma once
//for using not accurate timer (1ms) if your system is not supported. Usually using more accurate timer(ns).
//#define USE_TIME_CHRONO

#include "image.h"
#include <vector>
#include <ctime>

#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif


//Image util functions
//read PXM
int readPXM(const char* name, Image_8U& dst);

//write PXM
void writePXM(const char* name, const Image_8U& src);

//color (RGB) to gray
void cvtColorGray(const Image_8U& src, Image_8U& dest);
void cvtColorGray(const Image_32F& src, Image_32F& dest);

//copy make border
void copyMakeBorder(const Image_8U& src, Image_8U& dest, const int top, const int bottom, const int left, const int right);
void copyMakeBorder(const Image_32F& src, Image_32F& dest, const int top, const int bottom, const int left, const int right);

//split image
void split(const Image_8U& src, Image_8U* dest);
void split(const Image_32F& src, Image_32F* dest);

//marge image
void merge(const Image_8U* src, const int channel, Image_8U& dest);
void merge(const Image_32F* src, const int channel, Image_32F& dest);

//init (zero)
void image_zero(Image_8U& m);
void image_zero(Image_16S& m);
void image_zero(Image_32S& m);
void image_zero(Image_32F& m);
void image_zero(Image_64F& m);

//init (one)
void image_one(Image_8U& m);
void image_one(Image_16S& m);
void image_one(Image_32S& m);
void image_one(Image_32F& m);
void image_one(Image_64F& m);

//init (rand)
void image_rand(Image_8U& m, const unsigned char rand_min, const unsigned char rand_max);
void image_rand(Image_16S& m, const short rand_min, const short rand_max);
void image_rand(Image_32S& m, const int rand_min, const int rand_max);
void image_rand(Image_32F& m, const float rand_min, const float rand_max);
void image_rand(Image_64F& m, const double rand_min, const double rand_max);

//info
void image_info(const Image_8U m);
void image_info(const Image_16S m);
void image_info(const Image_32S m);
void image_info(const Image_32F m);
void image_info(const Image_64F m);

//rand
unsigned char rand_8u(const unsigned char rand_min, const unsigned char rand_max);
short rand_16s(const short rand_min, const short rand_max);
int rand_32s(const int rand_min, const int rand_max);
float rand_32f(const float rand_min, const float rand_max);
double rand_64f(const double rand_min, const double rand_max);

//calc PSNR
double calcPSNR(const Image_8U& src1, const Image_8U& src2);
double calcPSNR(const Image_32F& src1, const Image_32F& src2);

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
