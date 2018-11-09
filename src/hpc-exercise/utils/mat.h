#pragma once


struct Mat_8U;
struct Mat_16S;
struct Mat_32S;
struct Mat_32F;
struct Mat_64F;

struct Mat_8U
{
	unsigned char* data = nullptr;
	int type;
	int rows;
	int cols;

	Mat_8U(const unsigned char* data, const int rows, const int cols);
	Mat_8U(const int rows, const int cols);
	Mat_8U(const Mat_8U& m);
	Mat_8U(const Mat_16S& m);
	Mat_8U(const Mat_32S& m);
	Mat_8U(const Mat_32F& m);
	Mat_8U(const Mat_64F& m);
	~Mat_8U();
	Mat_8U& operator=(const Mat_8U& m);
	void show() const;
	int index(const int row, const int col) const
	{
		return row*cols + col;
	}
};

struct Mat_16S
{
	short* data = nullptr;
	int type;
	int rows;
	int cols;

	Mat_16S(const short* data, const int rows, const int cols);
	Mat_16S(const int rows, const int cols);
	Mat_16S(const Mat_16S& m);
	Mat_16S(const Mat_8U& m);
	Mat_16S(const Mat_32S& m);
	Mat_16S(const Mat_32F& m);
	Mat_16S(const Mat_64F& m);
	~Mat_16S();
	Mat_16S& operator=(const Mat_16S& m);
	void show() const;
	int index(const int row, const int col) const
	{
		return row*cols + col;
	}
};

struct Mat_32S
{
	int* data = nullptr;
	int type;
	int rows;
	int cols;

	Mat_32S(const int* data, const int rows, const int cols);
	Mat_32S(const int rows, const int cols);
	Mat_32S(const Mat_32S& m);
	Mat_32S(const Mat_8U& m);
	Mat_32S(const Mat_16S& m);
	Mat_32S(const Mat_32F& m);
	Mat_32S(const Mat_64F& m);
	~Mat_32S();
	Mat_32S& operator=(const Mat_32S& m);
	void show() const;
	int index(const int row, const int col) const
	{
		return row*cols + col;
	}
};

struct Mat_32F
{
	float* data = nullptr;
	int type;
	int rows;
	int cols;

	Mat_32F(const float* data, const int rows, const int cols);
	Mat_32F(const int rows, const int cols);
	Mat_32F(const Mat_32F& m);
	Mat_32F(const Mat_8U& m);
	Mat_32F(const Mat_16S& m);
	Mat_32F(const Mat_32S& m);
	Mat_32F(const Mat_64F& m);
	~Mat_32F();
	Mat_32F& operator=(const Mat_32F& m);
	void show() const;
	int index(const int row, const int col) const
	{
		return row*cols + col;
	}
};

struct Mat_64F
{
	double* data = nullptr;
	int type;
	int rows;
	int cols;

	Mat_64F(const double* data, const int rows, const int cols);
	Mat_64F(const int rows, const int cols);
	Mat_64F(const Mat_64F& m);
	Mat_64F(const Mat_8U& m);
	Mat_64F(const Mat_16S& m);
	Mat_64F(const Mat_32S& m);
	Mat_64F(const Mat_32F& m);
	~Mat_64F();
	Mat_64F& operator=(const Mat_64F& m);
	void show() const;
	int index(const int row, const int col) const
	{
		return row*cols + col;
	}
};

