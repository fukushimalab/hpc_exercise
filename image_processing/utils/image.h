#pragma once


struct Image_8U;
struct Image_16S;
struct Image_32S;
struct Image_32F;
struct Image_64F;

struct Image_8U
{
	unsigned char* data = nullptr;
	int type;
	int rows;
	int cols;
	int channels;

	Image_8U(){};
	Image_8U(const unsigned char* data, const int rows, const int cols, const int channels);
	Image_8U(const int rows, const int cols, const int channels);
	Image_8U(const Image_8U& m);
	Image_8U(const Image_16S& m);
	Image_8U(const Image_32S& m);
	Image_8U(const Image_32F& m);
	Image_8U(const Image_64F& m);
	~Image_8U();
	Image_8U& operator=(const Image_8U& m);
	void show() const;
};

struct Image_16S
{
	short* data = nullptr;
	int type;
	int rows;
	int cols;
	int channels;

	Image_16S(){};
	Image_16S(const short* data, const int rows, const int cols, const int channels);
	Image_16S(const int rows, const int cols, const int channels);
	Image_16S(const Image_16S& m);
	Image_16S(const Image_8U& m);
	Image_16S(const Image_32S& m);
	Image_16S(const Image_32F& m);
	Image_16S(const Image_64F& m);
	~Image_16S();
	Image_16S& operator=(const Image_16S& m);
	void show() const;
};

struct Image_32S
{
	int* data = nullptr;
	int type;
	int rows;
	int cols;
	int channels;

	Image_32S(){};
	Image_32S(const int* data, const int rows, const int cols, const int channels);
	Image_32S(const int rows, const int cols, const int channels);
	Image_32S(const Image_32S& m);
	Image_32S(const Image_8U& m);
	Image_32S(const Image_16S& m);
	Image_32S(const Image_32F& m);
	Image_32S(const Image_64F& m);
	~Image_32S();
	Image_32S& operator=(const Image_32S& m);
	void show() const;
};

struct Image_32F
{
	float* data = nullptr;
	int type;
	int rows;
	int cols;
	int channels;

	Image_32F(){};
	Image_32F(const float* data, const int rows, const int cols, const int channels);
	Image_32F(const int rows, const int cols, const int channels);
	Image_32F(const Image_32F& m);
	Image_32F(const Image_8U& m);
	Image_32F(const Image_16S& m);
	Image_32F(const Image_32S& m);
	Image_32F(const Image_64F& m);
	~Image_32F();
	Image_32F& operator=(const Image_32F& m);
	void show() const;
};

struct Image_64F
{
	double* data = nullptr;
	int type;
	int rows;
	int cols;
	int channels;

	Image_64F(){};
	Image_64F(const double* data, const int rows, const int cols, const int channels);
	Image_64F(const int rows, const int cols, const int channels);
	Image_64F(const Image_64F& m);
	Image_64F(const Image_8U& m);
	Image_64F(const Image_16S& m);
	Image_64F(const Image_32S& m);
	Image_64F(const Image_32F& m);
	~Image_64F();
	Image_64F& operator=(const Image_64F& m);
	void show() const;
};

