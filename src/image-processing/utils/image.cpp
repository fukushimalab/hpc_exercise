#include "image.h"
#ifdef __GNUC__
#include <x86intrin.h>
#elif _MSC_VER
#include <intrin.h>
#endif
#include <cstring>
#include <iostream>


//Image_8U
Image_8U::Image_8U(const unsigned char* data, const int rows, const int cols, const int channels)
	: rows(rows), cols(cols), type(sizeof(unsigned char)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Image_8U::Image_8U(const int rows, const int cols, const int channels)
	: rows(rows), cols(cols), type(sizeof(unsigned char)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (unsigned char*)_mm_malloc(size, 32);
}

Image_8U::Image_8U(const Image_8U& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char)), channels(m.channels)
{
	const int size = type * rows * cols * channels;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Image_8U::Image_8U(const Image_16S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char)), channels(m.channels)
{
	const int size = type * rows * cols * channels;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols*channels*channels;j++)
	{
		data[j]=(unsigned char)std::min(std::max(m.data[j], (short)0), (short)255);
	}
	memcpy((void*)this->data, (void*)m.data, size);
}

Image_8U::Image_8U(const Image_32S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char)), channels(m.channels)
{
	const int size = type * rows * cols * channels;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols*channels;j++)
	{
		data[j]=(unsigned char)std::min(std::max(m.data[j], (int)0), (int)255);
	}
}

Image_8U::Image_8U(const Image_32F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char)), channels(m.channels) 
{
	const int size = type * rows * cols * channels;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols*channels;j++)
	{
		data[j]=(unsigned char)std::min(std::max(m.data[j], (float)0), (float)255);
	}
}

Image_8U::Image_8U(const Image_64F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char)), channels(m.channels) 
{
	const int size = type * rows * cols * channels;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols*channels;j++)
	{
		data[j]=(unsigned char)std::min(std::max(m.data[j], (double)0), (double)255);
	}
}

Image_8U::~Image_8U()
{
	_mm_free(data);
}

Image_8U& Image_8U::operator=(const Image_8U& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(unsigned char);
		this->channels = m.channels;

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols * channels;
		this->data = (unsigned char*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

//Image_16S
Image_16S::Image_16S(const short* data, const int rows, const int cols, const int channels)
	: rows(rows), cols(cols), type(sizeof(short)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (short*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Image_16S::Image_16S(const int rows, const int cols, const int channels)
	: rows(rows), cols(cols), type(sizeof(short)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (short*)_mm_malloc(size, 32);
}

Image_16S::Image_16S(const Image_16S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (short*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Image_16S::Image_16S(const Image_8U& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (short*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols*channels;j++)
	{
		data[j]=m.data[j];
	}
}

Image_16S::Image_16S(const Image_32S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (short*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols*channels;j++)
	{
		data[j]=m.data[j];
	}
}

Image_16S::Image_16S(const Image_32F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (short*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols*channels;j++)
	{
		data[j]=m.data[j];
	}
}

Image_16S::Image_16S(const Image_64F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (short*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols*channels*channels;j++)
	{
		data[j]=m.data[j];
	}
}

Image_16S::~Image_16S()
{
	_mm_free(data);
}

Image_16S& Image_16S::operator=(const Image_16S& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(short);
		this->channels = m.channels;

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols * channels;
		this->data = (short*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

//Image_32S
Image_32S::Image_32S(const int* data, const int rows, const int cols, const int channels)
	: rows(rows), cols(cols), type(sizeof(int)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (int*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Image_32S::Image_32S(const int rows, const int cols, const int channels)
	: rows(rows), cols(cols), type(sizeof(int)), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (int*)_mm_malloc(size, 32);
}

Image_32S::Image_32S(const Image_32S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(int)), channels(m.channels)
{
	const int size = type * rows * cols * channels;
	this->data = (int*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Image_32S::Image_32S(const Image_8U& m)
	: rows(m.rows), cols(m.cols), type(sizeof(int)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (int*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_32S::Image_32S(const Image_16S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(int)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (int*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_32S::Image_32S(const Image_32F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(int)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (int*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_32S::Image_32S(const Image_64F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(int)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (int*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_32S::~Image_32S()
{
	_mm_free(data);
}

Image_32S& Image_32S::operator=(const Image_32S& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(int);
		this->channels = m.channels;

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols * channels;
		this->data = (int*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

//Image_32F
Image_32F::Image_32F(const float* data, const int rows, const int cols, const int channels)
	: type(sizeof(float)), rows(rows), cols(cols), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (float*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Image_32F::Image_32F(const int rows, const int cols, const int channels)
	: type(sizeof(float)), rows(rows), cols(cols), channels(channels)
{
	const int size = type * rows * cols * channels; 
	this->data = (float*)_mm_malloc(size, 32);
}

Image_32F::Image_32F(const Image_32F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(float)), channels(m.channels)
{
	const int size = type * rows * cols * channels;
	this->data = (float*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Image_32F::Image_32F(const Image_8U& m)
	: rows(m.rows), cols(m.cols), type(sizeof(float)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (float*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_32F::Image_32F(const Image_16S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(float)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (float*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_32F::Image_32F(const Image_32S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(float)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (float*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_32F::Image_32F(const Image_64F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(float)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (float*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_32F::~Image_32F()
{
	_mm_free(data);
}

Image_32F& Image_32F::operator=(const Image_32F& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(float);
		this->channels = m.channels;

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols * channels; 
		this->data = (float*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

//Image_64F
Image_64F::Image_64F(const double* data, const int rows, const int cols, const int channels)
	: type(sizeof(double)), rows(rows), cols(cols), channels(channels)
{
	const int size = type * rows * cols * channels;
	this->data = (double*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Image_64F::Image_64F(const int rows, const int cols, const int channels)
	: type(sizeof(double)), rows(rows), cols(cols), channels(channels)
{
	const int size = type * rows * cols * channels; 
	this->data = (double*)_mm_malloc(size, 32);
}

Image_64F::Image_64F(const Image_64F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(double)), channels(m.channels)
{
	const int size = type * rows * cols * channels;
	this->data = (double*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Image_64F::Image_64F(const Image_8U& m)
	: rows(m.rows), cols(m.cols), type(sizeof(double)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (double*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_64F::Image_64F(const Image_16S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(double)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (double*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_64F::Image_64F(const Image_32S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(double)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (double*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_64F::Image_64F(const Image_32F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(double)), channels(m.channels)
{
        const int size = type * rows * cols * channels;
        this->data = (double*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols*channels;j++)
        {
                data[j]=m.data[j];
        }
}

Image_64F::~Image_64F()
{
	_mm_free(data);
}

Image_64F& Image_64F::operator=(const Image_64F& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(double);
		this->channels = m.channels;

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols * channels; 
		this->data = (double*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

