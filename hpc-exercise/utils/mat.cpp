#include "mat.h"
#include <x86intrin.h>
#include <cstring>
#include <iostream>


//Mat_8U
Mat_8U::Mat_8U(const unsigned char* data, const int rows, const int cols)
	: rows(rows), cols(cols), type(sizeof(unsigned char))
{
	const int size = type * rows * cols;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Mat_8U::Mat_8U(const int rows, const int cols)
	: rows(rows), cols(cols), type(sizeof(unsigned char))
{
	const int size = type * rows * cols;
	this->data = (unsigned char*)_mm_malloc(size, 32);
}

Mat_8U::Mat_8U(const Mat_8U& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char))
{
	const int size = type * rows * cols;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Mat_8U::Mat_8U(const Mat_16S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char))
{
	const int size = type * rows * cols;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols;j++)
	{
		data[j]=m.data[j];
	}
	memcpy((void*)this->data, (void*)m.data, size);
}

Mat_8U::Mat_8U(const Mat_32S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char))
{
	const int size = type * rows * cols;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols;j++)
	{
		data[j]=m.data[j];
	}
}

Mat_8U::Mat_8U(const Mat_32F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char))
{
	const int size = type * rows * cols;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols;j++)
	{
		data[j]=m.data[j];
	}
}

Mat_8U::Mat_8U(const Mat_64F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(unsigned char))
{
	const int size = type * rows * cols;
	this->data = (unsigned char*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols;j++)
	{
		data[j]=m.data[j];
	}
}

Mat_8U::~Mat_8U()
{
	_mm_free(data);
}

Mat_8U& Mat_8U::operator=(const Mat_8U& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(unsigned char);

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols;
		this->data = (unsigned char*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

void Mat_8U::show() const
{
	std::cout << "[ " << rows << " x " << cols << " ] 8U" << std::endl;

	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < cols; i++)
		{
			std::cout << (int)data[j*cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

//Mat_16S
Mat_16S::Mat_16S(const short* data, const int rows, const int cols)
	: rows(rows), cols(cols), type(sizeof(short))
{
	const int size = type * rows * cols;
	this->data = (short*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Mat_16S::Mat_16S(const int rows, const int cols)
	: rows(rows), cols(cols), type(sizeof(short))
{
	const int size = type * rows * cols;
	this->data = (short*)_mm_malloc(size, 32);
}

Mat_16S::Mat_16S(const Mat_16S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short))
{
	const int size = type * rows * cols;
	this->data = (short*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Mat_16S::Mat_16S(const Mat_8U& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short))
{
	const int size = type * rows * cols;
	this->data = (short*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols;j++)
	{
		data[j]=m.data[j];
	}
}

Mat_16S::Mat_16S(const Mat_32S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short))
{
	const int size = type * rows * cols;
	this->data = (short*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols;j++)
	{
		data[j]=m.data[j];
	}
}

Mat_16S::Mat_16S(const Mat_32F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short))
{
	const int size = type * rows * cols;
	this->data = (short*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols;j++)
	{
		data[j]=m.data[j];
	}
}

Mat_16S::Mat_16S(const Mat_64F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(short))
{
	const int size = type * rows * cols;
	this->data = (short*)_mm_malloc(size, 32);
	for(int j=0;j<rows*cols;j++)
	{
		data[j]=m.data[j];
	}
}

Mat_16S::~Mat_16S()
{
	_mm_free(data);
}

Mat_16S& Mat_16S::operator=(const Mat_16S& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(short);

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols;
		this->data = (short*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

void Mat_16S::show() const
{
	std::cout << "[ " << rows << " x " << cols << " ] 16S" << std::endl;

	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < cols; i++)
		{
			std::cout << (int)data[j*cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}


//}Mat_32S
Mat_32S::Mat_32S(const int* data, const int rows, const int cols)
	: rows(rows), cols(cols), type(sizeof(int))
{
	const int size = type * rows * cols;
	this->data = (int*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Mat_32S::Mat_32S(const int rows, const int cols)
	: rows(rows), cols(cols), type(sizeof(int))
{
	const int size = type * rows * cols;
	this->data = (int*)_mm_malloc(size, 32);
}

Mat_32S::Mat_32S(const Mat_32S& m)
	: rows(m.rows), cols(m.cols), type(sizeof(int))
{
	const int size = type * rows * cols;
	this->data = (int*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Mat_32S::Mat_32S(const Mat_8U& m)
        : rows(m.rows), cols(m.cols), type(sizeof(int))
{
        const int size = type * rows * cols;
        this->data = (int*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_32S::Mat_32S(const Mat_16S& m)
        : rows(m.rows), cols(m.cols), type(sizeof(int))
{
        const int size = type * rows * cols;
        this->data = (int*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_32S::Mat_32S(const Mat_32F& m)
        : rows(m.rows), cols(m.cols), type(sizeof(int))
{
        const int size = type * rows * cols;
        this->data = (int*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_32S::Mat_32S(const Mat_64F& m)
        : rows(m.rows), cols(m.cols), type(sizeof(int))
{
        const int size = type * rows * cols;
        this->data = (int*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_32S::~Mat_32S()
{
	_mm_free(data);
}

Mat_32S& Mat_32S::operator=(const Mat_32S& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(int);

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols;
		this->data = (int*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

void Mat_32S::show() const
{
	std::cout << "[ " << rows << " x " << cols << " ] 32S" << std::endl;

	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < cols; i++)
		{
			std::cout << (int)data[j*cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

//Mat_32F
Mat_32F::Mat_32F(const float* data, const int rows, const int cols)
	: type(sizeof(float)), rows(rows), cols(cols)
{
	const int size = type * rows * cols;
	this->data = (float*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Mat_32F::Mat_32F(const int rows, const int cols)
	: type(sizeof(float)), rows(rows), cols(cols)
{
	const int size = type * rows * cols; 
	this->data = (float*)_mm_malloc(size, 32);
}

Mat_32F::Mat_32F(const Mat_32F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(float))
{
	const int size = type * rows * cols;
	this->data = (float*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Mat_32F::Mat_32F(const Mat_8U& m)
        : rows(m.rows), cols(m.cols), type(sizeof(float))
{
        const int size = type * rows * cols;
        this->data = (float*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_32F::Mat_32F(const Mat_16S& m)
        : rows(m.rows), cols(m.cols), type(sizeof(float))
{
        const int size = type * rows * cols;
        this->data = (float*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_32F::Mat_32F(const Mat_32S& m)
        : rows(m.rows), cols(m.cols), type(sizeof(float))
{
        const int size = type * rows * cols;
        this->data = (float*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_32F::Mat_32F(const Mat_64F& m)
        : rows(m.rows), cols(m.cols), type(sizeof(float))
{
        const int size = type * rows * cols;
        this->data = (float*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_32F::~Mat_32F()
{
	_mm_free(data);
}

Mat_32F& Mat_32F::operator=(const Mat_32F& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(float);

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols; 
		this->data = (float*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

void Mat_32F::show() const
{
	std::cout << "[ " << rows << " x " << cols << " ] 32F" << std::endl;

	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < cols; i++)
		{
			std::cout << data[j*cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

//Mat_64F
Mat_64F::Mat_64F(const double* data, const int rows, const int cols)
	: type(sizeof(double)), rows(rows), cols(cols)
{
	const int size = type * rows * cols;
	this->data = (double*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)data, size);
}

Mat_64F::Mat_64F(const int rows, const int cols)
	: type(sizeof(double)), rows(rows), cols(cols)
{
	const int size = type * rows * cols; 
	this->data = (double*)_mm_malloc(size, 32);
}

Mat_64F::Mat_64F(const Mat_64F& m)
	: rows(m.rows), cols(m.cols), type(sizeof(double))
{
	const int size = type * rows * cols;
	this->data = (double*)_mm_malloc(size, 32);
	memcpy((void*)this->data, (void*)m.data, size);
}

Mat_64F::Mat_64F(const Mat_8U& m)
        : rows(m.rows), cols(m.cols), type(sizeof(double))
{
        const int size = type * rows * cols;
        this->data = (double*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_64F::Mat_64F(const Mat_16S& m)
        : rows(m.rows), cols(m.cols), type(sizeof(double))
{
        const int size = type * rows * cols;
        this->data = (double*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_64F::Mat_64F(const Mat_32S& m)
        : rows(m.rows), cols(m.cols), type(sizeof(double))
{
        const int size = type * rows * cols;
        this->data = (double*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_64F::Mat_64F(const Mat_32F& m)
        : rows(m.rows), cols(m.cols), type(sizeof(double))
{
        const int size = type * rows * cols;
        this->data = (double*)_mm_malloc(size, 32);
        for(int j=0;j<rows*cols;j++)
        {
                data[j]=m.data[j];
        }
}

Mat_64F::~Mat_64F()
{
	_mm_free(data);
}

Mat_64F& Mat_64F::operator=(const Mat_64F& m)
{
	if (this != &m)
	{
		this->rows = m.rows;
		this->cols = m.cols;
		this->type = sizeof(double);

		if (data != nullptr)
		{
			_mm_free(this->data);
		}
		const int size = type * rows * cols; 
		this->data = (double*)_mm_malloc(size, 32);
		memcpy((void*)this->data, (void*)m.data, size);
	}
	return *this;
}

void Mat_64F::show() const
{
	std::cout << "[ " << rows << " x " << cols << " ] 64F" << std::endl;

	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < cols; i++)
		{
			std::cout << data[j*cols + i] << ", ";
		}
		std::cout << std::endl;
	}
}

