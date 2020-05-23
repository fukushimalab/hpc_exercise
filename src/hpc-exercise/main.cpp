#include "utils/mat.h"
#include "utils/mat_util.h"
#include "utils/simd_util.h"
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>

void inline _mm256_transpose_8x8_ps(__m256* dst, const __m256* src);
void inline rot(double a, double b, double& x, double& y, double radian);

void loofline_test(const int size, const int iteration)
{
	//FLOPS計算
	//x+1: (1*x+1)
	//x*x+x+1: x*(x+1)+1
	//x*x*x+x*x+x+1: x*(x*(x+1))+1

	std::cout << "loofline test" << std::endl;
	const int loop = iteration;

	CalcTime t;

	Mat_32F x(1, size);
	Mat_32F y(1, size);
	mat_rand(x, 1.f, 1.2f);
	mat_zero(y);

	int simdsize8 = size / 8;
	int simdsize16 = size / 16;

	printf("size %d, iteration\n", iteration);
	printf("order, GFLOPS, FLOPS/BYTE\n");

	int n = 0;

	//---------------------------------------------
	n = 1;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 2;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	
	//---------------------------------------------
	n = 3;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			
			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 4;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 5;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 6;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 7;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 8;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 9;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 10;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 11;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 12;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 14;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 14;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 15;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 16;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 17;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16
			ret = _mm256_fmadd_ps(mx, ret, mones);//17

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
			ret = _mm256_fmadd_ps(mx, ret, mones);//17
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 18;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16
			ret = _mm256_fmadd_ps(mx, ret, mones);//17
			ret = _mm256_fmadd_ps(mx, ret, mones);//18

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
			ret = _mm256_fmadd_ps(mx, ret, mones);//17
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
			ret = _mm256_fmadd_ps(mx, ret, mones);//18
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 19;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16
			ret = _mm256_fmadd_ps(mx, ret, mones);//17
			ret = _mm256_fmadd_ps(mx, ret, mones);//18
			ret = _mm256_fmadd_ps(mx, ret, mones);//19

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
			ret = _mm256_fmadd_ps(mx, ret, mones);//17
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
			ret = _mm256_fmadd_ps(mx, ret, mones);//18
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18
			ret = _mm256_fmadd_ps(mx, ret, mones);//19
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//19
			
			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));

	//---------------------------------------------
	n = 20;
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize8; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16
			ret = _mm256_fmadd_ps(mx, ret, mones);//17
			ret = _mm256_fmadd_ps(mx, ret, mones);//18
			ret = _mm256_fmadd_ps(mx, ret, mones);//19
			ret = _mm256_fmadd_ps(mx, ret, mones);//20

			_mm256_store_ps(py, ret);

			px += 8;
			py += 8;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	for (int j = 0; j < loop; j++)
	{
		t.start();
		float* px = x.data;
		float* py = y.data;
		const __m256 mones = _mm256_set1_ps(1.f);
		for (int i = simdsize16; i != 0; i--)
		{
			const __m256 mx = _mm256_load_ps(px);
			const __m256 mx2 = _mm256_load_ps(px + 8);

			__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
			__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
			ret = _mm256_fmadd_ps(mx, ret, mones);//2
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
			ret = _mm256_fmadd_ps(mx, ret, mones);//3
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
			ret = _mm256_fmadd_ps(mx, ret, mones);//4
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
			ret = _mm256_fmadd_ps(mx, ret, mones);//5
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
			ret = _mm256_fmadd_ps(mx, ret, mones);//6
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
			ret = _mm256_fmadd_ps(mx, ret, mones);//7
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
			ret = _mm256_fmadd_ps(mx, ret, mones);//8
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
			ret = _mm256_fmadd_ps(mx, ret, mones);//9
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
			ret = _mm256_fmadd_ps(mx, ret, mones);//10
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10
			ret = _mm256_fmadd_ps(mx, ret, mones);//11
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//11
			ret = _mm256_fmadd_ps(mx, ret, mones);//12
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//12
			ret = _mm256_fmadd_ps(mx, ret, mones);//13
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//13
			ret = _mm256_fmadd_ps(mx, ret, mones);//14
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//14
			ret = _mm256_fmadd_ps(mx, ret, mones);//15
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//15
			ret = _mm256_fmadd_ps(mx, ret, mones);//16
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//16
			ret = _mm256_fmadd_ps(mx, ret, mones);//17
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//17
			ret = _mm256_fmadd_ps(mx, ret, mones);//18
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//18
			ret = _mm256_fmadd_ps(mx, ret, mones);//19
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//19
			ret = _mm256_fmadd_ps(mx, ret, mones);//20
			ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//20

			_mm256_store_ps(py, ret);
			_mm256_store_ps(py + 8, ret2);

			px += 16;
			py += 16;
		}
		t.end();
	}
	printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
}

int main(const int argc, const char** argv)
{
	loofline_test(8*1024, 1000000); return 0;
	//課題1
	//行列積和演算AX+Bを計算するプログラムにおいて，行列積と和それぞれの実行時間をタイマーを挟むことで測定せよ．
	if (false)
	{
		std::cout << "exercise 1" << std::endl;
		const int loop = 10;
		const int row = 3;
		const int col = 3;
		Mat_32F a(row, col);
		mat_rand(a, 0, 100);
		Mat_32F x(row, col);
		mat_rand(x, 0, 100);
		Mat_32F b(row, col);
		mat_rand(b, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		CalcTime t;
		for (int i = 0; i < loop; i++)
		{
			t.start();	//時間計測開始
			ret = mat_add(mat_mul(a, x), b);
			t.end();	//時間計測終了
			std::cout << "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;

		a.show();
		x.show();
		b.show();
		ret.show();

		return 0;
	}

	//課題2
	//行列積を計算するプログラムにおいて，コンパイラオプションを変えて計算速度の計測し，その違いを観測せよ．
	if (false)
	{
		std::cout << "exercise 2" << std::endl;
		const int loop = 1000;
		const int row = 64;
		const int col = 64;
		Mat_32F a(row, col);
		mat_rand(a, 0, 100);
		Mat_32F b(row, col);
		mat_rand(b, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		CalcTime t;
		for (int i = 0; i < loop; i++)
		{
			t.start();
			ret = mat_mul(a, b);
			t.end();
			std::cout << "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題3
	//小さな行列に対して，各要素を
	//3x^6+3x^5+3*x^4+3*x^3+3
	//する計算するプログラムを作成し，乗算回数を削減する前と後で計算速度を比較せよ．
	if (false)
	{
		std::cout << "exercise 3" << std::endl;
		const int loop = 100000;
		const int row = 64;
		const int col = 64;
		Mat_32F x(row, col);
		mat_rand(x, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		//before
		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			//そのまま
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				const float v = x.data[i];
				ret.data[i] = 3.f * v * v * v * v * v * v
					+ 3.f * v * v * v * v * v
					+ 3.f * v * v * v * v
					+ 3.f;
			}

			t.end();
			//std::cout << "before: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "before: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//after
		for (int k = 0; k < loop; k++)
		{
			//ホーナー法つかった場合
			//3 * (x * x * x * (x * (x * (x + 1) + 1) + 1));
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				const float v = x.data[i];
				//ret.data[i] = XXXXXXX;
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "after : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題4----ok
	//小さな行列に対して，各要素を下記の定数倍するプログラムを作成し，数式の展開前後で計算速度を比較せよ．
	//(2π+sqrt(5)+0.5^2)x
	if (false)
	{
		std::cout << "exercise 4" << std::endl;
		const int loop = 10000;
		const int row = 64;
		const int col = 64;
		Mat_32F x(row, col);
		mat_rand(x, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		CalcTime t;

		//before
		for (int k = 0; k < loop; k++)
		{
			//毎回計算する場合
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				//XXXX
			}

			t.end();
			//std::cout << "before: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "before: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//after
		for (int k = 0; k < loop; k++)
		{
			//先に計算する場合
			t.start();
			const int size = x.rows * x.cols;
			//XXXXXXXXXXXX //定数値を計算
			for (int i = 0; i < size; i++)
			{
				//計算
				//ret.data[i] = XXXXXXXX;
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "after : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題4.5
	//浮動小数点と整数の加算と積の変換：多くのCPUで効果ないので没．SIMD化などが必要なのと，演算強度が足りない
	if (false)
	{
		std::cout << "exercise 4.5" << std::endl;
		const int loop = 100000000;
		const int row = 64;
		const int col = 64;
		Mat_32F x(row, col);
		mat_rand(x, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		CalcTime t;

		//float演算
		for (int k = 0; k < loop; k++)
		{
			//浮動小数点の乗算
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i += 8)
			{
				//計算
				ret.data[i] = 2.f * x.data[i];
			}

			t.end();
		}
		std::cout << "float mul: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//浮動小数点の加算
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				ret.data[i] = x.data[i] + x.data[i];
			}
			t.end();
		}
		std::cout << "float add : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//int演算
		Mat_32S xi(row, col);
		mat_rand(xi, 0, 100);
		Mat_32S reti(row, col);
		mat_zero(reti);

		for (int k = 0; k < loop; k++)
		{
			//整数の乗算
			t.start();
			const int size = xi.rows * xi.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				reti.data[i] = 2 * xi.data[i];
			}

			t.end();
		}
		std::cout << "int   mul: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//浮動小数点の加算
			t.start();
			const int size = xi.rows * xi.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				reti.data[i] = xi.data[i] + xi.data[i];
			}
			t.end();
		}
		std::cout << "int   add : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題5
	//小さな行列に対して，各要素を3.141592で除算する計算するプログラムを作成し，除算を削減する前と後で計算速度を比較せよ．
	//大きな行列で行うと，効果が少ない可能性があるため注意すること．
	if (false)
	{
		std::cout << "exercise 5" << std::endl;
		const int loop = 100;
		const int row = 64;
		const int col = 64;
		Mat_32F x(row, col);
		mat_rand(x, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);

		CalcTime t;
		//before
		for (int k = 0; k < loop; k++)
		{
			//除算の場合
			t.start();
			for (int j = 0; j < x.rows; j++)
			{
				for (int i = 0; i < x.cols; i++)
				{
					//計算
					//XXXX
				}
			}
			t.end();
			//std::cout << "div: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "div: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//after
		for (int k = 0; k < loop; k++)
		{
			//乗算にした場合
			t.start();
			for (int j = 0; j < x.rows; j++)
			{
				for (int i = 0; i < x.cols; i++)
				{
					//計算
					//XXXX
				}
			}
			t.end();
			//std::cout << "mul: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "mul: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//1重ループの場合
		for (int k = 0; k < loop; k++)
		{
			//1重ループの除算の場合
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				//XXXXXXXXXX
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "div1: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//1重ループの乗算にした場合
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				//XXXXXXXXXX
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "mul1: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題6
	//小さな４つの行列A,B,C,Dに対して，行列の各要素ごとに`(a/b)*(c/d)`を計算するプログラムを作成し，除算を削減する前と後で計算速度を比較せよ．
	if (false)
	{
		std::cout << "exercise 6" << std::endl;
		const int loop = 100;
		const int row = 64;
		const int col = 64;
		Mat_32F a(row, col);
		mat_rand(a, 0, 100);
		Mat_32F b(row, col);
		mat_rand(b, 0, 100);
		Mat_32F c(row, col);
		mat_rand(c, 0, 100);
		Mat_32F d(row, col);
		mat_rand(d, 0, 100);
		Mat_32F ret(row, col);
		mat_zero(ret);


		CalcTime t;
		//before
		for (int k = 0; k < loop; k++)
		{
			//普通に計算
			t.start();
			for (int j = 0; j < ret.rows; j++)
			{
				for (int i = 0; i < ret.cols; i++)
				{
					//計算
					//XXXX
				}
			}
			t.end();
			//std::cout << "before: time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "before: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//after
		for (int k = 0; k < loop; k++)
		{
			//除算を削減した場合
			t.start();
			for (int j = 0; j < ret.rows; j++)
			{
				for (int i = 0; i < ret.cols; i++)
				{
					//計算
					//XXXX
				}
			}
			t.end();
			//std::cout << "after : time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "after: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題7
	//行列の各要素を2乗，3乗，4乗．．．n乗としたときに，mulで作ったものとpowで作ったものの速度を比較せよ．
	//また，nがいくつの時にmulのほうが速くなるのか（それとも常時powのほうが遅い・速いのか）比較せよ．
	if (false)
	{
		std::cout << "exercise 7" << std::endl;
		const int loop = 1000;
		const int row = 64;
		const int col = 64;
		Mat_64F x(row, col);
		mat_rand(x, 0.0, 100.0);
		Mat_64F ret(row, col);
		mat_zero(ret);

		CalcTime t;

		//2乗をpow計算
		int pow_n = 2;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			//powで計算
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				ret.data[i] = pow(x.data[i], pow_n);
			}
			t.end();
		}
		std::cout << pow_n << " : pow time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//2乗をmulで計算
		for (int j = 0; j < loop; j++)
		{
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				ret.data[i] = x.data[i] * x.data[i];
			}
			t.end();
		}
		std::cout << pow_n << " : mul time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//3乗をpow計算
		pow_n = 3;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			//powで計算
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				ret.data[i] = pow(x.data[i], pow_n);
			}
			t.end();
		}
		std::cout << pow_n << " : pow time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//3乗をmulで計算
		for (int j = 0; j < loop; j++)
		{
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				ret.data[i] = x.data[i] * x.data[i] * x.data[i];
			}
			t.end();
		}
		std::cout << pow_n << " : mul time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//4乗をpow計算
		pow_n = 4;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			//powで計算
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				ret.data[i] = pow(x.data[i], pow_n);
			}
			t.end();
		}
		std::cout << pow_n << " : pow time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//4乗をmulで計算
		for (int j = 0; j < loop; j++)
		{
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				ret.data[i] = x.data[i] * x.data[i] * x.data[i] * x.data[i];
			}
			t.end();
		}
		std::cout << pow_n << " : mul time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//n乗をpow計算(nを変える問題．powは64をサンプルとして入力している．)
		pow_n = 64;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			//powで計算
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				ret.data[i] = pow(x.data[i], pow_n);
			}
			t.end();
		}
		std::cout << pow_n << " : pow time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//n乗をmulで計算（nは任意）
		for (int j = 0; j < loop; j++)
		{
			t.start();
			const int size = x.rows * x.cols;
			for (int i = 0; i < size; i++)
			{
				//計算
				ret.data[i] = x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i]
					* x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i]
					* x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i]
					* x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i]
					* x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i]
					* x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i]
					* x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i] * x.data[i]
					;
			}
			t.end();
		}
		std::cout << pow_n << " : mul time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題8
	//2つの行列の和を`unsigned char, short, int, float, double`で計算しそれぞれ比較せよ．
	//なお，大きい行列サイズでないと，効果がでない場合がある．
	if (false)
	{
		std::cout << "exercise 8" << std::endl;
		const int loop = 500;
		const int row = 1024;
		const int col = 1024;

		//unsigend char
		Mat_8U a_8u(row, col);
		mat_rand(a_8u, 0, 100);
		Mat_8U b_8u(row, col);
		mat_rand(b_8u, 0, 100);
		Mat_8U ret_8u(row, col);
		mat_zero(ret_8u);

		CalcTime t;
		for (int i = 0; i < loop; i++)
		{
			t.start();
			//unsigned char
			//XXXX //hint mat_add
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << " 8U : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//short
		Mat_16S a_16s(row, col);
		mat_rand(a_16s, 0, 100);
		Mat_16S b_16s(row, col);
		mat_rand(b_16s, 0, 100);
		Mat_16S ret_16s(row, col);
		mat_zero(ret_16s);

		for (int i = 0; i < loop; i++)
		{
			t.start();
			//short
			//XXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "16S : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//int
		Mat_32S a_32s(row, col);
		mat_rand(a_32s, 0, 100);
		Mat_32S b_32s(row, col);
		mat_rand(b_32s, 0, 100);
		Mat_32S ret_32s(row, col);
		mat_zero(ret_32s);

		for (int i = 0; i < loop; i++)
		{
			t.start();
			//int
			//XXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "32S : time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//float
		Mat_32F a_32f(row, col);
		mat_rand(a_32f, 0, 100);
		Mat_32F b_32f(row, col);
		mat_rand(b_32f, 0, 100);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);

		for (int i = 0; i < loop; i++)
		{
			t.start();
			//float
			//XXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "32F : time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//double
		Mat_64F a_64f(row, col);
		mat_rand(a_64f, 0, 100);
		Mat_64F b_64f(row, col);
		mat_rand(b_64f, 0, 100);
		Mat_64F ret_64f(row, col);
		mat_zero(ret_64f);

		for (int i = 0; i < loop; i++)
		{
			t.start();
			//double
			//XXXX
			t.end();
			//std::cout<< "time: " << t.getLastTime() << " ms" << std::endl;
		}
		std::cout << "64F : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題9
	//intの行列を整数で2倍，浮動小数点で2.f倍,整数を１ビットだけビットシフトすることで2倍する場合の計算速度を比較せよ．
	//また，intの行列を整数で2で除算する場合，浮動小数点で2で除算する場合，浮動小数点の0.5で乗算する場合，１ビットだけビットシフトすることで1/2倍する場合の速度を比較せよ．
	//加えて，floatの行列で，2.0で除算する場合と0.5で乗算する場合を比較せよ．
	//なお，浮動小数点で乗算する場合は整数の場合よりも遅い． 
	//また，大きい行列サイズでないと，効果がでない場合がある．
	if (false)
	{
		std::cout << "exercise 9" << std::endl;
		const int loop = 10000;
		const int row = 1024;
		const int col = 1024;

		//int
		Mat_32S x_32s(row, col);
		mat_rand(x_32s, 0, 100);
		Mat_32S ret_32s(row, col);
		mat_zero(ret_32s);

		CalcTime t;
		//2x mul
		for (int k = 0; k < loop; k++)
		{
			//2倍 乗算
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "2x mul (int)   : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//2.0x mul
		for (int k = 0; k < loop; k++)
		{
			//2.0倍 乗算(double)
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "2x mul (double): time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//2x bit shift
		for (int k = 0; k < loop; k++)
		{
			//2倍 ビットシフト
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "2x bit shift   : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//1/2 div int
		for (int k = 0; k < loop; k++)
		{
			//1/2 除算
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "1/2 div (int)  : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//1/2 div float
		for (int k = 0; k < loop; k++)
		{
			//1/2 除算
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "1/2 div(double): time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//1/2 -> mul 0.5
		for (int k = 0; k < loop; k++)
		{
			//1/2 0.5乗算で実現
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "1/2 mul(double): time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//1/2->bit shift
		for (int k = 0; k < loop; k++)
		{
			//1/2 ビットシフト
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "1/2 bit shift  : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//float
		Mat_32F x_32f(row, col);
		mat_rand(x_32f, 0, 100);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);

		//1/2 div
		for (int k = 0; k < loop; k++)
		{
			//1/2 除算
			t.start();
			const int size = ret_32f.cols * ret_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "float: 1/2 div : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//1/2 -> mul 0.5
		for (int k = 0; k < loop; k++)
		{
			//1/2 0.5乗算
			t.start();
			const int size = ret_32f.cols * ret_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "float: 1/2 mul : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題10
	//floatの行列を3.141倍する場合と，intの行列を3.141倍する場合と，intの行列を3.141倍を固定小数点で行う場合とshortの行列を3.141倍を固定小数点で行う場合で計算し比較せよ．
	//ただし，shortの配列ではオーバーフローに注意せよ．
	//例えば，3.141を10ビットシフトで表現する場合，3.141 * 1024 = 3216であり，short maxは32768であるため，入力の値は最大10までしかとることができない．課題のコードでは0～100の乱数であるため，適宜シフトの量を工夫せよ．
	if (false)
	{
		std::cout << "exercise 10" << std::endl;
		const int loop = 10000;
		const int row = 1024;
		const int col = 1024;

		CalcTime t;

		Mat_32F x_32f(row, col);
		mat_rand(x_32f, 0, 100);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);
		//float
		for (int k = 0; k < loop; k++)
		{
			//浮動小数点
			t.start();
			const int size = ret_32f.cols * ret_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "float mul(float): time (avg): " << t.getAvgTime() << " ms" << std::endl;

		Mat_32S x_32s(row, col);
		mat_rand(x_32s, 0, 100);
		Mat_32S ret_32s(row, col);
		mat_zero(ret_32s);
		//int floating point mul
		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "int   mul(float): time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//int fix
		for (int k = 0; k < loop; k++)
		{
			//固定小数点
			t.start();
			const int size = ret_32s.cols * ret_32s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "int   mul(fixed): time (avg): " << t.getAvgTime() << " ms" << std::endl;

		Mat_16S x_16s(row, col);
		mat_rand(x_16s, 0, 100);
		Mat_16S ret_16s(row, col);
		mat_zero(ret_16s);
		//short fix
		for (int k = 0; k < loop; k++)
		{
			//固定小数点
			t.start();
			const int size = ret_16s.cols * ret_16s.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "short mul(fixed): time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題11
	//floatの行列への定数値(3.141592f)の四則演算(加算，乗算，除算)と，`sqrt sin, cos, exp, log`関数の適用した場合と計算時間を比較せよ．
	//また，`sin, cos, exp, log, sqrt`計算はテーブル参照も作成した場合についても比較せよ．
	//なお，環境によっては，演算したほうが速い演算もある可能性がある．
	if (false)
	{
		std::cout << "exercise 11" << std::endl;
		const int loop = 10000;
		const int row = 128;
		const int col = 128;

		CalcTime t;

		Mat_32F x_32f(row, col);
		mat_rand(x_32f, 0.f, 255.f);
		Mat_32F ret_32f(row, col);
		mat_zero(ret_32f);

		//四則演算
		for (int k = 0; k < loop; k++)
		{
			//加算
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "add: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		for (int k = 0; k < loop; k++)
		{
			//乗算
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "mul: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//除算
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "div: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//sqrt
		for (int k = 0; k < loop; k++)
		{
			//sqrt関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "sqr: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//sin
		for (int k = 0; k < loop; k++)
		{
			//sin関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "sin: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//cos
		for (int k = 0; k < loop; k++)
		{
			//cos関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "cos: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//exp
		for (int k = 0; k < loop; k++)
		{
			//exp関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "exp: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//log
		for (int k = 0; k < loop; k++)
		{
			//log関数
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "log: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		std::cout << std::endl;

		float LUT[256];
		//sqrt LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			LUT[i] = (float)sqrt(i);// 以下の課題ようのヒントのために，埋めてあります．
		}
		for (int k = 0; k < loop; k++)
		{
			//sqrt LUT
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				ret_32f.data[i] = LUT[(int)x_32f.data[i]];// 以下の課題ようのヒントのために，埋めてあります．
			}
			t.end();
		}
		std::cout << "sqr LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//sin LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXX
		}
		for (int k = 0; k < loop; k++)
		{
			t.start();
			//sin LUT
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX				
			}
			t.end();
		}
		std::cout << "sin LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//cos LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//cos LUT
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "cos LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//exp LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//exp LUT
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "exp LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//log LUT
		for (int i = 0; i < 256; i++)
		{
			//LUT作成
			//XXXX
		}
		for (int k = 0; k < loop; k++)
		{
			//log LUT
			t.start();
			const int size = x_32f.cols * x_32f.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "log LUT: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題(省略可)
	//ループ構造変換の課題にあるので省略した
	if (false)
	{
		std::cout << "exercise (省略可1)" << std::endl;
		CalcTime t;
		const int loop = 1000;
		const int size = 1024;
		Mat_32F a(size, size);
		Mat_32F b(size, size);

		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < size; j++)
			{
				for (int i = 0; i < size; i++)
				{
					if (i == j)
						a.data[size * j + i] = 1.f;
					else
						a.data[size * j + i] = b.data[size * j + i];
				}
			}
			t.end();
		}
		std::cout << "simple loop: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < size; j++)
			{
				for (int i = 0; i < size; i++)
				{
					a.data[size * j + i] = b.data[size * j + i];
				}
			}
			for (int i = 0; i < size; i++)
			{
				a.data[size * i + i] = 1.f;
			}
			t.end();
		}

		std::cout << "loop unswitching: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		return 0;
	}

	//課題12
	//小さな行列A,Bの各要素を任意のradianだけ回転させて，x,yにして格納するプログラムを記述し，inline展開の有無で速度がどのように変わるか計測せよ．
	//また，関数をべた書きした場合とも比較せよ．
	//ただし，-O3のオプションを付けると強制的にinline展開される可能性がある．
	//inline void rot(double a, double b, double &x, double &y, double radian)
	//{
	//	x = a * cos(radian);
	//	y = b * sin(radian);
	//}
	if (false)
	{
		std::cout << "exercise 12" << std::endl;
		const int loop = 10000;
		const int row = 64;
		const int col = 64;

		CalcTime t;

		Mat_64F a(row, col);
		mat_rand(a, 0.0, 100.0);
		Mat_64F b(row, col);
		mat_rand(b, 0.0, 100.0);
		Mat_64F x(row, col);
		mat_zero(x);
		Mat_64F y(row, col);
		mat_zero(y);

		const float radian = 2.2f;

		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = a.cols * a.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
			}
			t.end();
		}
		std::cout << "func time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//べた書き（関数の中身をループないに書く）
		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = a.cols * a.rows;
			for (int i = 0; i < size; i++)
			{
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "loop time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題13
	//行列A，Bの各要素の乗算を行うときに，結果を行列Cに格納する場合と行列Aに上書きする場合との計算時間をせよ．
	if (false)
	{
		std::cout << "exercise 13" << std::endl;
		const int loop = 10000;
		const int row = 64;
		const int col = 64;
		//c = a x b
		{
			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				//C=A*B
				Mat_64F a(row, col);
				Mat_64F b(row, col);
				Mat_64F c(row, col);
				const int size = a.rows * a.cols;
				for (int i = 0; i < size; i++)
				{
					//XXXX
				}
				t.end();
			}
			std::cout << "c=a*b: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		}

		//a = a x b
		{

			CalcTime t;
			for (int k = 0; k < loop; k++)
			{
				t.start();
				//A=AxB
				Mat_64F a(row, col);
				Mat_64F b(row, col);
				const int size = a.rows * a.cols;
				for (int i = 0; i < size; i++)
				{
					//XXXX
				}
				t.end();
			}
			std::cout << "a=a*b: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		}

		return 0;
	}

	//課題14
	//上記の0で初期化するコードをループの順序を変えてどちらが速いか計測して検証せよ．
	//また，行列積のコードのループの順序を変えてどれが速いか計測して検証せよ．
	//なお，行列積のコードは記述済みであり，どれも違いがないことを正解との差分の二乗和（MSE）で評価している．
	if (false)
	{
		std::cout << "exercise 14" << std::endl;
		const int loop = 1000;
		const int row = 128;
		const int col = 128;

		CalcTime t;

		Mat_32F x(row, col);

		//col, row
		for (int k = 0; k < loop; k++)
		{
			int i = 0, j = 0;
			t.start();
			//xxxxxx	hint: i loop row
			{
				//xxxxxx hint:j loop col
				{
					x.data[col * i + j] = 0.f;
				}
			}
			t.end();
		}
		std::cout << "col-row: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//row, col
		for (int k = 0; k < loop; k++)
		{
			int i = 0, j = 0;
			t.start();
			//xxxxxx hint:j loop col
			{
				//xxxxxx	hint: i loop row
				{
					x.data[col * i + j] = 0.f;
				}
			}
			t.end();
		}
		std::cout << "row-col: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		//行列積
		const int size = 128;
		Mat_32F a(size, size);
		Mat_32F b(size, size);
		Mat_32F c(size, size);
		Mat_32F ans(size, size);
		mat_rand(a, 0, 10);
		mat_rand(b, 0, 10);

		//for anser
		mat_zero(ans);
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < size; ++j)
			{
				for (int k = 0; k < size; ++k)
				{
					//ans[i][j] += a[i][k] * b[k][j];
					ans.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
				}
			}
		}

		//i, j, k
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int i = 0; i < size; ++i)
			{
				float* pc = c.data + i * size;
				float* pa = a.data + i * size;
				for (int j = 0; j < size; ++j)
				{
					float* pb = b.data + j;
					float* pcj = &pc[j];
					for (int k = 0; k < size; ++k)
					{
						//c[i][j] += a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						*pcj += pa[k] * pb[k * size];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "i-j-k: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//i, k, j
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int i = 0; i < size; ++i)
			{
				float* pc = c.data + i * size;
				float* pa = a.data + i * size;
				for (int k = 0; k < size; ++k)
				{
					float* pb = b.data + k * size;
					const float pak = pa[k];
					for (int j = 0; j < size; ++j)
					{
						//c[i][j] += a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						//pc[j] += (pak * pb[j]);//コンパイラ最適化でfmaがかかる
						pc[j] += pa[k] * pb[j];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "i-k-j: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//j, i, k
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int j = 0; j < size; ++j)
			{
				float* pb = b.data + j;
				for (int i = 0; i < size; ++i)
				{
					float* pc = c.data + i * size;
					float* pa = a.data + i * size;
					for (int k = 0; k < size; ++k)
					{
						//c[i][j] = c[i][j] + a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						pc[j] += pa[k] * pb[k * size];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "j-i-k: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//j, k, i
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int j = 0; j < size; ++j)
			{
				float* pc = c.data + j;
				float* pb = b.data + j;
				for (int k = 0; k < size; ++k)
				{
					float* pa = a.data + k;
					const float pbj = pb[k * size];
					for (int i = 0; i < size; ++i)
					{
						//c[i][j] = c[i][j] + a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						pc[i * size] += pa[i * size] * pbj;
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "j-k-i: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//k, i, j
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int k = 0; k < size; ++k)
			{
				float* pb = b.data + k * size;
				for (int i = 0; i < size; ++i)
				{
					float* pc = c.data + i * size;
					float* pa = a.data + i * size;
					const float pak = pa[k];
					for (int j = 0; j < size; ++j)
					{
						//c[i][j] = c[i][j] + a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						//pc[j] += pak * pb[j];//コンパイル最適化でFMAがかかる
						pc[j] += pa[k] * pb[j];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "k-i-j: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		//k, j, i
		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			for (int k = 0; k < size; ++k)
			{
				float* pa = a.data + k;
				float* pb = b.data + k * size;
				for (int j = 0; j < size; ++j)
				{
					float* pc = c.data + j;
					const float pbj = pb[j];
					for (int i = 0; i < size; ++i)
					{
						//c[i][j] = c[i][j] + a[i][k] * b[k][j];
						//c.data[i * size + j] += a.data[i * size + k] * b.data[k * size + j];
						pc[i * size] += pa[i * size] * pbj;
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "k-j-i: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		std::cout << "diff from ans: " << mat_diff(ans, c) << std::endl;

		return 0;
	}

	//課題15
	//アンローリングの段数を2,4,8,16,32,...と変更することで，速度がどのように変わるか計測せよ．
	if (false)
	{
		std::cout << "exercise 15" << std::endl;
		const int loop = 1000;
		const int size = 65535;

		CalcTime t;

		float* x = new float[size];
		float* y = new float[size];
		const float a = 2.f;
		const float b = 1.f;

		for (int i = 0; i < size; i++)
		{
			x[i] = rand_32f(0.f, 100.f);
		}

		//unrolling 1
		for (int j = 0; j < loop; j++)
		{
			//unrolling 1
			t.start();
			for (int i = 0; i < size; i++)
			{
				y[i + 0] = a * x[i + 0] + b;
			}
			t.end();
		}
		std::cout << "no unrolling: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 2
		for (int j = 0; j < loop; j++)
		{
			//unrolling 2
			t.start();
			for (int i = 0; i < size; i += 2)
			{
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "unrolling  2: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 4
		for (int j = 0; j < loop; j++)
		{
			//unrolling 4
			t.start();
			for (int i = 0; i < size; i += 4)
			{
				//XXXX
				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "unrolling  4: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 8
		for (int j = 0; j < loop; j++)
		{
			//unrolling 8
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "unrolling  8: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 16
		for (int j = 0; j < loop; j++)
		{
			//unrolling 16
			t.start();
			for (int i = 0; i < size; i += 16)
			{
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "unrolling 16: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 32
		for (int j = 0; j < loop; j++)
		{
			//unrolling 32
			t.start();
			for (int i = 0; i < size; i += 32)
			{
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "unrolling 32: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//unrolling 64
		for (int j = 0; j < loop; j++)
		{
			//unrolling 64
			t.start();
			for (int i = 0; i < size; i += 64)
			{
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "unrolling 64: time (avg): " << t.getAvgTime() << " ms" << std::endl;
		delete[] x;
		delete[] y;
		return 0;
	}

	//課題16
	//上記のプログラムを実装し，ループピーリングの有無で速度がどのように変わるか計測せよ．
	if (false)
	{
		std::cout << "exercise 16" << std::endl;
		const int loop = 100;
		const int size = 65535;
		float* x = new float[size];
		float* y = new float[size];

		CalcTime t;

		//init
		for (int i = 0; i < size; i++)
		{
			x[i] = rand_32f(0, 100);
			y[i] = rand_32f(0, 100);
		}

		//Original
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; ++i)
			{
				if (i == 0)
				{
					y[i] = (x[i] + x[i + 1]) / 2.f;
				}
				else if (i == size - 1)
				{
					y[i] = (x[i - 1] + x[i]) / 2.f;
				}
				else
				{
					y[i] = (x[i - 1] + x[i] + x[i + 1]) / 3.f;
				}
			}
			t.end();
		}
		std::cout << "original: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//Loop peeling
		for (int j = 0; j < loop; j++)
		{
			t.start();
			{
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "loop peeling: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		delete[] x;
		delete[] y;
		return 0;
	}

	//課題17
	//Mat_32Sの配同士の加算(a=a+b)をループつぶしをするか否かで計算時間を比較せよ．
	if (false)
	{
		std::cout << "exercise 17" << std::endl;
		const int loop = 100000;
		const int width = 128;
		const int height = 128;
		//heightとwidthを大きくする場合は，注意すること．
		//下記にa[height][width]という宣言があり，大きすぎる配列はヒープ領域からメモリを取ってこれないため．

		CalcTime t;

		Mat_32S ma(height, width);
		Mat_32S mb(height, width);

		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < height; j++)
			{
				for (int i = 0; i < width; i++)
				{
					ma.data[j * width + i] = ma.data[j * width + i] + mb.data[j * width + i];
				}
			}
			t.end();
		}
		std::cout << "mat before: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//下記の2次元配列の場合のサンプルを参考にするとよい．
		for (int k = 0; k < loop; k++)
		{
			t.start();
			const int size = width * height;
			//XXXXXXX
			//XXXXXXX
			for (int i = 0; i < size; i++)
			{
				//XXXXXXX
			}
			t.end();
		}
		std::cout << "mat after : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		/*
		//2次元配列の場合のsample code
		//以下はコード記述済み．コンパイル時最適化のせいで，おそらく違いはない．興味がある人はgcc main.cpp -sでアセンブラ出力して確認すること．
		int a[height][width];
		int b[height][width];
		// before
		for (int k = 0; k < loop; k++)
		{
			t.start();
			for (int j = 0; j < height; j++)
			{
				for (int i = 0; i < width; i++)
				{
					a[j][i] = a[j][i] + b[j][i];
				}
			}
			t.end();
		}
		std::cout << "array before: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//after
		const int size = width * height;
		for (int k = 0; k < loop; k++)
		{
			t.start();
			int* pa = &a[0][0];
			int* pb = &b[0][0];
			for (int i = 0; i < size; i++)
			{
				*pa++ = *pa + *pb++;
			}
			t.end();
		}
		std::cout << "array after : time (avg): " << t.getAvgTime() << " ms" << std::endl;
		*/
		return 0;
	}

	//課題18
	//上記のコードを実行し，並列に動作していることを確認せよ．
	//また，並列化を有効にする場合としない場合の計算時間を比較せよ．
	if (false)
	{
		std::cout << "exercise 18" << std::endl;
		//XXXX
		for (int i = 0; i < 100; i++)
		{
			std::cout << i << std::endl; //並列化したい処理
		}
		return 0;
	}

	//課題19
	//総和を計算するコードで，reduction指定子を使用する場合としない場合で計算結果がどのようになるか比較せよ．
	if (false)
	{
		std::cout << "exercise 19" << std::endl;
		//並列化なし（逐次実行）
		std::cout << "no parallelization" << std::endl;
		for (int j = 0; j < 10; j++)
		{
			int sum = 0;
			for (int i = 0; i < 100; i++)
			{
				sum += i;
			}
			std::cout << sum << std::endl;
		}

		//omp使用
		std::cout << "parallelization" << std::endl;
		for (int j = 0; j < 10; j++)
		{
			int sum = 0;
			//XXXX
			for (int i = 0; i < 100; i++)
			{
				sum += i;
			}
			std::cout << sum << std::endl;
		}

		//omp reduction使用
		std::cout << "parallelization with reduction" << std::endl;
		for (int j = 0; j < 10; j++)
		{
			int sum = 0;
			//XXXX
			for (int i = 0; i < 100; i++)
			{
				sum += i;
			}
			std::cout << sum << std::endl;
		}
		return 0;
	}

	//課題20
	//二つの行列の各要素の積を計算するコードで，スレッド数を変更して，計算時間がどのように推移するのかを確認せよ．
	//なお，スレッド数は，計算機のコア数以上の物まで指定せよ．
	//8コア16スレッドのPCでは，16コアよりも大きいスレッド数（例えば32までなど）までを指定せよ．
	if (false)
	{
		std::cout << "exercise 20" << std::endl;
		const int loop = 1000;
		const int size = 128;

		CalcTime t;

		Mat_32F a(size, size);
		mat_rand(a, 0, 100);
		Mat_32F b(size, size);
		mat_rand(b, 0, 100);
		Mat_32F c(size, size);

		for (int l = 0; l < loop; l++)
		{
			mat_zero(c);
			t.start();
			//#pragma omp parallel for num_threads(n)で並列化，nに任意の整数を入れる
			//XXXX
			for (int i = 0; i < size; ++i)
			{
				for (int k = 0; k < size; ++k)
				{
					for (int j = 0; j < size; ++j)
					{
						c.data[i * c.cols + j] = a.data[i * a.cols + k] * b.data[k * b.cols + j];
					}
				}
			}
			t.end();
			//std::cout << "time : " << t.getLastTime() << std::endl;
		}
		std::cout << "time (avg): " << t.getAvgTime() << " ms" << std::endl;
		return 0;
	}

	//課題21
	//四則演算のコードを書いてprintfデバッグで確認せよ．
	if (false)
	{
		std::cout << "exercise 21" << std::endl;
		const __m256 a = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
		const __m256 b = _mm256_set_ps(15, 14, 13, 12, 11, 10, 9, 8);
		//逆順で入力(setr)こちらのほうが使いやすい
		const __m256 c = _mm256_setr_ps(0, 1, 2, 3, 4, 5, 6, 7);
		const __m256 d = _mm256_setr_ps(8, 9, 10, 11, 12, 13, 14, 15);
		__m256 e = _mm256_setzero_ps();
		//e = a + b
		//XXXXXX
		std::cout << "add a b: ";
		print_m256(e);
		//e = c + d
		//XXXXXX
		std::cout << "add c d: ";
		print_m256(e);//出力は上と同じのはず．

		//print_m256の関数を使わない場合はこうやって出力
		float temp[8];
		_mm256_store_ps(&temp[0], e);
		std::cout << "cout ex: ";
		for (int i = 0; i < 8; i++) std::cout << temp[i] << ", ";
		std::cout << std::endl;

		//減算
		//e = a - b
		//XXXXXX
		std::cout << "sub a b: ";
		print_m256(e);
		//e = c - d
		//XXXXXX
		std::cout << "sub c d: ";
		print_m256(e);

		//乗算
		//e = a * b
		//XXXXXX
		std::cout << "mul a b: ";
		print_m256(e);
		//e = c * d
		//XXXXXX
		std::cout << "mul c d: ";
		print_m256(e);

		//除算
		//e = a / b
		//XXXXXX
		std::cout << "div a b: ";
		print_m256(e);
		//e = c / d
		//XXXXXX
		std::cout << "div c d: ";
		print_m256(e);

		return 0;
	}

	//課題22
	//配列a,x,bに対して，`(((a*x+b)*x+b)*x+b)*x+b `の計算を配列ｃに格納するコードをmul/addで記述するものとFMAを使うもので記述し，FMAが速くなることを示せ．
	//なお，上記の関数は以下に等しい．
	//a=_mm256_fmadd_ps(a,b,c);
	//a=_mm256_fmadd_ps(a,b,c);
	//a=_mm256_fmadd_ps(a,b,c);
	//a=_mm256_fmadd_ps(a,b,c);
	//これは，単純にFMAが1度だとメモリで律速するこのコードでは計算速度の差が出にくいためである．差が小さければ，より演算を増やせば良い．
	//また，シングルコアのGFLOPSと演算強度[FLOPS/BYTE]をFMA命令で計算するサンプルが書いてある．これをプロットし，グラフとして図示せよ．
	//if (false)
	{
		std::cout << "exercise 22" << std::endl;
		const int loop = 10000000;
		const int col = 32;
		const int row = 32;
		const int size = col * row;

		CalcTime t;

		Mat_32F a(row, col);
		Mat_32F x(row, col);
		Mat_32F b(row, col);
		Mat_32F c(row, col);
		mat_rand(a, 0, 1);
		mat_rand(b, 0, 1);
		mat_rand(x, 0, 1);
		mat_zero(c);

		/*
		//mul, add
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* pa = a.data;
			float* px = x.data;
			float* pb = b.data;
			float* pc = c.data;
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(pa + i);
				const __m256 mx = _mm256_load_ps(px + i);
				const __m256 mb = _mm256_load_ps(pb + i);

				//mul,addを使って
				__m256 temp;
				temp = _mm256_add_ps(_mm256_mul_ps(ma, mx), mb);
				temp = _mm256_add_ps(_mm256_mul_ps(temp, mx), mb);
				temp = _mm256_add_ps(_mm256_mul_ps(temp, mx), mb);
				temp = _mm256_add_ps(_mm256_mul_ps(temp, mx), mb);
				_mm256_store_ps(pc + i, temp);
			}
			t.end();
		}
		std::cout << "mul-add: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//fma
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* pa = a.data;
			float* px = x.data;
			float* pb = b.data;
			float* pc = c.data;
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(pa + i);
				const __m256 mx = _mm256_load_ps(px + i);
				const __m256 mb = _mm256_load_ps(pb + i);

				//fmaを使って
				__m256 temp;
				temp = _mm256_fmadd_ps(ma, mx, mb);
				temp = _mm256_fmadd_ps(temp, mx, mb);
				temp = _mm256_fmadd_ps(temp, mx, mb);
				temp = _mm256_fmadd_ps(temp, mx, mb);
				_mm256_store_ps(pc + i, temp);
			}
			t.end();
		}
		std::cout << "fma    : time (avg): " << t.getAvgTime() << " ms" << std::endl;
		*/
		std::cout << "======" << std::endl;
		std::cout << "GFLOPS" << std::endl;
		std::cout << "======" << std::endl;
		//FLOPS計算
		//x+1: (1*x+1)
		//x*x+x+1: x*(x+1)+1
		//x*x*x+x*x+x+1: x*(x*(x+1))+1
		int simdsize = size / 8;

		//n=...
		std::cout << "..." << std::endl;
		int n = 1;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x.data;
			float* pc = c.data;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				//1次
				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1

				_mm256_store_ps(pc, ret);

				px += 8;
				pc += 8;
			}
			t.end();
		}
		std::cout << "fma " << n << "th, " << n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000) << ", GFLOPS, " << n * 2.0 / (4 * 4) << ", FLOPS/BYTE" << std::endl;

		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x.data;
			float* pc = c.data;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize / 2; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				//1次
				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1

				_mm256_store_ps(pc, ret);
				_mm256_store_ps(pc + 8, ret2);

				px += 16;
				pc += 16;
			}
			t.end();
		}
		std::cout << "fma " << n << "th, " << n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000) << ", GFLOPS, " << n * 2.0 / (4 * 4) << ", FLOPS/BYTE" << std::endl;

		n = 2;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x.data;
			float* pc = c.data;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				//1次
				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				_mm256_store_ps(pc, ret);

				px += 8;
				pc += 8;
			}
			t.end();
		}
		std::cout << "fma " << n << "th, " << n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000) << ", GFLOPS, " << n * 2.0 / (4 * 4) << ", FLOPS/BYTE" << std::endl;

		n = 5;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x.data;
			float* pc = c.data;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				//20次
				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5

				_mm256_store_ps(pc, ret);

				px += 8;
				pc += 8;
			}
			t.end();
		}
		std::cout << "fma " << n << "th, " << n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000) << ", GFLOPS, " << n * 2.0 / (4 * 4) << ", FLOPS/BYTE" << std::endl;

		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x.data;
			float* pc = c.data;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize / 2; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				//20次
				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5

				_mm256_store_ps(pc, ret);
				_mm256_store_ps(pc + 8, ret2);

				px += 16;
				pc += 16;
			}
			t.end();
		}
		std::cout << "fma " << n << "th, " << n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000) << ", GFLOPS, " << n * 2.0 / (4 * 4) << ", FLOPS/BYTE" << std::endl;

		n = 10;
		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x.data;
			float* pc = c.data;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);

				//20次
				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10

				_mm256_store_ps(pc, ret);

				px += 8;
				pc += 8;
			}
			t.end();
		}
		std::cout << "fma " << n << "th, " << n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000) << ", GFLOPS, " << n * 2.0 / (4 * 4) << ", FLOPS/BYTE" << std::endl;

		for (int j = 0; j < loop; j++)
		{
			t.start();
			float* px = x.data;
			float* pc = c.data;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = simdsize / 2; i != 0; i--)
			{
				const __m256 mx = _mm256_load_ps(px);
				const __m256 mx2 = _mm256_load_ps(px + 8);

				//20次
				__m256 ret = _mm256_fmadd_ps(mones, mx, mones);//1
				__m256 ret2 = _mm256_fmadd_ps(mones, mx2, mones);//1
				ret = _mm256_fmadd_ps(mx, ret, mones);//2
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//2
				ret = _mm256_fmadd_ps(mx, ret, mones);//3
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//3
				ret = _mm256_fmadd_ps(mx, ret, mones);//4
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//4
				ret = _mm256_fmadd_ps(mx, ret, mones);//5
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//5
				ret = _mm256_fmadd_ps(mx, ret, mones);//6
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//6
				ret = _mm256_fmadd_ps(mx, ret, mones);//7
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//7
				ret = _mm256_fmadd_ps(mx, ret, mones);//8
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//8
				ret = _mm256_fmadd_ps(mx, ret, mones);//9
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//9
				ret = _mm256_fmadd_ps(mx, ret, mones);//10
				ret2 = _mm256_fmadd_ps(mx2, ret2, mones);//10

				_mm256_store_ps(pc, ret);
				_mm256_store_ps(pc + 8, ret2);

				px += 16;
				pc += 16;
			}
			t.end();
		}
		std::cout << "fma " << n << "th, " << n * 2.0 * size / (t.getAvgTime() * 0.001) / (1000 * 1000 * 1000) << ", GFLOPS, " << n * 2.0 / (4 * 4) << ", FLOPS/BYTE" << std::endl;

		return 0;
	}

	//課題23
	//divとrcp,sqrtとrsqrtの実行速度を比較せよ．
	if (false)
	{
		std::cout << "exercise 23" << std::endl;
		const int loop = 10000;
		const int size = 1024;
#ifdef __GNUC__
		float __attribute__((aligned(32))) a[size];
		float __attribute__((aligned(32))) b[size];
#elif _MSC_VER
		__declspec(align(32)) float a[size];
		__declspec(align(32)) float b[size];
#endif

		for (int i = 0; i < size; i++)
		{
			a[i] = rand_32f(0, 100);
			b[i] = rand_32f(0, 100);
		}

		CalcTime t;
		//div
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mb = _mm256_load_ps(b + i);

				__m256 temp;
				//divを使って
				//XXXX
			}
			t.end();
		}
		std::cout << "div: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//rcp
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mb = _mm256_load_ps(b + i);

				__m256 temp;
				//rcpをつかって
				//XXXX
			}
			t.end();
		}
		std::cout << "rcp: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//sqrt
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mb = _mm256_load_ps(b + i);

				__m256 temp;
				//sqrtを使って
				//XXXX
			}
			t.end();
		}
		std::cout << " sqrt: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		//rsqrt
		for (int j = 0; j < loop; j++)
		{
			t.start();
			for (int i = 0; i < size; i += 8)
			{
				const __m256 ma = _mm256_load_ps(a + i);
				const __m256 mb = _mm256_load_ps(b + i);

				__m256 temp;
				//rsqrtを使って
				//XXXX
			}
			t.end();
		}
		std::cout << "rsqrt: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}


	//課題24
	//haddとdpで要素の総和を取るプログラムを作成せよ．
	//また，それぞれの計算時間を比較せよ ．
	//（この課題は，後にリダクションの最適化でもう一度登場します．）
	if (false)
	{
		std::cout << "exercise 24" << std::endl;
		const int loop = 100000;
#ifdef __GNUC__
		float __attribute__((aligned(32))) a[8];
		float __attribute__((aligned(32))) b[8];
#elif _MSC_VER
		__declspec(align(32)) float a[8];
		__declspec(align(32)) float b[8];
#endif

		for (int i = 0; i < 8; i++)
		{
			a[i] = rand_32f(0, 100);
			b[i] = rand_32f(0, 100);
		}

		CalcTime t;
		//hadd
		for (int i = 0; i < loop; i++)
		{
			t.start();
			__m256 ma = _mm256_load_ps(a);
			//haddを使って
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			t.end();
		}
		std::cout << "hadd: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		const __m256 one = _mm256_set1_ps(1);
		//dp
		for (int i = 0; i < loop; i++)
		{
			t.start();
			__m256 mb = _mm256_load_ps(b);
			//dpを使って
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			//XXXX
			t.end();
		}
		std::cout << "dp: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題25
	//行列aにおいて要素の値があるしきい値以上の場合だけ3乗し，それ以外は何もしない処理をベクトル化実装せよ．
	//if(a[i]>=threshold) a[i]=a[i]*a[i]*a[i];
	if (false)
	{
		std::cout << "exercise 25" << std::endl;
		const int size = 8;
		const float th = 50;

#ifdef __GNUC__
		float __attribute__((aligned(32))) a[size];
#elif _MSC_VER
		__declspec(align(32)) float a[size];
#endif

		for (int i = 0; i < size; i++)
		{
			a[i] = rand_32f(0, 100);
		}

		const __m256 mth = _mm256_set1_ps(th);
		for (int i = 0; i < size; i += 8)
		{
			const __m256 ma = _mm256_load_ps(a + i);

			__m256 temp;
			//cmp, mul, blendvを使って
			//XXXX
			//XXXX
			//XXXX
			_mm256_store_ps(a + i, temp);
		}
		for (int i = 0; i < size; i++)
		{
			std::cout << a[i] << ", ";
			if (i != 0 && i % size == 0)
				std::cout << std::endl;
		}
		std::cout << std::endl;
		return 0;
	}

	//課題26
	//上記のコードのように，SIMD命令を使う場合におけるループアンローリングを8，16，32，64と行い，計算時間を比較せよ．
	if (false)
	{
		std::cout << "exercise 26" << std::endl;
		const int loop = 1000;
		const int size = 1024 * 3;
#ifdef __GNUC__
		__attribute__((aligned(32))) float a[size];
		__attribute__((aligned(32))) float b[size];
		__attribute__((aligned(32))) float c[size];
#elif _MSC_VER
		__declspec(align(32)) float a[size];
		__declspec(align(32)) float b[size];
		__declspec(align(32)) float c[size];
#endif

		//init
		for (int i = 0; i < size; i++)
		{
			a[i] = rand_32f(0, 100);
			b[i] = rand_32f(0, 100);
			c[i] = 0;
		}

		CalcTime t;
		// unrolling 8
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 8
			for (int i = 0; i < size; i += 8)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
				_mm256_store_ps(c + i, _mm256_add_ps(ma, mb));
			}
			t.end();
		}
		std::cout << "  8: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		// unrolling 16
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 16
			for (int i = 0; i < size; i += 16)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
				//XXXX

				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << " 16: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		// unrolling 32
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 32
			for (int i = 0; i < size; i += 32)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << " 32: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		// unrolling 64
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 64
			for (int i = 0; i < size; i += 64)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << " 64: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		// unrolling 128
		for (int j = 0; j < loop; j++)
		{
			t.start();
			// unrolling 128
			for (int i = 0; i < size; i += 128)
			{
				__m256 ma = _mm256_load_ps(a + i);
				__m256 mb = _mm256_load_ps(b + i);
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX

				//XXXX
				//XXXX
				//XXXX
			}
			t.end();
		}
		std::cout << "128: time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	//課題27
	//上の関数を用いて，データ構造の相互変換を確認せよ． ４ｘ４のdouble型データの転置を手書きで考えてみよ．
	//この問題は，自力で答えを埋めると時間がかかり過ぎるためプログラムの答えはすでに埋めてある．
	//そのため，この動作を追って転置ができていることを手書きで確認するだけで良い．
	if (false)
	{
		std::cout << "exercise 27" << std::endl;
		const int size = 64;
#ifdef __GNUC__
		__attribute__((aligned(32))) float a[size];
		__attribute__((aligned(32))) float b[size];
#elif _MSC_VER
		__declspec(align(32)) float a[size];
		__declspec(align(32)) float b[size];
#endif
		for (int i = 0; i < size; i++)
		{
			a[i] = i;
			b[i] = 0;
		}
		__m256 ma[8], mb[8];
		for (int i = 0; i < 8; i++)
		{
			ma[i] = _mm256_load_ps((float*)(&a[i * 8]));
			mb[i] = _mm256_setzero_ps();
		}

		for (int i = 0; i < 8; i++)
		{
			print_m256(ma[i]);
		}

		//転置
//XXXX

		for (int i = 0; i < 8; i++)
		{
			print_m256(mb[i]);
		}
		return 0;
	}

	//課題28
	//__m256i（int）型を_m256（float）型に変換せよ．
	if (false)
	{
		std::cout << "exercise 28" << std::endl;
		__m256i m32i = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
		__m256 m32f = _mm256_setzero_ps();

		print_m256(m32f);
		//cvtを使って
//XXXX
		print_m256(m32f);
	}


	//課題29
	//上記のコードを実行し，計算時間を計測せよ．
	//また，スカラ実装，スカラ実装＋並列化，SIMD実装のみを作成し，計算時間を比較せよ．
	if (false)
	{
		std::cout << "exercise 29" << std::endl;
		const int loop = 1000;
		const int size = 128;
		Mat_32F a(size, size);
		Mat_32F b(size, size);
		Mat_32F c(size, size);

		//init
		mat_rand(a, 0, 100);
		mat_rand(b, 0, 100);
		mat_zero(c);

		CalcTime t;
		for (int k = 0; k < loop; k++)
		{
			//スカラー実装
			t.start();
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "scalar    : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		for (int k = 0; k < loop; k++)
		{
			//スカラー，並列化実装
			t.start();
			//XXXX
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i++)
				{
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "scalar+omp: time (avg): " << t.getAvgTime() << " ms" << std::endl;


		for (int k = 0; k < loop; k++)
		{
			//SIMD実装
			t.start();
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i += 8)
				{
					//XXXX
					//XXXX
					//XXXX
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "SIMD      : time (avg): " << t.getAvgTime() << " ms" << std::endl;


		for (int k = 0; k < loop; k++)
		{
			//SIMD，並列化実装
			t.start();
			//XXXX
			for (int j = 0; j < size; ++j)
			{
				for (int i = 0; i < size; i += 8)
				{
					//XXXX
					//XXXX
					//XXXX
					//XXXX
				}
			}
			t.end();
		}
		std::cout << "SIMD+omp  : time (avg): " << t.getAvgTime() << " ms" << std::endl;

		return 0;
	}

	std::cout << "no select" << std::endl;
	return 0;
}

//課題12用
inline void rot(double a, double b, double& x, double& y, double radian)
{
	x = a * cos(radian);
	y = b * sin(radian);
}

//課題27用
inline void _mm256_transpose_8x8_ps(__m256* dst, const __m256* src)
{
	__m256  tmp[8], tmpp[8];

	for (int i = 0; i < 8; i += 2)
	{
		tmp[i + 0] = _mm256_unpacklo_ps(src[i], src[i + 1]);
		tmp[i + 1] = _mm256_unpackhi_ps(src[i], src[i + 1]);
	}
	for (int i = 0; i < 8; i += 4)
	{
		tmpp[i + 0] = _mm256_shuffle_ps(tmp[i + 0], tmp[i + 2], _MM_SHUFFLE(1, 0, 1, 0));
		tmpp[i + 1] = _mm256_shuffle_ps(tmp[i + 0], tmp[i + 2], _MM_SHUFFLE(3, 2, 3, 2));
	}
	for (int i = 0; i < 8; i += 4)
	{
		tmpp[i + 2] = _mm256_shuffle_ps(tmp[i + 1], tmp[i + 3], _MM_SHUFFLE(1, 0, 1, 0));
		tmpp[i + 3] = _mm256_shuffle_ps(tmp[i + 1], tmp[i + 3], _MM_SHUFFLE(3, 2, 3, 2));
	}
	for (int i = 0; i < 4; i++)
	{
		dst[i + 0] = _mm256_permute2f128_ps(tmpp[i], tmpp[i + 4], 0x20);
		dst[i + 4] = _mm256_permute2f128_ps(tmpp[i], tmpp[i + 4], 0x31);
	}
}