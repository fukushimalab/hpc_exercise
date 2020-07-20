#pragma once

#include "mat_util.h"
#include "simd_util.h"
#include <iostream>
#include <algorithm>
#include <cfloat>

#include <omp.h>

template<int size>
void loofline_test_cpp(const int iteration, const int num_thread = -1)
{
	//x+1: (1*x+1)
	//x*x+x+1: x*(x+1)+1
	//x*x*x+x*x+x+1: x*(x*(x+1))+1
	const int thread_max = (num_thread == -1) ? omp_get_num_procs() : num_thread;
	omp_set_num_threads(thread_max);

	std::cout << "loofline test: " << thread_max << "threads" << std::endl;
	const int loop = iteration;

	CalcTime t;

	float* x = (float*)_mm_malloc(sizeof(float) * size * thread_max, 32);
	float* y = (float*)_mm_malloc(sizeof(float) * size * thread_max, 32);

	float* ptr = x;
	const float rand_max = 1.f;
	const float rand_min = 0.f;
	const float v = (float)(rand_max - rand_min) / static_cast<float>(RAND_MAX);
	for (int i = 0; i < size * thread_max; i++)
	{
		*ptr++ = rand_min + (rand() * v);
	}

	printf("size %d KBYTE, iteration %d\n", size * (int)(sizeof(float)) / 1024, iteration);
	printf("order, GFLOPS, FLOPS/BYTE\n");

	int n = 0;
	{
		//---------------------------------------------
		n = 1;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				py[i] = px[i] + 1.f;//1
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 2;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				py[i] = px[i] * temp + 1.f;//2
			}

		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 3;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				py[i] = px[i] * temp + 1.f;//3
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 4;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				py[i] = px[i] * temp + 1.f;//4
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 5;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				py[i] = px[i] * temp + 1.f;//5
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 6;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				py[i] = px[i] * temp + 1.f;//6
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 7;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				py[i] = px[i] * temp + 1.f;//7
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 8;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				py[i] = px[i] * temp + 1.f;//8
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 9;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				py[i] = px[i] * temp + 1.f;//9
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 10;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				py[i] = px[i] * temp + 1.f;//10
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 11;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				py[i] = px[i] * temp + 1.f;//11
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 12;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				temp = px[i] * temp + 1.f;//11
				py[i] = px[i] * temp + 1.f;//12
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 13;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				temp = px[i] * temp + 1.f;//11
				temp = px[i] * temp + 1.f;//12
				py[i] = px[i] * temp + 1.f;//13
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 14;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				temp = px[i] * temp + 1.f;//11
				temp = px[i] * temp + 1.f;//12
				temp = px[i] * temp + 1.f;//13
				py[i] = px[i] * temp + 1.f;//14
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 15;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				temp = px[i] * temp + 1.f;//11
				temp = px[i] * temp + 1.f;//12
				temp = px[i] * temp + 1.f;//13
				temp = px[i] * temp + 1.f;//14
				py[i] = px[i] * temp + 1.f;//15
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 16;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				temp = px[i] * temp + 1.f;//11
				temp = px[i] * temp + 1.f;//12
				temp = px[i] * temp + 1.f;//13
				temp = px[i] * temp + 1.f;//14
				temp = px[i] * temp + 1.f;//15
				py[i] = px[i] * temp + 1.f;//16
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 17;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				temp = px[i] * temp + 1.f;//11
				temp = px[i] * temp + 1.f;//12
				temp = px[i] * temp + 1.f;//13
				temp = px[i] * temp + 1.f;//14
				temp = px[i] * temp + 1.f;//15
				temp = px[i] * temp + 1.f;//16
				py[i] = px[i] * temp + 1.f;//17
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 18;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				temp = px[i] * temp + 1.f;//11
				temp = px[i] * temp + 1.f;//12
				temp = px[i] * temp + 1.f;//13
				temp = px[i] * temp + 1.f;//14
				temp = px[i] * temp + 1.f;//15
				temp = px[i] * temp + 1.f;//16
				temp = px[i] * temp + 1.f;//17
				py[i] = px[i] * temp + 1.f;//18
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 19;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				temp = px[i] * temp + 1.f;//11
				temp = px[i] * temp + 1.f;//12
				temp = px[i] * temp + 1.f;//13
				temp = px[i] * temp + 1.f;//14
				temp = px[i] * temp + 1.f;//15
				temp = px[i] * temp + 1.f;//16
				temp = px[i] * temp + 1.f;//17
				temp = px[i] * temp + 1.f;//18
				py[i] = px[i] * temp + 1.f;//19
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 20;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			for (int i = 0; i < size; i++)
			{
				float temp = px[i] + 1.f;//1
				temp = px[i] * temp + 1.f;//2
				temp = px[i] * temp + 1.f;//3
				temp = px[i] * temp + 1.f;//4
				temp = px[i] * temp + 1.f;//5
				temp = px[i] * temp + 1.f;//6
				temp = px[i] * temp + 1.f;//7
				temp = px[i] * temp + 1.f;//8
				temp = px[i] * temp + 1.f;//9
				temp = px[i] * temp + 1.f;//10
				temp = px[i] * temp + 1.f;//11
				temp = px[i] * temp + 1.f;//12
				temp = px[i] * temp + 1.f;//13
				temp = px[i] * temp + 1.f;//14
				temp = px[i] * temp + 1.f;//15
				temp = px[i] * temp + 1.f;//16
				temp = px[i] * temp + 1.f;//17
				temp = px[i] * temp + 1.f;//18
				temp = px[i] * temp + 1.f;//19
				py[i] = px[i] * temp + 1.f;//20
			}
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));
	}

	_mm_free(x);
	_mm_free(y);
}

#define USE_DAZ
#define USE_FTZ
template<int size>
void loofline_test(const int iteration, const int num_thread = -1)
{

#ifdef USE_DAZ
#include <pmmintrin.h>
	_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
	
#ifdef USE_FTZ
#include <xmmintrin.h>
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif
	//show_mxcsr(true, true, false); //check status

	//x+1: (1*x+1)
	//x*x+x+1: x*(x+1)+1
	//x*x*x+x*x+x+1: x*(x*(x+1))+1

	const int thread_max = (num_thread == -1) ? omp_get_num_procs() : num_thread;
	omp_set_num_threads(thread_max);

	std::cout << "loofline test: " << thread_max << "threads" << std::endl;
	const int loop = iteration;

	CalcTime t;

	float* x = (float*)_mm_malloc(sizeof(float) * size * thread_max, 32);
	float* y = (float*)_mm_malloc(sizeof(float) * size * thread_max, 32);

	float* ptr = x;
	const float rand_max = 1.f;
	const float rand_min = 0.f;
	const float v = (float)(rand_max - rand_min) / static_cast<float>(RAND_MAX);
	for (int i = 0; i < size * thread_max; i++)
	{
		*ptr++ = std::max(FLT_MIN, rand_min + (rand() * v));
	}

	int simdsize16 = size / 16;

	printf("size %d KBYTE, iteration %d\n", size * (int)(sizeof(float)) / 1024, iteration);
	printf("order, FLOPS/BYTE, GFLOPS\n");

	int n = 0;
	{
		//---------------------------------------------
		n = 1;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 2;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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

		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 3;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 4;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 5;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 6;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 7;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 8;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 9;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 10;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 11;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 12;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 13;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 14;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 15;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 16;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 17;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 18;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 19;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));

		//---------------------------------------------
		n = 20;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			const int v = omp_get_thread_num();
			float* px = x + v * size;
			float* py = y + v * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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
		}
		t.end();
		printf("%02d, %f, %f\n", n, n * 2.0 / (2 * 4), n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000));
	}

	_mm_free(x);
	_mm_free(y);
}

template<int size>
void timer_test(const int iteration, const int num_thread = -1)
{
	const int thread_max = (num_thread == -1) ? omp_get_max_threads() : num_thread;
	omp_set_num_threads(thread_max);

	std::cout << "timer test: " << thread_max << "threads" << std::endl;
	const int loop = iteration;

	CalcTime t;

	float* x = (float*)_mm_malloc(sizeof(float) * size * thread_max, 32);
	float* y = (float*)_mm_malloc(sizeof(float) * size * thread_max, 32);

	float* ptr = x;
	const float rand_max = 1.f;
	const float rand_min = 0.f;
	const float v = (float)(rand_max - rand_min) / static_cast<float>(RAND_MAX);
	for (int i = 0; i < size * thread_max; i++)
	{
		*ptr++ = rand_min + (rand() * v);
	}

	int simdsize16 = size / 16;

	printf("size %d KBYTE, iteration %d\n", size * (int)(sizeof(float)) / 1024, iteration);
	printf("order, GFLOPS, FLOPS/BYTE\n");

	int n = 0;
	{
		n = 20;
		t.start();
#pragma omp parallel for
		for (int j = 0; j < loop; j++)
		{
			float* px = x + omp_get_thread_num() * size;
			float* py = y + omp_get_thread_num() * size;
			const __m256 mones = _mm256_set1_ps(1.f);
			for (int i = 0; i < simdsize16; i++)
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

		}
		t.end();
		printf("%f sec \n", t.getLastTime() / 1000);
		printf("%02d, %f, %f\n", n, n * 2.0 * size / (t.getLastTime() / loop * 0.001) / (1000 * 1000 * 1000), n * 2.0 / (2 * 4));
	}

	_mm_free(x);
	_mm_free(y);
}