#pragma once
#ifdef __GNUC__
#include <x86intrin.h>
#elif _MSC_VER
#include <intrin.h>
#endif
#include <iostream>

void print_m256(const __m256 a)
{
#ifdef __GNUC__
	std::cout << ((float*)&a)[0] << ", "
		<< ((float*)&a)[1] << ", "
		<< ((float*)&a)[2] << ", "
		<< ((float*)&a)[3] << ", "
		<< ((float*)&a)[4] << ", "
		<< ((float*)&a)[5] << ", "
		<< ((float*)&a)[6] << ", "
		<< ((float*)&a)[7] << std::endl;
#elif _MSC_VER
	std::cout << a.m256_f32[0] << ", "
		<< a.m256_f32[1] << ", "
		<< a.m256_f32[2] << ", "
		<< a.m256_f32[3] << ", "
		<< a.m256_f32[4] << ", "
		<< a.m256_f32[5] << ", "
		<< a.m256_f32[6] << ", "
		<< a.m256_f32[7] << std::endl;
#endif
}

void print_m256d(const __m256d a)
{
#ifdef __GNUC__
	std::cout << ((double*)&a)[0] << ", "
		<< ((double*)&a)[1] << ", "
		<< ((double*)&a)[2] << ", "
		<< ((double*)&a)[3] << std::endl;
#elif _MSC_VER
	std::cout << a.m256d_f64[0] << ", "
		<< a.m256d_f64[1] << ", "
		<< a.m256d_f64[2] << ", "
		<< a.m256d_f64[3] << std::endl;
#endif
}

void print_m256i_i8(const __m256i a)
{
#ifdef __GNUC__
	std::cout << (int)((char*)&a)[0] << ", "
		<< (int)((char*)&a)[1] << ", "
		<< (int)((char*)&a)[2] << ", "
		<< (int)((char*)&a)[3] << ", "
		<< (int)((char*)&a)[4] << ", "
		<< (int)((char*)&a)[5] << ", "
		<< (int)((char*)&a)[6] << ", "
		<< (int)((char*)&a)[7] << ", " << std::endl;
	std::cout << (int)((char*)&a)[8] << ", "
		<< (int)((char*)&a)[9] << ", "
		<< (int)((char*)&a)[10] << ", "
		<< (int)((char*)&a)[11] << ", "
		<< (int)((char*)&a)[12] << ", "
		<< (int)((char*)&a)[13] << ", "
		<< (int)((char*)&a)[14] << ", "
		<< (int)((char*)&a)[15] << ", " << std::endl;
	std::cout << (int)((char*)&a)[16] << ", "
		<< (int)((char*)&a)[17] << ", "
		<< (int)((char*)&a)[18] << ", "
		<< (int)((char*)&a)[19] << ", "
		<< (int)((char*)&a)[20] << ", "
		<< (int)((char*)&a)[21] << ", "
		<< (int)((char*)&a)[22] << ", "
		<< (int)((char*)&a)[23] << ", " << std::endl;
	std::cout << (int)((char*)&a)[24] << ", "
		<< (int)((char*)&a)[25] << ", "
		<< (int)((char*)&a)[26] << ", "
		<< (int)((char*)&a)[27] << ", "
		<< (int)((char*)&a)[28] << ", "
		<< (int)((char*)&a)[29] << ", "
		<< (int)((char*)&a)[30] << ", "
		<< (int)((char*)&a)[31] << ", " << std::endl;
#elif _MSC_VER
	std::cout << (int)a.m256i_i8[0] << ", "
		<< (int)a.m256i_i8[1] << ", "
		<< (int)a.m256i_i8[2] << ", "
		<< (int)a.m256i_i8[3] << ", "
		<< (int)a.m256i_i8[4] << ", "
		<< (int)a.m256i_i8[5] << ", "
		<< (int)a.m256i_i8[6] << ", "
		<< (int)a.m256i_i8[7] << ", " << std::endl;
	std::cout << (int)a.m256i_i8[8] << ", "
		<< (int)a.m256i_i8[9] << ", "
		<< (int)a.m256i_i8[10] << ", "
		<< (int)a.m256i_i8[11] << ", "
		<< (int)a.m256i_i8[12] << ", "
		<< (int)a.m256i_i8[13] << ", "
		<< (int)a.m256i_i8[14] << ", "
		<< (int)a.m256i_i8[15] << ", " << std::endl;
	std::cout << (int)a.m256i_i8[16] << ", "
		<< (int)a.m256i_i8[17] << ", "
		<< (int)a.m256i_i8[18] << ", "
		<< (int)a.m256i_i8[19] << ", "
		<< (int)a.m256i_i8[20] << ", "
		<< (int)a.m256i_i8[21] << ", "
		<< (int)a.m256i_i8[22] << ", "
		<< (int)a.m256i_i8[23] << ", " << std::endl;
	std::cout << (int)a.m256i_i8[24] << ", "
		<< (int)a.m256i_i8[25] << ", "
		<< (int)a.m256i_i8[26] << ", "
		<< (int)a.m256i_i8[27] << ", "
		<< (int)a.m256i_i8[28] << ", "
		<< (int)a.m256i_i8[29] << ", "
		<< (int)a.m256i_i8[30] << ", "
		<< (int)a.m256i_i8[31] << ", " << std::endl;
#endif
}

void print_m256i_u8(const __m256i a)
{
#ifdef __GNUC__
	std::cout << (int)((unsigned char*)&a)[0] << ", "
		<< (int)((unsigned char*)&a)[1] << ", "
		<< (int)((unsigned char*)&a)[2] << ", "
		<< (int)((unsigned char*)&a)[3] << ", "
		<< (int)((unsigned char*)&a)[4] << ", "
		<< (int)((unsigned char*)&a)[5] << ", "
		<< (int)((unsigned char*)&a)[6] << ", "
		<< (int)((unsigned char*)&a)[7] << ", " << std::endl;
	std::cout << (int)((unsigned char*)&a)[8] << ", "
		<< (int)((unsigned char*)&a)[9] << ", "
		<< (int)((unsigned char*)&a)[10] << ", "
		<< (int)((unsigned char*)&a)[11] << ", "
		<< (int)((unsigned char*)&a)[12] << ", "
		<< (int)((unsigned char*)&a)[13] << ", "
		<< (int)((unsigned char*)&a)[14] << ", "
		<< (int)((unsigned char*)&a)[15] << ", " << std::endl;
	std::cout << (int)((unsigned char*)&a)[16] << ", "
		<< (int)((unsigned char*)&a)[17] << ", "
		<< (int)((unsigned char*)&a)[18] << ", "
		<< (int)((unsigned char*)&a)[19] << ", "
		<< (int)((unsigned char*)&a)[20] << ", "
		<< (int)((unsigned char*)&a)[21] << ", "
		<< (int)((unsigned char*)&a)[22] << ", "
		<< (int)((unsigned char*)&a)[23] << ", " << std::endl;
	std::cout << (int)((unsigned char*)&a)[24] << ", "
		<< (int)((unsigned char*)&a)[25] << ", "
		<< (int)((unsigned char*)&a)[26] << ", "
		<< (int)((unsigned char*)&a)[27] << ", "
		<< (int)((unsigned char*)&a)[28] << ", "
		<< (int)((unsigned char*)&a)[29] << ", "
		<< (int)((unsigned char*)&a)[30] << ", "
		<< (int)((unsigned char*)&a)[31] << ", " << std::endl;
#elif MSC_VER
	std::cout << (int)a.m256i_u8[0] << ", "
		<< (int)a.m256i_u8[1] << ", "
		<< (int)a.m256i_u8[2] << ", "
		<< (int)a.m256i_u8[3] << ", "
		<< (int)a.m256i_u8[4] << ", "
		<< (int)a.m256i_u8[5] << ", "
		<< (int)a.m256i_u8[6] << ", "
		<< (int)a.m256i_u8[7] << ", " << std::endl;
	std::cout << (int)a.m256i_u8[8] << ", "
		<< (int)a.m256i_u8[9] << ", "
		<< (int)a.m256i_u8[10] << ", "
		<< (int)a.m256i_u8[11] << ", "
		<< (int)a.m256i_u8[12] << ", "
		<< (int)a.m256i_u8[13] << ", "
		<< (int)a.m256i_u8[14] << ", "
		<< (int)a.m256i_u8[15] << ", " << std::endl;
	std::cout << (int)a.m256i_u8[16] << ", "
		<< (int)a.m256i_u8[17] << ", "
		<< (int)a.m256i_u8[18] << ", "
		<< (int)a.m256i_u8[19] << ", "
		<< (int)a.m256i_u8[20] << ", "
		<< (int)a.m256i_u8[21] << ", "
		<< (int)a.m256i_u8[22] << ", "
		<< (int)a.m256i_u8[23] << ", " << std::endl;
	std::cout << (int)a.m256i_u8[24] << ", "
		<< (int)a.m256i_u8[25] << ", "
		<< (int)a.m256i_u8[26] << ", "
		<< (int)a.m256i_u8[27] << ", "
		<< (int)a.m256i_u8[28] << ", "
		<< (int)a.m256i_u8[29] << ", "
		<< (int)a.m256i_u8[30] << ", "
		<< (int)a.m256i_u8[31] << ", " << std::endl;
#endif
}

void print_m256i_i16(const __m256i a)
{
#ifdef __GNUC__
	std::cout << ((short*)&a)[0] << ", "
		<< ((short*)&a)[1] << ", "
		<< ((short*)&a)[2] << ", "
		<< ((short*)&a)[3] << ", "
		<< ((short*)&a)[4] << ", "
		<< ((short*)&a)[5] << ", "
		<< ((short*)&a)[6] << ", "
		<< ((short*)&a)[7] << ", " << std::endl;
	std::cout << (int)((short*)&a)[8] << ", "
		<< ((short*)&a)[9] << ", "
		<< ((short*)&a)[10] << ", "
		<< ((short*)&a)[11] << ", "
		<< ((short*)&a)[12] << ", "
		<< ((short*)&a)[13] << ", "
		<< ((short*)&a)[14] << ", "
		<< ((short*)&a)[15] << ", " << std::endl;
#elif MSC_VER
	std::cout << a.m256i_i16[0] << ", "
		<< a.m256i_i16[1] << ", "
		<< a.m256i_i16[2] << ", "
		<< a.m256i_i16[3] << ", "
		<< a.m256i_i16[4] << ", "
		<< a.m256i_i16[5] << ", "
		<< a.m256i_i16[6] << ", "
		<< a.m256i_i16[7] << ", " << std::endl;
	std::cout << a.m256i_i16[8] << ", "
		<< a.m256i_i16[9] << ", "
		<< a.m256i_i16[10] << ", "
		<< a.m256i_i16[11] << ", "
		<< a.m256i_i16[12] << ", "
		<< a.m256i_i16[13] << ", "
		<< a.m256i_i16[14] << ", "
		<< a.m256i_i16[15] << ", " << std::endl;
#endif
}

void print_m256i_u16(const __m256i a)
{
#ifdef __GNUC__
	std::cout << ((unsigned short*)&a)[0] << ", "
		<< ((unsigned short*)&a)[1] << ", "
		<< ((unsigned short*)&a)[2] << ", "
		<< ((unsigned short*)&a)[3] << ", "
		<< ((unsigned short*)&a)[4] << ", "
		<< ((unsigned short*)&a)[5] << ", "
		<< ((unsigned short*)&a)[6] << ", "
		<< ((unsigned short*)&a)[7] << ", " << std::endl;
	std::cout << (int)((unsigned short*)&a)[8] << ", "
		<< ((unsigned short*)&a)[9] << ", "
		<< ((unsigned short*)&a)[10] << ", "
		<< ((unsigned short*)&a)[11] << ", "
		<< ((unsigned short*)&a)[12] << ", "
		<< ((unsigned short*)&a)[13] << ", "
		<< ((unsigned short*)&a)[14] << ", "
		<< ((unsigned short*)&a)[15] << ", " << std::endl;
#elif MSC_VER
	std::cout << a.m256i_u16[0] << ", "
		<< a.m256i_u16[1] << ", "
		<< a.m256i_u16[2] << ", "
		<< a.m256i_u16[3] << ", "
		<< a.m256i_u16[4] << ", "
		<< a.m256i_u16[5] << ", "
		<< a.m256i_u16[6] << ", "
		<< a.m256i_u16[7] << ", " << std::endl;
	std::cout << a.m256i_u16[8] << ", "
		<< a.m256i_u16[9] << ", "
		<< a.m256i_u16[10] << ", "
		<< a.m256i_u16[11] << ", "
		<< a.m256i_u16[12] << ", "
		<< a.m256i_u16[13] << ", "
		<< a.m256i_u16[14] << ", "
		<< a.m256i_u16[15] << ", " << std::endl;
#endif
}

void print_m256i_i32(const __m256i a)
{
#ifdef __GNUC__
	std::cout << ((int*)&a)[0] << ", "
		<< ((int*)&a)[1] << ", "
		<< ((int*)&a)[2] << ", "
		<< ((int*)&a)[3] << ", "
		<< ((int*)&a)[4] << ", "
		<< ((int*)&a)[5] << ", "
		<< ((int*)&a)[6] << ", "
		<< ((int*)&a)[7] << ", " << std::endl;
#elif MSC_VER
	std::cout << a.m256i_i32[0] << ", "
		<< a.m256i_i32[1] << ", "
		<< a.m256i_i32[2] << ", "
		<< a.m256i_i32[3] << ", "
		<< a.m256i_i32[4] << ", "
		<< a.m256i_i32[5] << ", "
		<< a.m256i_i32[6] << ", "
		<< a.m256i_i32[7] << ", " << std::endl;
#endif
}

void print_m256i_u32(const __m256i a)
{
#ifdef __GNUC__
	std::cout << ((unsigned int*)&a)[0] << ", "
		<< ((unsigned int*)&a)[1] << ", "
		<< ((unsigned int*)&a)[2] << ", "
		<< ((unsigned int*)&a)[3] << ", "
		<< ((unsigned int*)&a)[4] << ", "
		<< ((unsigned int*)&a)[5] << ", "
		<< ((unsigned int*)&a)[6] << ", "
		<< ((unsigned int*)&a)[7] << ", " << std::endl;
#elif MSC_VER
	std::cout << a.m256i_u32[0] << ", "
		<< a.m256i_u32[1] << ", "
		<< a.m256i_u32[2] << ", "
		<< a.m256i_u32[3] << ", "
		<< a.m256i_u32[4] << ", "
		<< a.m256i_u32[5] << ", "
		<< a.m256i_u32[6] << ", "
		<< a.m256i_u32[7] << ", " << std::endl;
#endif
}

void print_m256i_i64(const __m256i a)
{
#ifdef __GNUC__
	std::cout << ((long long*)&a)[0] << ", "
		<< ((long long*)&a)[1] << ", "
		<< ((long long*)&a)[2] << ", "
		<< ((long long*)&a)[3] << std::endl;
#elif MSC_VER
	std::cout << a.m256i_i64[0] << ", "
		<< a.m256i_i64[1] << ", "
		<< a.m256i_i64[2] << ", "
		<< a.m256i_i64[3] << std::endl;
#endif
}

void print_m256i_u64(const __m256i a)
{
#ifdef __GNUC__
	std::cout << ((unsigned long long*)&a)[0] << ", "
		<< ((unsigned long long*)&a)[1] << ", "
		<< ((unsigned long long*)&a)[2] << ", "
		<< ((unsigned long long*)&a)[3] << std::endl;
#elif MSC_VER
	std::cout << a.m256i_u64[0] << ", "
		<< a.m256i_u64[1] << ", "
		<< a.m256i_u64[2] << ", "
		<< a.m256i_u64[3] << std::endl;
#endif
}

