#pragma once

#ifdef __GNUC__
#include <x86intrin.h>
//#include <intrin.h>
#elif _MSC_VER
#include <intrin.h>
#pragma warning(disable: 4996)
#endif

#include <iostream>

inline void print_m256(const __m256 a)
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

inline void print_m256d(const __m256d a)
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

inline void print_m256i_i8(const __m256i a)
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

inline void print_m256i_u8(const __m256i a)
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
#elif _MSC_VER
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

inline void print_m256i_i16(const __m256i a)
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
#elif _MSC_VER
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

inline void print_m256i_u16(const __m256i a)
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
#elif _MSC_VER
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

inline void print_m256i_i32(const __m256i a)
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
#elif _MSC_VER
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

inline void print_m256i_u32(const __m256i a)
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
#elif _MSC_VER
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

inline void print_m256i_i64(const __m256i a)
{
#ifdef __GNUC__
	std::cout << ((long long*)&a)[0] << ", "
		<< ((long long*)&a)[1] << ", "
		<< ((long long*)&a)[2] << ", "
		<< ((long long*)&a)[3] << std::endl;
#elif _MSC_VER
	std::cout << a.m256i_i64[0] << ", "
		<< a.m256i_i64[1] << ", "
		<< a.m256i_i64[2] << ", "
		<< a.m256i_i64[3] << std::endl;
#endif
}

inline void print_m256i_u64(const __m256i a)
{
#ifdef __GNUC__
	std::cout << ((unsigned long long*) & a)[0] << ", "
		<< ((unsigned long long*) & a)[1] << ", "
		<< ((unsigned long long*) & a)[2] << ", "
		<< ((unsigned long long*) & a)[3] << std::endl;
#elif _MSC_VER
	std::cout << a.m256i_u64[0] << ", "
		<< a.m256i_u64[1] << ", "
		<< a.m256i_u64[2] << ", "
		<< a.m256i_u64[3] << std::endl;
#endif
}


inline void print_m128(const __m128 a)
{
#ifdef __GNUC__
	std::cout << ((float*)&a)[0] << ", "
		<< ((float*)&a)[1] << ", "
		<< ((float*)&a)[2] << ", "
		<< ((float*)&a)[3] << std::endl;
#elif _MSC_VER
	std::cout
		<< a.m128_f32[0] << ", "
		<< a.m128_f32[1] << ", "
		<< a.m128_f32[2] << ", "
		<< a.m128_f32[3] << std::endl;
#endif
}

inline void print_m128i_u8(const __m128i a)
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
#elif _MSC_VER
	std::cout
		<< (int)a.m128i_u8[0] << ", "
		<< (int)a.m128i_u8[1] << ", "
		<< (int)a.m128i_u8[2] << ", "
		<< (int)a.m128i_u8[3] << ", "
		<< (int)a.m128i_u8[4] << ", "
		<< (int)a.m128i_u8[5] << ", "
		<< (int)a.m128i_u8[6] << ", "
		<< (int)a.m128i_u8[7] << ", " //<< std::endl;
	//std::cout 
		<< (int)a.m128i_u8[8] << ", "
		<< (int)a.m128i_u8[9] << ", "
		<< (int)a.m128i_u8[10] << ", "
		<< (int)a.m128i_u8[11] << ", "
		<< (int)a.m128i_u8[12] << ", "
		<< (int)a.m128i_u8[13] << ", "
		<< (int)a.m128i_u8[14] << ", "
		<< (int)a.m128i_u8[15] << ", " << std::endl;
#endif
}

//unsigned charx16->floatx16
inline void _mm256_cvtepu8_psx2(__m128i src, __m256& dest0, __m256& dest1)
{
	dest0 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(src));
	dest1 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_shuffle_epi32(src, _MM_SHUFFLE(1, 0, 3, 2))));
}

//unsigned charx8 ->float
inline void _mm256_load_epu8cvtpsx2(const __m128i* P, __m256& dest0, __m256& dest1)
{
	__m128i src = _mm_load_si128(P);
	dest0 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(src));
	dest1 = _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_shuffle_epi32(src, _MM_SHUFFLE(1, 0, 3, 2))));
}

//floatx16->unsigned charx16
inline __m128i _mm256_cvtpsx2_epu8(const __m256 v0, const __m256 v1)
{
	return _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(_mm256_packus_epi16(_mm256_packs_epi32(_mm256_cvtps_epi32(v0), _mm256_cvtps_epi32(v1)), _mm256_setzero_si256()), _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7)));
}

//unsigned charx8 ->float
inline __m256 _mm256_load_epu8cvtps(const __m128i* P)
{
	return _mm256_cvtepi32_ps(_mm256_cvtepu8_epi32(_mm_loadl_epi64((__m128i*)P)));
}

//float->unsigned charx8
__m128i inline _mm256_cvtps_epu8(__m256 src)
{
	__m256i srci = _mm256_cvtps_epi32(src);//float-> int
	__m256i src16 = _mm256_packs_epi32(srci, _mm256_setzero_si256());//int->short:0Ç∆packÇ≈shortÇ…Ç∑ÇÈÅD
	__m256i src8 = _mm256_packus_epi16(src16, _mm256_setzero_si256());//short->uchar:0Ç∆pack"us"Ç≈unsigned charÇ…Ç∑ÇÈÅDuÇ™Ç»Ç¢Ç∆charÇ…Ç»ÇÈÅD
	__m256i src8perm = _mm256_permutevar8x32_epi32(src8, _mm256_setr_epi32(0, 4, 1, 5, 2, 6, 3, 7));//permutevarÇ≈ï¿Ç◊ë÷Ç¶
	return _mm256_castsi256_si128(src8perm);
}

inline void _mm256_load_cvtps_bgr2planar_ps(const float* ptr, __m256& a, __m256& b, __m256& c)
{
	__m256 bgr0 = _mm256_loadu_ps(ptr);
	__m256 bgr1 = _mm256_loadu_ps(ptr + 8);
	__m256 bgr2 = _mm256_loadu_ps(ptr + 16);

	__m256 s02_low = _mm256_permute2f128_ps(bgr0, bgr2, 0 + 2 * 16);
	__m256 s02_high = _mm256_permute2f128_ps(bgr0, bgr2, 1 + 3 * 16);

	__m256 b0 = _mm256_blend_ps(_mm256_blend_ps(s02_low, s02_high, 0x24), bgr1, 0x92);
	__m256 g0 = _mm256_blend_ps(_mm256_blend_ps(s02_high, s02_low, 0x92), bgr1, 0x24);
	__m256 r0 = _mm256_blend_ps(_mm256_blend_ps(bgr1, s02_low, 0x24), s02_high, 0x92);

	a = _mm256_shuffle_ps(b0, b0, 0x6c);
	b = _mm256_shuffle_ps(g0, g0, 0xb1);
	c = _mm256_shuffle_ps(r0, r0, 0xc6);
}

void inline _mm256_storeu_ps_color(void* dst, const __m256 rsrc, const __m256 gsrc, const __m256 bsrc)
{
	const int smask1 = _MM_SHUFFLE(1, 2, 3, 0);
	const int smask2 = _MM_SHUFFLE(2, 3, 0, 1);
	const int smask3 = _MM_SHUFFLE(3, 0, 1, 2);
	const int bmask1 = 0x44;
	const int bmask2 = 0x22;
	const int pmask1 = 0x20;
	const int pmask2 = 0x30;
	const int pmask3 = 0x31;
	const __m256 aa = _mm256_shuffle_ps(rsrc, rsrc, smask1);
	const __m256 bb = _mm256_shuffle_ps(gsrc, gsrc, smask2);
	const __m256 cc = _mm256_shuffle_ps(bsrc, bsrc, smask3);
	__m256 bval = _mm256_blend_ps(_mm256_blend_ps(aa, cc, bmask1), bb, bmask2);
	__m256 gval = _mm256_blend_ps(_mm256_blend_ps(cc, bb, bmask1), aa, bmask2);
	__m256 rval = _mm256_blend_ps(_mm256_blend_ps(bb, aa, bmask1), cc, bmask2);
	_mm256_storeu_ps((float*)dst + 0, _mm256_permute2f128_ps(bval, rval, pmask1));
	_mm256_storeu_ps((float*)dst + 8, _mm256_permute2f128_ps(gval, bval, pmask2));
	_mm256_storeu_ps((float*)dst + 16, _mm256_permute2f128_ps(rval, gval, pmask3));
}

void inline _mm256_stream_ps_color(void* dst, const __m256 rsrc, const __m256 gsrc, const __m256 bsrc)
{
	const int smask1 = _MM_SHUFFLE(1, 2, 3, 0);
	const int smask2 = _MM_SHUFFLE(2, 3, 0, 1);
	const int smask3 = _MM_SHUFFLE(3, 0, 1, 2);
	const int bmask1 = 0x44;
	const int bmask2 = 0x22;
	const int pmask1 = 0x20;
	const int pmask2 = 0x30;
	const int pmask3 = 0x31;
	const __m256 aa = _mm256_shuffle_ps(rsrc, rsrc, smask1);
	const __m256 bb = _mm256_shuffle_ps(gsrc, gsrc, smask2);
	const __m256 cc = _mm256_shuffle_ps(bsrc, bsrc, smask3);
	__m256 bval = _mm256_blend_ps(_mm256_blend_ps(aa, cc, bmask1), bb, bmask2);
	__m256 gval = _mm256_blend_ps(_mm256_blend_ps(cc, bb, bmask1), aa, bmask2);
	__m256 rval = _mm256_blend_ps(_mm256_blend_ps(bb, aa, bmask1), cc, bmask2);
	_mm256_stream_ps((float*)dst + 0, _mm256_permute2f128_ps(bval, rval, pmask1));
	_mm256_stream_ps((float*)dst + 8, _mm256_permute2f128_ps(gval, bval, pmask2));
	_mm256_stream_ps((float*)dst + 16, _mm256_permute2f128_ps(rval, gval, pmask3));
}

inline int get_simd_ceil(const int val, const int simdwidth)
{
	return (val % simdwidth == 0) ? val : (val / simdwidth + 1) * simdwidth;
}

inline int get_simd_floor(const int val, const int simdwidth)
{
	return (val / simdwidth) * simdwidth;
}