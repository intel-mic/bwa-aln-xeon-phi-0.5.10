#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifdef __MIC__
#include <immintrin.h>
#include <zmmintrin.h>
#endif



typedef unsigned int uint32_t;
typedef unsigned char ubyte_t;
typedef uint32_t bwtint_t;

#define count 100000

#define OCC_INTERVAL 0x80	// 2^7 = 128

#define CHR_PTR_OPT			// get byte by char pointer instead of >>& 
#define RM_MEMCPY_OPT		// remove memcpy
#define SIMP_OPT			// simplify loops


typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	// occurance array, separated to two parts
	uint32_t cnt_table[256];
	// suffix array
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;



void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {	// 8 bits: 2^8 = 256, 1 0000 0000, all possibility of 1 byte
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)	// ACTG, 4bp, 2 bits per bp, 1 bytes = 8 bits = (4 2bits)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);	// (0 ~ 3) << 0/8/16/24
		bwt->cnt_table[i] = x;
		// index: all possibility of 1 byte
		// value: (0:7): number of (bp0), (9:15): number of (bp1), (16:23): number of (bp2), (24:31), number of (bp3)
	}
}



#define bwt_occ_intv(b, k) ((b)->bwt + (k)/OCC_INTERVAL*12)

#define __occ_aux4(bwt, b)											\
	((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
	 + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

// 4 bytes: (byte0 4bp counts) + (byte1 4bp counts) + (byte2 4bp count) + (byte3 4bp counts)
// input: int32, output: 4 8bits(4 bp counts)
// 只有一个位置k的计算
inline void bwt_occ4_opt(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
	bwtint_t l, j, x;
	uint32_t *p;
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, 4 * sizeof(bwtint_t));
		return;
	}
	if (k >= bwt->primary) --k; // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	//p = bwt->bwt + k/OCC_INTERVAL*12;
#ifndef RM_MEMCPY_OPT	
	memcpy(cnt, p, 16);
#else
	uint32_t *p0 = p;
	_mm_prefetch((char*)p0, _MM_HINT_T0);
#endif

	p += 4;
	
#ifdef SIMP_OPT

	int ii;
	j = k & 0xfffffff0;				// zero lower 4 bits，1111 0000, k >> 4 << 4
	int start = k & 0xffffff80;		// zero lower 7 bits, 1000 0000, k / OCC_INTERVAL * OCC_INTERVAL
	int end   = j;
	int len	  = (end - start)>>4;	// /16 = 2^4

	x = 0;
	//printf("--------start, len = %d, p = %x\n", len, p);

#if 1

#pragma novector
	for (ii = 0; ii < len; ii ++){
		unsigned char* p_chr = (unsigned char*)&p[ii];
		x += bwt->cnt_table[p_chr[0]] + bwt->cnt_table[p_chr[1]] + bwt->cnt_table[p_chr[2]] + bwt->cnt_table[p_chr[3]];
	}

#else

	unsigned char* p_ptr = (unsigned char*)p;
#pragma simd
	for (ii = 0; ii < len*4; ii ++){
		//printf("p = %x\n", p);
		x += bwt->cnt_table[p_ptr[ii]];// + bwt->cnt_table[p_ptr[ii+1]] + bwt->cnt_table[p_ptr[ii+2]] + bwt->cnt_table[p_ptr[ii+3]];
	}
	
#endif

	p += len;
	//printf("--------end, p = %x\n", p);	
	
#else	// not SIMP_OPT
	
	j = k >> 4 << 4;
	
#pragma novector
	for (l = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; l < j; l += 16, ++p){
#ifndef CHR_PTR_OPT
		x += __occ_aux4(bwt, *p);	// opt-report: gather was generated for the variable bwt
#else
		unsigned char* p_chr = (unsigned char*)p;
		x += bwt->cnt_table[p_chr[0]] + bwt->cnt_table[p_chr[1]] + bwt->cnt_table[p_chr[2]] + bwt->cnt_table[p_chr[3]];
#endif
	}
	
#endif 	// SIMP_OPT

	x += __occ_aux4(bwt, *p & ~((1U<<((~k&15)<<1)) - 1)) - (~k&15);
	
#ifndef RM_MEMCPY_OPT	
	cnt[0] += x&0xff;
	cnt[1] += x>>8&0xff;
	cnt[2] += x>>16&0xff;
	cnt[3] += x>>24;
#else
	unsigned char *x_chr = (unsigned char *)&x;
	cnt[0] = p0[0] + x_chr[0]; cnt[1] = p0[1] + x_chr[1]; cnt[2] = p0[2] + x_chr[2]; cnt[3] = p0[3] + x_chr[3];
#endif

	//cnt[0] += x&0xff+(bwt_occ_intv(bwt,k))[0];
	//cnt[1] += x>>8&0xff + (bwt_occ_intv(bwt,k))[1];
	//cnt[2] += x>>16&0xff + (bwt_occ_intv(bwt,k))[2];
	//cnt[3] += x>>24 + (bwt_occ_intv(bwt,k))[2];
}

// #define OCC_INTERVAL 0x80: 128, 2^7
// typedef uint32_t bwtint_t;
// an analogy to bwt_occ4_opt() but more efficient, requiring k <= l
// 既有位置k又有位置l的计算
inline void bwt_2occ4_opt(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
	bwtint_t _k, _l;
	if (k == l) {					
		bwt_occ4_opt(bwt, k, cntk);				// opt-report: inlined bwt_occ4，bwt_occ4，k==l，当前只有一个位置
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
		return;
	}
	_k = (k >= bwt->primary)? k-1 : k;			// primary, 初始位置
	_l = (l >= bwt->primary)? l-1 : l;			
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		bwt_occ4_opt(bwt, k, cntk);				// 只有一个位置
		bwt_occ4_opt(bwt, l, cntl);
	} else {
		bwtint_t i, j, x, y;
		uint32_t *p;
		int cl[4];
		if (k >= bwt->primary) {--k; /*fprintf(stderr, "------- --k\n");*/	}// because $ is not in bwt 
		if (l >= bwt->primary) --l;
		cl[0] = cl[1] = cl[2] = cl[3] = 0;
		
		p = bwt_occ_intv(bwt, k);			// bwt: ref中，找到k的位置

#ifndef RM_MEMCPY_OPT
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
#else
		uint32_t *p0 = p;
		_mm_prefetch((char*)p0, _MM_HINT_T0);
#endif

		p += 4;
		// prepare cntk[]

#ifdef SIMP_OPT	
	
	int ii;
	j = k & 0xfffffff0;				// zero lower 4 bits，1111 0000, k >> 4 << 4
	int start = k & 0xffffff80;		// zero lower 7 bits, 1000 0000, k / OCC_INTERVAL * OCC_INTERVAL
	int end   = j;
	int len	  = (end - start)>>4;	// /16 = 2^4
	x = 0;
#pragma novector
	for (ii = 0; ii < len; ii ++){
		unsigned char* p_chr = (unsigned char*)&p[ii];
		x += bwt->cnt_table[p_chr[0]] + bwt->cnt_table[p_chr[1]] + bwt->cnt_table[p_chr[2]] + bwt->cnt_table[p_chr[3]];
	}
	p += len;
	i = start + len<<4;				// len * 16

#else	// not SIMP_OPT

		j = k >> 4 << 4;

#pragma novector
		// k, 128位倍数区间段内的数字读出来, j把k值后面几位去掉，i += 16每个取16个碱基
		for (i = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; i < j; i += 16, ++p){	
#ifndef CHR_PTR_OPT
			x += __occ_aux4(bwt, *p);			// gather，多少个碱基
#else
			unsigned char* p_chr = (unsigned char*)p;
			x += bwt->cnt_table[p_chr[0]] + bwt->cnt_table[p_chr[1]] + bwt->cnt_table[p_chr[2]] + bwt->cnt_table[p_chr[3]];			
#endif
		}
		
#endif 	// SIMP_OPT

		y = x;															// y用来做l的值
		x += __occ_aux4(bwt, *p & ~((1U<<((~k&15)<<1)) - 1)) - (~k&15);	// 利用前面结果+结果
		// calculate cntl[] and finalize cntk[]

#if 0//def SIMP_OPT	

	j = l & 0xfffffff0;				// zero lower 4 bits，1111 0000, k >> 4 << 4
	start = i;						// zero lower 7 bits, 1000 0000, k / OCC_INTERVAL * OCC_INTERVAL
	end   = j;
	len	  = (end - start)>>4;		// /16 = 2^4

#pragma novector
	for (ii = start; ii < end; ii += 16, ++p) {							// y累加计算l
	//for (ii = 0; ii < len; ii ++) {
		unsigned char* p_chr = (unsigned char*)p;
		//unsigned char* p_chr = (unsigned char*)&p[ii];
		y += bwt->cnt_table[p_chr[0]] + bwt->cnt_table[p_chr[1]] + bwt->cnt_table[p_chr[2]] + bwt->cnt_table[p_chr[3]];		
	}
	//p += len;

#else	// not SIMP_OPT		
		
		j = l >> 4 << 4;
		
#pragma novector
		for (; i < j; i += 16, ++p) {									// y累加计算l
#ifndef CHR_PTR_OPT
			y += __occ_aux4(bwt, *p);	// opt-report: prefetch
#else
			unsigned char* p_chr = (unsigned char*)p;
			y += bwt->cnt_table[p_chr[0]] + bwt->cnt_table[p_chr[1]] + bwt->cnt_table[p_chr[2]] + bwt->cnt_table[p_chr[3]];			
#endif			
		}
		
#endif 	// SIMP_OPT		
		
		y += __occ_aux4(bwt, *p & ~((1U<<((~l&15)<<1)) - 1)) - (~l&15);
		
#ifndef RM_MEMCPY_OPT
		memcpy(cntl, cntk, 16);
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;	// 0~k有多少个碱基
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;	// 0~l有多少个碱基
#else
		unsigned char *x_chr, *y_chr;
		x_chr = (unsigned char *)&x;
		y_chr = (unsigned char *)&y;
		cntk[0] = p0[0] + x_chr[0]; cntk[1] = p0[1] + x_chr[1]; cntk[2] = p0[2] + x_chr[2]; cntk[3] = p0[3] + x_chr[3];
		cntl[0] = p0[0] + y_chr[0]; cntl[1] = p0[1] + y_chr[1]; cntl[2] = p0[2] + y_chr[2]; cntl[3] = p0[3] + y_chr[3];
#endif
	}
}


// Ori
inline void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
	bwtint_t l, j, x;
	uint32_t *p;
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, 4 * sizeof(bwtint_t));
		return;
	}
	if (k >= bwt->primary) --k; // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	//p = bwt->bwt + k/OCC_INTERVAL*12;
	memcpy(cnt, p, 16);
	p += 4;
	j = k >> 4 << 4;
#pragma novector
	for (l = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; l < j; l += 16, ++p)
		x += __occ_aux4(bwt, *p);
	x += __occ_aux4(bwt, *p & ~((1U<<((~k&15)<<1)) - 1)) - (~k&15);
	cnt[0] += x&0xff;
	cnt[1] += x>>8&0xff;
	cnt[2] += x>>16&0xff;
	cnt[3] += x>>24;
	//cnt[0] += x&0xff+(bwt_occ_intv(bwt,k))[0];
	//cnt[1] += x>>8&0xff + (bwt_occ_intv(bwt,k))[1];
	//cnt[2] += x>>16&0xff + (bwt_occ_intv(bwt,k))[2];
	//cnt[3] += x>>24 + (bwt_occ_intv(bwt,k))[2];
}

// #define OCC_INTERVAL 0x80: 128, 2^7
// typedef uint32_t bwtint_t;
// an analogy to bwt_occ4() but more efficient, requiring k <= l
// 既有位置k又有位置l的计算
inline void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
	bwtint_t _k, _l;
	if (k == l) {
		bwt_occ4(bwt, k, cntk);
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
		return;
	}
	_k = (k >= bwt->primary)? k-1 : k;
	_l = (l >= bwt->primary)? l-1 : l;
	if (_l/OCC_INTERVAL != _k/OCC_INTERVAL || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		bwtint_t i, j, x, y;
		uint32_t *p;
		int cl[4];
		if (k >= bwt->primary) --k; 	// because $ is not in bwt
		if (l >= bwt->primary) --l;
		cl[0] = cl[1] = cl[2] = cl[3] = 0;
		p = bwt_occ_intv(bwt, k);		//#define bwt_occ_intv(b, k) ((b)->bwt + (k)/OCC_INTERVAL*12)
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
		p += 4;
		// prepare cntk[]
		j = k >> 4 << 4;
#pragma novector
		for (i = k / OCC_INTERVAL * OCC_INTERVAL, x = 0; i < j; i += 16, ++p)
			x += __occ_aux4(bwt, *p);
		y = x;
		x += __occ_aux4(bwt, *p & ~((1U<<((~k&15)<<1)) - 1)) - (~k&15);
		// calculate cntl[] and finalize cntk[]
		j = l >> 4 << 4;
#pragma novector
		for (; i < j; i += 16, ++p) y += __occ_aux4(bwt, *p);
		y += __occ_aux4(bwt, *p & ~((1U<<((~l&15)<<1)) - 1)) - (~l&15);
		memcpy(cntl, cntk, 16);
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
}



double mysecond()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

int main()
{
	double stime, etime, time0, time1;
	int i, j;

	bwt_t bwt;
	bwt.primary = 156168822;
	int len = 1000000000/OCC_INTERVAL*12;

	bwt.bwt = (uint32_t*)malloc(len * sizeof(uint32_t));
	//printf("bwt = %x, len = %d\n", bwt.bwt, len);
		
	for(i = 0; i < len; ++i){
		bwt.bwt[i] = rand();
	}
	
	bwt_gen_cnt_table(&bwt);
	
	int k_array[36] = {119886139, 119886142, 119886160, 13935521, 48936885, 48936898, 48937043, 102311548, 71820894, 30790181, 16904028,  16904184, 51781826, 119886136, 84431695, 119886149, 135287117, 40002292, 127798639, 84431864, 119886139, 119886142, 119886126, 84431780, 112458781, 143782083, 7831661, 7833497, 95767560, 95767898, 37563944, 92282015, 119886132, 112456522, 20032937, 20033008};
	
	int l_array[36] = {119886142, 119886160, 48936898, 13935531, 119886146, 48937047, 119886149, 119886142, 20033008, 25268493, 95140018, 91859897, 119886145, 119886145, 84431935, 84431935, 132918865, 92396349, 104975060, 28858425, 157277902, 99193687, 1162671, 489074, 119886145, 255, 66, 101532506, 101532506, 28349555, 46569004, 48937106, 99193687, 119886132, 62186957, 28858425};
	
	bwtint_t k, l, cntk[4], cntl[4], cntk_opt[4], cntl_opt[4];
	
#if 1
	// Ori
	stime = mysecond();
	for(i = 0; i < count; ++i){
		for(j = 0; j < 36; ++j){
			k = k_array[j];
			l = l_array[j];
			bwt_2occ4(&bwt, k, l, cntk, cntl);
		}
	}
	etime = mysecond();
	time0 = etime - stime;
	printf("Ori time = %f\n", time0);
	fflush(stdout);
	
	// Opt
	stime = mysecond();
	for(i = 0; i < count; ++i){
		for(j = 0; j < 36; ++j){
			k = k_array[j];
			l = l_array[j];
			bwt_2occ4_opt(&bwt, k, l, cntk_opt, cntl_opt);
		}
	}
	etime = mysecond();
	time1 = etime - stime;
	printf("Opt time = %f\n", time1);
	printf("Speedup  = %fx\n", time0/time1);
	fflush(stdout);
#endif

	int errcnt = 0;
	for(i = 0; i < 4; ++i){
		if(cntk[i] != cntk_opt[i]){
			printf("cntk[i] = %10d, cntk_opt[i] = %10d\n", cntk[i], cntk_opt[i]);
			errcnt ++;
		}
		if(cntl[i] != cntl_opt[i]){ 
			printf("cntl[i] = %10d, cntl_opt[i] = %10d\n", cntl[i], cntl_opt[i]);
			errcnt ++;
		}
	}
	printf("errors: %d/8\n", errcnt);
	
	free(bwt.bwt);	
	
	return 0;
}









