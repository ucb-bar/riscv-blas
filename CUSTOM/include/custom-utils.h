#ifndef CUSTOM_UTILS_H
#define CUSTOM_UTILS_H

// Header file for custom utility functions

#include <stdint.h>

#define MAX(a, b)  ((a>b)?a:b)
void hwacha_init();
int rdcycle();
int rdinstret();
void* safe_malloc(int size);
/* void printfloatmatrix(int channels, int width, int height, float* M); */
/* void printintmatrix(int channels, int width, int height, int* M); */
/* void printint16matrix(int channels, int width, int height, int16_t* M); */
/* void fill_seq_32(float* p, int n, int mode); */
/* void fill_seq_16(int16_t* p, int n, int mode); */
void setvcfg(int nd, int nw, int nh, int np);
int setvlen(int vlen);
int setvlen32(int vlen);

#define PRELOAD(label) ({                                       \
      void* out;                                                \
      asm volatile ("la %0, " label : "=r" (out));              \
      asm volatile ("lw t0, 0(%0)" : : "r" (out) : "t0");       \
      out;                                                      \
    })

#define MEMTOUCH(iaddr, type, bound) ({                           \
      volatile type* addr = iaddr;                                \
      volatile type t;                                            \
      t = (addr)[0];                                              \
      (addr)[0] = t;                                              \
      volatile type* tf = (type*) (((((uintptr_t) (addr)) >> 12) + 1) << 12);    \
      for (; tf - (addr) < bound; tf += (1 << 12) / sizeof(type)) {     \
        t = tf[0];                                                      \
        tf[0] = t;                                                      \
      }                                                                 \
    })

#define VF(label) {                             \
    void* dest;                                   \
    asm volatile ("la %0, " label : "=r" (dest)); \
    asm volatile ("vf 0(%0)" : : "r" (dest));     \
  }
/* void memcpy_16(int16_t* src, int16_t* dest, int len); */
/* void memcpy_32(float* src, float* dest, int len); */

/* void cvt_32_16(float* src, int16_t* dest, int len); */
/* void cvt_16_32(int16_t* src, float* dest, int len); */

/* void gather_16(const int* id, const int16_t* src, int16_t* dest, int len); */
/* void gather_32(const int* id, const float* src, float* dest, int len); */
/* void fill_16(int N, float ALPHA, int16_t* X); */
/* void fill_32(int N, float ALPHA, float* X); */
/* void normalize_16(int16_t *x, int16_t *mean, int16_t *variance, int filters, int spatial); */
/* void normalize_32(float *x, float *mean, float *variance, int filters, int spatial); */
/* void scale_16(int16_t* x, int16_t scale, int size); */
/* void scale_32(float* x, float scale, int size); */
/* void add_16(int16_t* x, int16_t y, int size); */
/* void add_32(float* x, float y, int size); */

/* void square_32(int N, float* x, float* dest); */

/* void axpy_32(int N, float A, float* X, float* Y); */

#endif
