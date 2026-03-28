#ifndef _ASM_INCLUDE_
#define _ASM_INCLUDE_


/*
    Copyright 2011-2014 Andrew V. Sutherland
    See LICENSE file for license details.
*/

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__x86_64__)

#define _asm_div_q_q(q,r,x,y)               asm ("divq %4" :"=a"(q) ,"=d"(r) : "0"(x), "1"(r), "rm"(y))
#define _asm_mult_1_1(z1,z0,x0,y0)      asm ("mulq %3" :"=a"(z0) ,"=d"(z1) : "a"(x0), "rm"(y0))
#define _asm_mult_2_2_1(z1,z0,x1,x0,y0) asm ("mulq %3" :"=a"(z0) ,"=d"(z1) : "a"(x0), "rm"(y0));(z1)+=(y0)*(x1)
#define _asm_addto_2_2(z1,z0,x1,x0)     asm ("addq %3,%0;adcq %5,%1":"=r"(z0),"=r"(z1): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1):"cc")
#define _asm_addto_2_1(z1,z0,x0)            asm ("addq %3,%0;adcq $0,%1":"=r"(z0),"=r"(z1): "0"(z0), "rim"(x0),  "1"(z1):"cc")
#define _asm_addto_3_3(z2,z1,z0,x2,x1,x0)   asm ("addq %4,%0;adcq %6,%1;adcq %8,%2":"=r"(z0),"=r"(z1),"=r"(z2): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1), "2"(z2), "rim"(x2) :"cc")
#define _asm_addto_3_2(z2,z1,z0,x1,x0)      asm ("addq %4,%0;adcq %6,%1;adcq 0,%2":"=r"(z0),"=r"(z1),"=r"(z2): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1), "2"(z2) :"cc")
#define _asm_subfrom_2_2(z1,z0,x1,x0)       asm ("subq %3,%0;sbbq %5,%1":"=r"(z0),"=r"(z1): "0"(z0), "rim"(x0),  "1"(z1), "rim"(x1):"cc")
// increment needs to propogate the carry - performance not critical here anyway
#define _asm_inc_2(z1,z0)               asm ("addq $1,%0;adcq $0,%1":"=r"(z0),"=r"(z1): "0"(z0), "1"(z1):"cc")
//#define _asm_shiftl_2(z1,z0)              asm ("shlq %0,1;rclq %1,1":"=r"(z0),"=r"(z1): "0"(z0), "1"(z1):"cc")
//#define _asm_shiftr_2(z1,z0)              asm ("shrq %1,1;rcrq %0,1":"=r"(z0),"=r"(z1): "0"(z0), "1"(z1):"cc")

#elif defined(__aarch64__)

/* 64x64 -> 128-bit multiply: (z1:z0) = x0 * y0.  Compiles to mul/umulh. */
#define _asm_mult_1_1(z1,z0,x0,y0) \
    do { __uint128_t _r = (__uint128_t)(unsigned long)(x0) * (unsigned long)(y0); \
         (z0) = (unsigned long)_r; (z1) = (unsigned long)(_r >> 64); } while(0)

/* 64x64 -> 128-bit multiply then accumulate: (z1:z0) = x0*y0, z1 += y0*x1 */
#define _asm_mult_2_2_1(z1,z0,x1,x0,y0) \
    do { __uint128_t _r = (__uint128_t)(unsigned long)(x0) * (unsigned long)(y0); \
         (z0) = (unsigned long)_r; (z1) = (unsigned long)(_r >> 64); \
         (z1) += (unsigned long)(y0) * (unsigned long)(x1); } while(0)

/* 128-bit add: (z1:z0) += (x1:x0).  Compiles to adds/adcs. */
#define _asm_addto_2_2(z1,z0,x1,x0) \
    do { __uint128_t _s = ((__uint128_t)(z1) << 64 | (z0)) + ((__uint128_t)(x1) << 64 | (x0)); \
         (z0) = (unsigned long)_s; (z1) = (unsigned long)(_s >> 64); } while(0)

/* Add 64-bit to 128-bit: (z1:z0) += x0 */
#define _asm_addto_2_1(z1,z0,x0) \
    do { __uint128_t _s = ((__uint128_t)(z1) << 64 | (z0)) + (unsigned long)(x0); \
         (z0) = (unsigned long)_s; (z1) = (unsigned long)(_s >> 64); } while(0)

/* 192-bit add: (z2:z1:z0) += (x2:x1:x0) */
#define _asm_addto_3_3(z2,z1,z0,x2,x1,x0) \
    do { __uint128_t _s = (__uint128_t)(z0) + (unsigned long)(x0); \
         (z0) = (unsigned long)_s; \
         _s = (__uint128_t)(z1) + (unsigned long)(x1) + (_s >> 64); \
         (z1) = (unsigned long)_s; \
         (z2) += (unsigned long)(x2) + (unsigned long)(_s >> 64); } while(0)

/* 192-bit add of 128-bit value: (z2:z1:z0) += (x1:x0) */
#define _asm_addto_3_2(z2,z1,z0,x1,x0) \
    do { __uint128_t _s = (__uint128_t)(z0) + (unsigned long)(x0); \
         (z0) = (unsigned long)_s; \
         _s = (__uint128_t)(z1) + (unsigned long)(x1) + (_s >> 64); \
         (z1) = (unsigned long)_s; \
         (z2) += (unsigned long)(_s >> 64); } while(0)

/* 128-bit subtract: (z1:z0) -= (x1:x0).  Compiles to subs/sbcs. */
#define _asm_subfrom_2_2(z1,z0,x1,x0) \
    do { __uint128_t _s = ((__uint128_t)(z1) << 64 | (z0)) - ((__uint128_t)(x1) << 64 | (x0)); \
         (z0) = (unsigned long)_s; (z1) = (unsigned long)(_s >> 64); } while(0)

/* Increment 128-bit value: (z1:z0) += 1 */
#define _asm_inc_2(z1,z0) \
    do { __uint128_t _s = ((__uint128_t)(z1) << 64 | (z0)) + 1; \
         (z0) = (unsigned long)_s; (z1) = (unsigned long)(_s >> 64); } while(0)

/* 128/64 division: q = (r:x) / y, r = (r:x) % y.  Not used outside asm.h. */
#define _asm_div_q_q(q,r,x,y) \
    do { __uint128_t _n = ((__uint128_t)(r) << 64) | (unsigned long)(x); \
         (q) = (unsigned long)(_n / (unsigned long)(y)); \
         (r) = (unsigned long)(_n % (unsigned long)(y)); } while(0)

#else
#error "asm.h: unsupported architecture (need x86_64 or aarch64)"
#endif


#define _asm_mult_3_2_1(z2,z1,z0,x1,x0,y0)  { register unsigned long __u; \
                                       _asm_mult_1_1 (__u,z0,x0,y0); \
                                       _asm_mult_1_1 (z2,z1,x1,y0); \
                                       _asm_addto_2_1 (z2,z1,__u); }

// This function assumes that x[1] and y[1] < 2^31
static inline void _asm_mult_3_2_2 (unsigned long z[3], unsigned long x[2], unsigned long y[2])
{
    register unsigned long U, V, R0,R1;
    
    R1 = 0;
    _asm_mult_1_1(R0,z[0],x[0],y[0]);
    _asm_mult_1_1(U,V,x[0],y[1]);
    _asm_addto_2_2(R1,R0,U,V);
    _asm_mult_1_1(U,V,x[1],y[0]);
    _asm_addto_2_2(R1,R0,U,V);
    z[1] = R0;
    z[2] = x[1]*y[1]+R1;
}

// This function assumes that x[1] and y[1] < 2^31
static inline void _asm_mult_3_2_2r (unsigned long z2, unsigned long z1, unsigned long z0, unsigned long x[2], unsigned long y[2])
{
    register unsigned long U, V, R0,R1;
    
    R1 = 0;
    _asm_mult_1_1(R0,z0,x[0],y[0]);
    _asm_mult_1_1(U,V,x[0],y[1]);
    _asm_addto_2_2(R1,R0,U,V);
    _asm_mult_1_1(U,V,x[1],y[0]);
    _asm_addto_2_2(R1,R0,U,V);
    z1 = R0;
    z2 = x[1]*y[1]+R1;
}


// This function assumes that x[1] < 2^31.  For no obvious reason, this is slower than multiplying?!
static inline void _asm_square_3_2 (unsigned long z[3], unsigned long x[2])
{
    register unsigned long U, V, R0,R1;
    
    R1 = 0;
    _asm_mult_1_1(R0,z[0],x[0],x[0]);
    _asm_mult_1_1(U,V,x[0],x[1]);
    _asm_addto_2_2(R1,R0,U,V);
    _asm_addto_2_2(R1,R0,U,V);
    z[1] = R0;
    z[2] = x[1]*x[1]+R1;
}

#ifdef __cplusplus
}
#endif

#endif
