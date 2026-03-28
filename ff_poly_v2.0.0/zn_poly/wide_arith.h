/*
    Copyright (C) 2007, 2008, David Harvey

    This file is part of ff_poly.

    ff_poly is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 2 of the License.

    ff_poly is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ff_poly.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
   This file defines macros:
   
      ZNP_MUL_HI(hi, a, b)
         computes the high word of a * b, where a and b are unsigned longs.
         
      ZNP_MUL_WIDE(hi, lo, a, b)
         computes the high and low words of a * b, where a and b are unsigned
         longs.
         
      ZNP_ADD_WIDE(s1, s0, a1, a0, b1, b0)
         computes B*s1 + s0 = (B*a1 + a0) + (B*b1 + b0), where all the
         variables are unsigned longs. Carry is discarded.

      ZNP_SUB_WIDE(s1, s0, a1, a0, b1, b0)
         computes B*s1 + s0 = (B*a1 + a0) - (B*b1 + b0), where all the
         variables are unsigned longs. Borrow is discarded.
   
   If using gcc, there are several inline assembly implementations.
   
*/

#ifndef ZNP_WIDE_ARITH_H
#define ZNP_WIDE_ARITH_H

#include <limits.h>
#if ULONG_MAX == 4294967295U
#define UNSIGNED_LONG_BITS 32
#elif ULONG_MAX == 18446744073709551615U
#define UNSIGNED_LONG_BITS 64
#else
#error wide_arith requires that unsigned long is either 32 bits or 64 bits
#endif
#if (defined (__GNUC__) && (__GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ >= 7)))

//  To simplify things, we require gcc v2.7 or higher.

// ------ POWERPC ------

#if defined (_ARCH_PPC) || defined (__powerpc__) || defined (__POWERPC__)     \
  || defined (__ppc__) || defined (__ppc64__)                                 \
  || (defined (PPC) && ! defined (CPU_FAMILY)) /* gcc 2.7.x GNU&SysV */       \
  || (defined (PPC) && defined (CPU_FAMILY)    /* VxWorks */                  \
         && CPU_FAMILY == PPC)

#if (UNSIGNED_LONG_BITS == 32)

#define ZNP_MUL_HI(hi, a, b)                                                  \
   __asm__ ("mulhwu %0,%1,%2" : "=r" (hi) : "%r" (a), "r" (b));
#define ZNP_ADD_WIDE(s1, s0, a1, a0, b1, b0)                                  \
  do {                                                                        \
    if (__builtin_constant_p (b1) && (b1) == 0)                               \
      __asm__ ("{a%I4|add%I4c} %1,%3,%4\n\t{aze|addze} %0,%2"                 \
        : "=r" (s1), "=&r" (s0) : "r" (a1), "%r" (a0), "rI" (b0));            \
    else if (__builtin_constant_p (b1) && (b1) == ~(unsigned long) 0)         \
      __asm__ ("{a%I4|add%I4c} %1,%3,%4\n\t{ame|addme} %0,%2"                 \
        : "=r" (s1), "=&r" (s0) : "r" (a1), "%r" (a0), "rI" (b0));            \
    else                                                                      \
      __asm__ ("{a%I5|add%I5c} %1,%4,%5\n\t{ae|adde} %0,%2,%3"                \
        : "=r" (s1), "=&r" (s0)                                               \
        : "r" (a1), "r" (b1), "%r" (a0), "rI" (b0));                          \
  } while (0)
#define ZNP_SUB_WIDE(s1, s0, a1, a0, b1, b0)                                  \
  do {                                                                        \
    if (__builtin_constant_p (a1) && (a1) == 0)                               \
      __asm__ ("{sf%I3|subf%I3c} %1,%4,%3\n\t{sfze|subfze} %0,%2"             \
          : "=r" (s1), "=&r" (s0) : "r" (b1), "rI" (a0), "r" (b0));           \
    else if (__builtin_constant_p (a1) && (a1) == ~(unsigned long) 0)         \
      __asm__ ("{sf%I3|subf%I3c} %1,%4,%3\n\t{sfme|subfme} %0,%2"             \
          : "=r" (s1), "=&r" (s0) : "r" (b1), "rI" (a0), "r" (b0));           \
    else if (__builtin_constant_p (b1) && (b1) == 0)                          \
      __asm__ ("{sf%I3|subf%I3c} %1,%4,%3\n\t{ame|addme} %0,%2"               \
          : "=r" (s1), "=&r" (s0) : "r" (a1), "rI" (a0), "r" (b0));           \
    else if (__builtin_constant_p (b1) && (b1) == ~(unsigned long) 0)         \
      __asm__ ("{sf%I3|subf%I3c} %1,%4,%3\n\t{aze|addze} %0,%2"               \
          : "=r" (s1), "=&r" (s0) : "r" (a1), "rI" (a0), "r" (b0));           \
    else                                                                      \
      __asm__ ("{sf%I4|subf%I4c} %1,%5,%4\n\t{sfe|subfe} %0,%3,%2"            \
          : "=r" (s1), "=&r" (s0)                                             \
          : "r" (a1), "r" (b1), "rI" (a0), "r" (b0));                         \
  } while (0)
#elif (UNSIGNED_LONG_BITS == 64)
#define ZNP_MUL_HI(hi, a, b)                                                  \
   __asm__ ("mulhdu %0,%1,%2" : "=r" (hi) : "%r" (a), "r" (b));
#define ZNP_ADD_WIDE(s1, s0, a1, a0, b1, b0)                                  \
  do {                                                                        \
    if (__builtin_constant_p (b1) && (b1) == 0)                               \
      __asm__ ("{a%I4|add%I4c} %1,%3,%4\n\t{aze|addze} %0,%2"                 \
             : "=r" (s1), "=&r" (s0) : "r" (a1), "%r" (a0), "rI" (b0));       \
    else if (__builtin_constant_p (b1) && (b1) == ~(unsigned long) 0)         \
      __asm__ ("{a%I4|add%I4c} %1,%3,%4\n\t{ame|addme} %0,%2"                 \
             : "=r" (s1), "=&r" (s0) : "r" (a1), "%r" (a0), "rI" (b0));       \
    else                                                                      \
      __asm__ ("{a%I5|add%I5c} %1,%4,%5\n\t{ae|adde} %0,%2,%3"                \
             : "=r" (s1), "=&r" (s0)                                          \
             : "r" (a1), "r" (b1), "%r" (a0), "rI" (b0));                     \
  } while (0)

#define ZNP_SUB_WIDE(s1, s0, a1, a0, b1, b0)                                  \
  do {                                                                        \
    if (__builtin_constant_p (a1) && (a1) == 0)                               \
      __asm__ ("{sf%I3|subf%I3c} %1,%4,%3\n\t{sfze|subfze} %0,%2"             \
               : "=r" (s1), "=&r" (s0) : "r" (b1), "rI" (a0), "r" (b0));      \
    else if (__builtin_constant_p (a1) && (a1) == ~(unsigned long) 0)         \
      __asm__ ("{sf%I3|subf%I3c} %1,%4,%3\n\t{sfme|subfme} %0,%2"             \
               : "=r" (s1), "=&r" (s0) : "r" (b1), "rI" (a0), "r" (b0));      \
    else if (__builtin_constant_p (b1) && (b1) == 0)                          \
      __asm__ ("{sf%I3|subf%I3c} %1,%4,%3\n\t{ame|addme} %0,%2"               \
               : "=r" (s1), "=&r" (s0) : "r" (a1), "rI" (a0), "r" (b0));      \
    else if (__builtin_constant_p (b1) && (b1) == ~(unsigned long) 0)         \
      __asm__ ("{sf%I3|subf%I3c} %1,%4,%3\n\t{aze|addze} %0,%2"               \
               : "=r" (s1), "=&r" (s0) : "r" (a1), "rI" (a0), "r" (b0));      \
    else                                                                      \
      __asm__ ("{sf%I4|subf%I4c} %1,%5,%4\n\t{sfe|subfe} %0,%3,%2"            \
               : "=r" (s1), "=&r" (s0)                                        \
               : "r" (a1), "r" (b1), "rI" (a0), "r" (b0));                    \
  } while (0)

#endif

#endif
// ------ ALPHA ------

#if (defined (__alpha) && UNSIGNED_LONG_BITS == 64)

#define ZNP_MUL_HI(hi, a, b)                                                  \
   __asm__ ("umulh %r1,%2,%0" : "=r" (hi) : "%rJ" (a), "rI" (b));

#endif
// ------ IA64 ------

#if (defined (__ia64) && UNSIGNED_LONG_BITS == 64)

#define ZNP_MUL_HI(hi, a, b)                                                  \
   __asm__ ("xma.hu %0 = %1, %2, f0" : "=f" (hi) : "f" (a), "f" (b));
#endif
//  ------ x86 ------

#if ((defined (__i386__) || defined (__i486__)) && UNSIGNED_LONG_BITS == 32)

#define ZNP_MUL_WIDE(hi, lo, a, b)                                            \
   __asm__ ("mull %3" : "=a" (lo), "=d" (hi) : "%0" (a), "rm" (b));
#define ZNP_ADD_WIDE(s1, s0, a1, a0, b1, b0)                                  \
   __asm__ ("addl %5,%k1\n\tadcl %3,%k0"                                      \
           : "=r" (s1), "=&r" (s0)                                            \
           : "0"  ((unsigned long)(a1)), "g" ((unsigned long)(b1)),           \
             "%1" ((unsigned long)(a0)), "g" ((unsigned long)(b0)))
#define ZNP_SUB_WIDE(s1, s0, a1, a0, b1, b0)                                  \
   __asm__ ("subl %5,%k1\n\tsbbl %3,%k0"                                      \
           : "=r" (s1), "=&r" (s0)                                            \
           : "0" ((unsigned long)(a1)), "g" ((unsigned long)(b1)),            \
             "1" ((unsigned long)(a0)), "g" ((unsigned long)(b0)))

#endif
// ------ x86-64 ------

#if (defined (__x86_64__) && UNSIGNED_LONG_BITS == 64)

#define ZNP_MUL_WIDE(hi, lo, a, b)                                            \
   __asm__ ("mulq %3" : "=a" (lo), "=d" (hi) : "%0" (a), "rm" (b));
#define ZNP_ADD_WIDE(s1, s0, a1, a0, b1, b0)                                  \
   __asm__ ("addq %5,%q1\n\tadcq %3,%q0"                                      \
      : "=r" (s1), "=&r" (s0)                                                 \
      : "0"  ((unsigned long)(a1)), "rme" ((unsigned long)(b1)),              \
        "%1" ((unsigned long)(a0)), "rme" ((unsigned long)(b0)))

#define ZNP_SUB_WIDE(s1, s0, a1, a0, b1, b0)                                  \
   __asm__ ("subq %5,%q1\n\tsbbq %3,%q0"                                      \
      : "=r" (s1), "=&r" (s0)                                                 \
      : "0" ((unsigned long)(a1)), "rme" ((unsigned long)(b1)),               \
        "1" ((unsigned long)(a0)), "rme" ((unsigned long)(b0)))
#endif
// ------ MIPS ------

#if (defined (__mips))

#if (UNSIGNED_LONG_BITS == 32)

#define ZNP_MUL_WIDE(hi, lo, a, b)                                            \
   __asm__ ("multu %2,%3" : "=l" (lo), "=h" (hi) : "d" (a), "d" (b));
#elif (UNSIGNED_LONG_BITS == 64)

#define ZNP_MUL_WIDE(hi, lo, a, b)                                            \
   __asm__ ("dmultu %2,%3" : "=l" (lo), "=h" (hi) : "d" (a), "d" (b));
#endif

#endif
//  -------- SPARC --------
#if (defined (__sparc__) && UNSIGNED_LONG_BITS == 32)

#if (defined (__sparc_v9__) || defined (__sparcv9) || \
     defined (__sparc_v8__) || defined (__sparcv8) || defined (__sparclite__))

#define ZNP_MUL_WIDE(hi, lo, a, b)                                            \
   __asm__ ("umul %2,%3,%1;rd %%y,%0"                                         \
            : "=r" (hi), "=r" (lo) : "r" (a), "r" (b));

#endif
#endif

#endif // __GNUC__ 
//  -------- AArch64 --------

#if defined(__aarch64__) && !defined(ZNP_MUL_WIDE)

/* 64x64 -> 128-bit multiply using __uint128_t.  AppleClang compiles this
   to the same mul/umulh instruction pair as hand-written assembly. */
#define ZNP_MUL_WIDE(hi, lo, a, b)                                            \
   do {                                                                       \
      __uint128_t _r = (__uint128_t)(unsigned long)(a)                        \
                     * (unsigned long)(b);                                     \
      (lo) = (unsigned long)_r;                                               \
      (hi) = (unsigned long)(_r >> 64);                                       \
   } while (0)

/* 128-bit add: (s1:s0) = (a1:a0) + (b1:b0).  Compiles to adds/adcs. */
#define ZNP_ADD_WIDE(s1, s0, a1, a0, b1, b0)                                  \
   do {                                                                       \
      __uint128_t _s = ((__uint128_t)(a1) << 64 | (a0))                       \
                     + ((__uint128_t)(b1) << 64 | (b0));                      \
      (s0) = (unsigned long)_s;                                               \
      (s1) = (unsigned long)(_s >> 64);                                       \
   } while (0)

/* 128-bit subtract: (s1:s0) = (a1:a0) - (b1:b0).  Compiles to subs/sbcs. */
#define ZNP_SUB_WIDE(s1, s0, a1, a0, b1, b0)                                  \
   do {                                                                       \
      __uint128_t _s = ((__uint128_t)(a1) << 64 | (a0))                       \
                     - ((__uint128_t)(b1) << 64 | (b0));                      \
      (s0) = (unsigned long)_s;                                               \
      (s1) = (unsigned long)(_s >> 64);                                       \
   } while (0)

#endif
//  -------- generic implementations --------
/*
   If neither ZNP_MUL_HI nor ZNP_MUL_WIDE has an assembly implementation,
   do it in C.

   (algorithm is from GMP's longlong.h)
*/
#if !(defined (ZNP_MUL_HI) || defined (ZNP_MUL_WIDE))

#warning No assembly implementation of wide multiplication available for this \
machine; using generic C code instead.

#define ZNP_MUL_WIDE(hi, lo, a, b)                                            \
   do {                                                                       \
      unsigned long __a = (a);                                                \
      unsigned long __b = (b);                                                \
      unsigned long __mask = (1UL << (UNSIGNED_LONG_BITS/2)) - 1;             \
                                                                              \
      unsigned long __a0 = __a & __mask;                                      \
      unsigned long __a1 = __a >> (UNSIGNED_LONG_BITS/2);                     \
      unsigned long __b0 = __b & __mask;                                      \
      unsigned long __b1 = __b >> (UNSIGNED_LONG_BITS/2);                     \
                                                                              \
      unsigned long __p00 = __a0 * __b0;                                      \
      unsigned long __p01 = __a0 * __b1;                                      \
      unsigned long __p10 = __a1 * __b0;                                      \
      unsigned long __p11 = __a1 * __b1;                                      \
                                                                              \
      __p01 += (__p00 >> (UNSIGNED_LONG_BITS/2));   /* no possible carry */   \
      __p01 += __p10;                                                         \
      if (__p01 < __p10)                                                      \
         __p11 += (1UL << (UNSIGNED_LONG_BITS/2));   /* propagate carry */    \
                                                                              \
      (hi) = __p11 + (__p01 >> (UNSIGNED_LONG_BITS/2));                       \
      (lo) = (__p00 & __mask) + (__p01 << (UNSIGNED_LONG_BITS/2));            \
   } while (0)
#endif
/*
   If only one of ZNP_MUL_HI and ZNP_MUL_WIDE is defined, define the other one
   in terms of that one.
*/
#if defined (ZNP_MUL_HI) && !defined (ZNP_MUL_WIDE)
#define ZNP_MUL_WIDE(hi, lo, a, b)                                            \
   do {                                                                       \
      unsigned long __a = (a), __b = (b);                                     \
      ZNP_MUL_HI ((hi), __a, __b);                                            \
      (lo) = __a * __b;                                                       \
   } while (0)
#endif
#if defined (ZNP_MUL_WIDE) && !defined (ZNP_MUL_HI)
#define ZNP_MUL_HI(hi, a, b)                                                  \
   do {                                                                       \
      unsigned long __dummy;                                                  \
      ZNP_MUL_WIDE ((hi), __dummy, (a), (b));                                 \
   } while (0)
#endif

#if !defined (ZNP_ADD_WIDE)
#define ZNP_ADD_WIDE(s1, s0, a1, a0, b1, b0)                                  \
   do {                                                                       \
      unsigned long __a0 = (a0);                                              \
      unsigned long __temp = __a0 + (b0);                                     \
      (s1) = (a1) + (b1) + (__temp < __a0);                                   \
      (s0) = __temp;                                                          \
  } while (0)
#endif
#if !defined (ZNP_SUB_WIDE)
#define ZNP_SUB_WIDE(s1, s0, a1, a0, b1, b0)                                  \
   do {                                                                       \
      unsigned long __a0 = (a0);                                              \
      unsigned long __b0 = (b0);                                              \
      unsigned long __temp = __a0 - __b0;                                     \
      (s1) = (a1) - (b1) - (__a0 < __b0);                                     \
      (s0) = __temp;                                                          \
  } while (0)
#endif
#endif
// end of file ****************************************************************
