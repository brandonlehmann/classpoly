/*
    Copyright 2007-2016 Andrew V. Sutherland

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

#ifndef _FFPOLYBIG_INCLUDE_
#define _FFPOLYBIG_INCLUDE_

#include "ff.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    ff_poly interface to zn_poly is via ffpolybig
*/
void ff_poly_zn_mod_setup (unsigned long p, int d); // called automatically
void ff_poly_zn_mod_clear ();                       // not called automatically

void ff_poly_square_big (ff_t h[], ff_t f[], int d);
void ff_poly_mult_big (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b);

// these are used for repeated middle products a*b1, a*b2, ... 
// setup1 returns ctx pointer that should be passed to execute1 for each middle product and then passed to clear1 to free up resources
void *ff_poly_middle_product_setup1 (ff_t a[], int d_a, int d_b);
void ff_poly_middle_product_execute1 (ff_t c[], ff_t b[], void *ctx);
void ff_poly_middle_product_clear1 (void *ctx);

// computes c[i] = sum_{j+k=d_b+i} a[j]*b[k] for i in [0,d_a-d_b]
void ff_poly_middle_product_big (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b);

void ff_poly_mult_bigx (ff_t c[], ff_t a[], int d_a, ff_t b[], int d_b);        // assumes integer format, no conversion out of montgomery
void ff_poly_square_bigx (ff_t h[], ff_t f[], int d);                           // assumes integer format, no conversion out of montgomery

#ifdef __cplusplus
}
#endif

#endif
