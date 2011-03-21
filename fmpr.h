/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#ifndef FMPR_H
#define FMPR_H

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h> 
#include "fmpz.h"


mp_limb_t mpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

void mpn_sqr_basecase(mp_ptr, mp_srcptr, mp_size_t);


#if FLINT64
#define LIMBS64 1
#define LIMBS128 2
#define LIMBS192 3
#define LIMBS256 4
#else
#define LIMBS64 2
#define LIMBS128 4
#define LIMBS192 6
#define LIMBS256 8
#endif

/*
SIGN_EXTENDED would indicate that the exponent is an fmpz.
*/

#define SIGN_NONE 0
#define SIGN_MINUS 1
#define SIGN_EXTENDED 2
#define SIGN_INF 4
#define SIGN_NAN 8

typedef struct { unsigned int sign; fmpz exp; mp_limb_t d[LIMBS64]; } fmpr64;
typedef struct { unsigned int sign; fmpz exp; mp_limb_t d[LIMBS128]; } fmpr128;
typedef struct { unsigned int sign; fmpz exp; mp_limb_t d[LIMBS192]; } fmpr192;
typedef struct { unsigned int sign; fmpz exp; mp_limb_t d[LIMBS256]; } fmpr256;

typedef fmpr64 fmpr64_t[1];
typedef fmpr128 fmpr128_t[1];
typedef fmpr192 fmpr192_t[1];
typedef fmpr256 fmpr256_t[1];

void fmpr64_init(fmpr64_t);
void fmpr128_init(fmpr128_t);
void fmpr192_init(fmpr192_t);
void fmpr256_init(fmpr256_t);

void fmpr64_clear(fmpr64_t);
void fmpr128_clear(fmpr128_t);
void fmpr192_clear(fmpr192_t);
void fmpr256_clear(fmpr256_t);

void fmpr64_debug(fmpr64_t);
void fmpr128_debug(fmpr128_t);
void fmpr192_debug(fmpr192_t);
void fmpr256_debug(fmpr256_t);

void fmpr64_set_ui(fmpr64_t, unsigned long);
void fmpr128_set_ui(fmpr128_t, unsigned long);
void fmpr192_set_ui(fmpr192_t, unsigned long);
void fmpr256_set_ui(fmpr256_t, unsigned long);

void fmpr64_mul(fmpr64_t, const fmpr64_t, const fmpr64_t);
void fmpr128_mul(fmpr128_t, const fmpr128_t, const fmpr128_t);
void fmpr192_mul(fmpr192_t, const fmpr192_t, const fmpr192_t);
void fmpr256_mul(fmpr256_t, const fmpr256_t, const fmpr256_t);

void fixed_print(mp_srcptr x, long limbs, long frac_limbs);

mp_size_t fixed_set_mpfr(mp_ptr y, mpfr_srcptr t, mp_size_t prec);

void fixed_get_mpfr(mpfr_t dest, mp_srcptr src, mp_size_t size, mp_size_t prec);



void
fixed_eval_series_1(mp_ptr y,
                    mp_srcptr coeffs,
                    mp_srcptr denoms,
                    int split_denoms,
                    mp_srcptr x, int limbs, int terms, int sums);


/* TODO: allow caches to be resized at runtime */
#define EXP_CACHE1_BITS 8
#define EXP_CACHE2_BITS 8
#define EXP_CACHE_PREC_BITS 320
#define EXP_CACHE_PREC_LIMBS (EXP_CACHE_PREC_BITS / FLINT_BITS)

extern mp_limb_t exp_cache_1[1 << EXP_CACHE1_BITS][EXP_CACHE_PREC_LIMBS + 1];
extern mp_limb_t exp_cache_2[1 << EXP_CACHE1_BITS][EXP_CACHE_PREC_LIMBS + 1];
extern int exp_cache_initialised;
void exp_cache_init();

void fixed_exp_cache(mp_ptr y, mp_srcptr x, mp_size_t limbs, long tol_bits);
void fixed_exp_mpfr(mp_ptr y, mp_srcptr x, mp_size_t limbs, long tol_bits);
void fixed_exp(mp_ptr y, mp_srcptr x, mp_size_t limbs, long tol_bits);


#define LOG_CACHE_BITS 9
#define LOG_CACHE_PREC_BITS 320
#define LOG_CACHE_PREC_LIMBS (LOG_CACHE_PREC_BITS / FLINT_BITS) 

extern mp_limb_t log_cache[1 << LOG_CACHE_BITS][LOG_CACHE_PREC_LIMBS];
extern int log_cache_initialised;
void log_cache_init();
void fixed_log_cache(mp_ptr y, mp_srcptr x, mp_size_t limbs, long tol_bits);
void fixed_log_mpfr(mp_ptr y, mp_srcptr x, mp_size_t limbs, long tol_bits);
void fixed_log(mp_ptr y, mp_srcptr x, mp_size_t limbs, long tol_bits);

#endif
