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


#endif
