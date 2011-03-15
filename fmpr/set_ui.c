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

#include "fmpr.h"

#define DEFINE(func, type, limbs)                                   \
void func(type x, unsigned long c)                                  \
{                                                                   \
    int i;                                                          \
                                                                    \
    if (x->sign & SIGN_EXTENDED)                                    \
        fmpz_clear(x->exp);                                         \
    x->sign = SIGN_NONE;                                            \
                                                                    \
    for (i = 0; i < limbs - 1; i++)                                 \
        x->d[i] = 0;                                                \
                                                                    \
    i = FLINT_BITS - FLINT_BIT_COUNT(c);                            \
    x->d[limbs-1] = (c << i);                                       \
    x->sign = SIGN_NONE;                                            \
    x->exp = -((long)i) - (limbs-1)*FLINT_BITS;                     \
}


DEFINE(fmpr64_set_ui, fmpr64_t, LIMBS64)
DEFINE(fmpr128_set_ui, fmpr128_t, LIMBS128)
DEFINE(fmpr192_set_ui, fmpr192_t, LIMBS192)
DEFINE(fmpr256_set_ui, fmpr256_t, LIMBS256)
