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

    2x2 mul code taken from MPFR 2.3.0
    (Copyright (C) 1991-2007 Free Software Foundation, Inc.)

******************************************************************************/

#include "fmpr.h"


mp_limb_t mpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

#define nn_mul_1x1 umul_ppmm

#define nn_mul_2x1(r2, r1, r0, a1, a0, b0)                  \
    do {                                                    \
        mp_limb_t t1;                                       \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
    } while (0)

#define nn_mul_2x2(r3, r2, r1, r0, a1, a0, b1, b0)          \
    do {                                                    \
        mp_limb_t t1, t2, t3;                               \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
        umul_ppmm(t1, t2, a0, b1);                          \
        umul_ppmm(r3, t3, a1, b1);                          \
        add_ssaaaa(r3, t1, r3, t1, 0, t3);                  \
        add_ssaaaa(r2, r1, r2, r1, t1, t2);                 \
        r3 += r2 < t1;                                      \
    } while (0)

#define nn_mul_nxn(d, a, b, n)                                              \
    do {                                                                    \
        switch(n)                                                           \
        {                                                                   \
            case 1:                                                         \
                nn_mul_1x1(d[1], d[0], a[0], b[0]);                         \
                break;                                                      \
            case 2:                                                         \
                nn_mul_2x2(d[3], d[2], d[1], d[0], a[1], a[0], b[1], b[0]); \
                break;                                                      \
            default:                                                        \
                mpn_mul_basecase(d, a, n, b, n);                            \
        }                                                                   \
    } while (0)

/*
    Set (d, dn) = (s, sn) rounded to dn limbs, with downward (truncating)
    rounding. Shift is 1 if result was multiplied by two, 0 otherwise.

    Assumes that (s, sn) is the result of multiplication by
    top-normalised floating-point mantissas, i.e. the second topmost
    bit in the top limb of s is set (unless the result is zero).

    Assumes sn >= dn and sn >= 1.
*/
#define nn_round_prod_down(shift, d, dn, s, sn)                     \
    do {                                                            \
        int i;                                                      \
        shift = !((s)[(sn)-1] >> (FLINT_BITS-1));                   \
        if (shift)                                                  \
        {                                                           \
            for (i = (dn) - 1; i >= 1; i--)                         \
                (d)[i] = ((s)[i+(sn)-(dn)] << 1)                    \
                       | ((s)[i+(sn)-(dn)-1] >> (FLINT_BITS-1));    \
            (d)[0] = ((s)[(sn)-(dn)] << 1);                         \
        }                                                           \
        else                                                        \
        {                                                           \
            for (i = (dn) - 1; i >= 0; i--)                         \
                (d)[i] = (s)[i+(sn)-(dn)];                          \
        }                                                           \
    } while (0)


#define DEFINE(func, type, limbs)                                   \
void func(type z, const type x, const type y)                       \
{                                                                   \
    int shift;                                                      \
    fmpz exp;                                                       \
    mp_limb_t tmp[2*limbs];                                         \
                                                                    \
    /* TODO: here we would handle special values and large exp */   \
    if ((x->sign | y->sign | z->sign) > 1)                          \
        abort();                                                    \
                                                                    \
    nn_mul_nxn(tmp, x->d, y->d, limbs);                             \
    nn_round_prod_down(shift, z->d, limbs, tmp, 2*limbs);           \
                                                                    \
    exp = x->exp + y->exp + limbs*FLINT_BITS + shift;               \
                                                                    \
    /* TODO: here we would promote exp and set the sign bit */      \
    if (exp < COEFF_MIN || exp > COEFF_MAX)                         \
        abort();                                                    \
                                                                    \
    z->sign = x->sign ^ y->sign;                                    \
    z->exp = exp;                                                   \
}

DEFINE(fmpr64_mul, fmpr64_t, LIMBS64)
DEFINE(fmpr128_mul, fmpr128_t, LIMBS128)
DEFINE(fmpr192_mul, fmpr192_t, LIMBS192)
DEFINE(fmpr256_mul, fmpr256_t, LIMBS256)
