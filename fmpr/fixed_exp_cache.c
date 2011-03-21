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

#include <mpir.h>
#include <mpfr.h>
#include "fmpr.h"


static const mp_limb_t exp_coeffs[] =
{
#if FLINT64
    2432902008176640000UL,
    2432902008176640000UL,
    1216451004088320000UL,
    405483668029440000UL,
    101370917007360000UL,
    20274183401472000UL,
    3379030566912000UL,
    482718652416000UL,
    60339831552000UL,
    6704425728000UL,
    670442572800UL,
    60949324800UL,
    5079110400UL,
    390700800UL,
    27907200UL,
    1860480UL,
    116280UL,
    6840UL,
    380UL,
    20UL,
    1UL,
#else
    479001600UL,
    479001600UL,
    239500800UL,
    79833600UL,
    19958400UL,
    3991680UL,
    665280UL,
    95040UL,
    11880UL,
    1320UL,
    132UL,
    12UL,
    1UL,
#endif
};

#define EXP_DIVFREE_MAXTERMS (sizeof(exp_coeffs) / sizeof(mp_limb_t))


static const unsigned char fac_bits[] =
{
    1, 1, 2, 3, 5, 7, 10, 13, 16, 19, 22, 26, 29, 33, 37, 41, 45, 49,
    53, 57, 62, 66, 70, 75, 80, 84, 89, 94, 98, 103, 108, 113, 118, 123,
    128, 133, 139, 144, 149, 154, 160, 165, 170, 176, 181, 187, 192, 198,
    203, 209, 215, 220, 226, 232, 238, 243, 249, 255
};


static __inline__ long
exp_needed_terms(long reduced, long tol_bits)
{
    int i;

    i = 2 + (tol_bits + 1) / (1 + reduced);

    while (reduced*i + fac_bits[i] - 1 > tol_bits) i--;

    return i + 1;
}

static __inline__ int
exp_use_sums(long terms)
{
    if (terms < 9) return 2;
    if (terms < 16) return 3;
    return 4;
}


void fixed_exp_cache(mp_ptr y, mp_srcptr x, mp_size_t limbs, long tol_bits)
{
    mp_limb_t top;
    long terms;
    int sums;

    int i1, i2;

    if (exp_cache_initialised == 0)
        exp_cache_init();

    top = x[limbs - 1];
    i1 = top >> (FLINT_BITS - EXP_CACHE1_BITS);
    i2 = (top >> (FLINT_BITS - (EXP_CACHE1_BITS + EXP_CACHE2_BITS))) & \
        ((1UL<<EXP_CACHE2_BITS) - 1);

    terms = exp_needed_terms(EXP_CACHE1_BITS + EXP_CACHE2_BITS,
                        tol_bits);
    sums = exp_use_sums(terms);

    if (terms <= EXP_DIVFREE_MAXTERMS)
    {
        mp_limb_t t[EXP_CACHE_PREC_LIMBS*2 + 2];
        mp_limb_t u[EXP_CACHE_PREC_LIMBS*2 + 2];

        mpn_copyi(t, x, limbs - 1);
        t[limbs - 1] = (top << (EXP_CACHE1_BITS + EXP_CACHE2_BITS)) >> \
                            (EXP_CACHE1_BITS + EXP_CACHE2_BITS);

        fixed_eval_series_1(u, exp_coeffs, t, limbs, terms, sums);

        /* exp(x1+x2+t) = exp(x1)*exp(x2)*exp(t) */
        mpn_mul_basecase(t, u, limbs + 1,
            exp_cache_1[i1] + (EXP_CACHE_PREC_LIMBS - limbs), limbs + 1);
        mpn_mul_basecase(u, t + limbs, limbs + 1,
            exp_cache_2[i2] + (EXP_CACHE_PREC_LIMBS - limbs), limbs + 1);

        mpn_copyi(y, u + limbs, limbs + 1);
        return;
    }

    printf("fixed_exp: not implemented: %ld terms\n", terms);
    abort();
}
