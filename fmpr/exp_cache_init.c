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


mp_limb_t exp_cache_1[1 << EXP_CACHE1_BITS][EXP_CACHE_PREC_LIMBS + 1];
mp_limb_t exp_cache_2[1 << EXP_CACHE1_BITS][EXP_CACHE_PREC_LIMBS + 1];

int exp_cache_initialised = 0;

void
exp_cache_init()
{
    int i, prec;
    mpfr_t h, exph;

    prec = EXP_CACHE_PREC_LIMBS * FLINT_BITS + 1;

    mpfr_init2(h, prec + 16);
    mpfr_init2(exph, prec + 16);

    mpfr_set_ui_2exp(h, 1, -EXP_CACHE1_BITS, GMP_RNDN);
    mpfr_exp(exph, h, GMP_RNDN);
    mpfr_set_ui(h, 1, GMP_RNDN);

    for (i = 0; i < (1 << EXP_CACHE1_BITS); i++)
    {
        fixed_set_mpfr(exp_cache_1[i], h, EXP_CACHE_PREC_LIMBS);
        mpfr_mul(h, h, exph, GMP_RNDN);
    }

    mpfr_set_ui_2exp(h, 1, -EXP_CACHE1_BITS-EXP_CACHE2_BITS, GMP_RNDN);
    mpfr_exp(exph, h, GMP_RNDN);
    mpfr_set_ui(h, 1, GMP_RNDN);

    for (i = 0; i < (1 << EXP_CACHE2_BITS); i++)
    {
        fixed_set_mpfr(exp_cache_2[i], h, EXP_CACHE_PREC_LIMBS);
        mpfr_mul(h, h, exph, GMP_RNDN);
    }

    mpfr_clear(h);
    mpfr_clear(exph);

    exp_cache_initialised = 1;
}
