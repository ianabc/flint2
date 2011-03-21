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


mp_limb_t log_cache[1 << LOG_CACHE_BITS][LOG_CACHE_PREC_LIMBS];

int log_cache_initialised = 0;

void
log_cache_init()
{
    int i, prec;
    mpfr_t h, logh;

    prec = LOG_CACHE_PREC_LIMBS * FLINT_BITS + 1;

    mpfr_init2(h, prec + 16);
    mpfr_init2(logh, prec + 16);

    for (i = 0; i < (1 << LOG_CACHE_BITS); i++)
    {
        mpfr_set_ui_2exp(h, i, -LOG_CACHE_BITS, GMP_RNDN);
        mpfr_add_ui(h, h, 1, GMP_RNDN);
        mpfr_log(logh, h, GMP_RNDN);
        fixed_set_mpfr(log_cache[i], logh, LOG_CACHE_PREC_LIMBS);
    }

    mpfr_clear(h);
    mpfr_clear(logh);

    log_cache_initialised = 1;
}
