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

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "fmpr.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("fixed_set/get_mpfr....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        mp_size_t size;
        mpfr_t t, u, v;
        mp_ptr x;
        int bits, prec, result;

        bits = n_randint(state, 1000) + 2;
        prec = n_randint(state, 2000) / FLINT_BITS + 1;

        mpfr_init2(t, bits);
        mpfr_init2(u, bits);
        mpfr_init2(v, bits);
        x = malloc(sizeof(mp_limb_t) * prec);

        mpfr_urandom(t, state->gmp_state, GMP_RNDD);

        size = fixed_set_mpfr(x, t, prec);
        fixed_get_mpfr(u, x, size, prec);

        mpfr_sub(v, t, u, GMP_RNDN);
        mpfr_abs(v, v, GMP_RNDN);

        result = (mpfr_cmp_ui_2exp(v, 1, -FLINT_BITS*prec) < 0);

        if (!result)
        {
            printf("FAIL\n");
            mpfr_printf("%.20Rg %.20Rg %.20Rg bits=%d prec=%d result=%d\n",
                t, u, v, bits, prec, result);
            abort();
        }

        mpfr_clear(t);
        mpfr_clear(u);
        mpfr_clear(v);
        free(x);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
