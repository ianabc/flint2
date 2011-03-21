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
#include "fmpr.h"


#define STACK_ALLOC 128

void
fixed_eval_series_1(mp_ptr y,
                    mp_srcptr coeffs,
                    mp_srcptr denoms,
                    int split_denoms,
                    mp_srcptr x, int limbs, int terms, int sums)
{
    mp_limb_t tmp_stack[STACK_ALLOC];
    mp_ptr tmp, xpow, t, u, swap, s;

    int i, j, tmp_size;
    int tsize;

    tmp_size =  2 * (2 * limbs);     /* temporary products t, u */
    tmp_size += sums * (limbs + 1);  /* sums */
    tmp_size += (1 + sums) * limbs;  /* x powers */

    if (tmp_size > STACK_ALLOC)
        tmp = malloc(tmp_size * sizeof(mp_limb_t));
    else
        tmp = tmp_stack;

    t = tmp;
    u = t + 2*limbs;
    s = u + 2*limbs;
    xpow = s + sums * (limbs + 1);
    mpn_zero(s, sums * (limbs + 1));

#define XPOW_WRITE(_i) (xpow + (sums - (_i) + 1)*limbs)   /* XXX: +1? */
#define XPOW_READ(_i) (XPOW_WRITE(_i) + limbs)
#define T_WRITE (t)
#define T_READ (t + limbs)
#define U_WRITE (u)
#define U_READ (u + limbs)
#define SUM(_i) (s + (_i) * (limbs + 1))

    /* Compute x^2, x^3, ..., x^sums */
    if (sums >= 2)
        mpn_sqr_basecase(XPOW_WRITE(2), x, limbs);
    else
        mpn_copyi(XPOW_READ(1), x, limbs);
    for (i = 3; i <= sums; i++)
        mpn_mul_basecase(XPOW_WRITE(i), XPOW_READ(i-1), limbs, x, limbs);

    mpn_copyi(t + limbs, XPOW_READ(sums), limbs);
    tsize = limbs;

    for (i = 0; i < sums; i++)
        SUM(i)[limbs] = coeffs[i];

    for (i = sums; i < terms; i += sums)
    {
        for (j = 0; j < sums && i + j < terms; j++)
            SUM(j)[limbs] += mpn_addmul_1(SUM(j), T_READ, tsize, coeffs[i + j]);

        if (i + sums != terms)
        {
            mpn_mul_basecase(U_WRITE, XPOW_READ(sums), limbs, T_READ, tsize);
            swap = u; u = t; t = swap;
        }
    }

    if (split_denoms)
    {
        for (i = 0; i < sums; i++)
            mpn_divrem_1(SUM(i), 0, SUM(i), limbs + 1, denoms[i]);

        for (i = 1; i < sums; i++)
        {
            mp_srcptr xp = (i == 1) ? x : XPOW_READ(i);
            mpn_mul_basecase(tmp, SUM(i), limbs + 1, xp, limbs);
            mpn_add(SUM(0), SUM(0), limbs + 1, tmp + limbs, limbs + 1);
        }

        mpn_copyi(y, SUM(0), limbs + 1);
    }
    else
    {
        for (i = 1; i < sums; i++)
        {
            mp_srcptr xp = (i == 1) ? x : XPOW_READ(i);
            mpn_mul_basecase(tmp, SUM(i), limbs + 1, xp, limbs);
            mpn_add(SUM(0), SUM(0), limbs + 1, tmp + limbs, limbs + 1);
        }
        mpn_divrem_1(y, 0, SUM(0), limbs + 1, denoms[0]);
    }

    if (tmp_size > STACK_ALLOC)
        free(tmp);
}
