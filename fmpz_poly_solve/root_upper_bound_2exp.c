/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

slong fmpz_poly_solve_root_upper_bound_2exp(const fmpz_poly_t F)
{
    slong q1, q2, p, i, j, d, len, ad_sgn;
    const fmpz * f;

    len = F->length;
    if (len == 0)
        return 0;

    d = len - 1;
    f = F->coeffs;
    
    ad_sgn = fmpz_sgn(f + d);

    q1 = WORD_MIN;
    for (i = 0; i < d; i++)
    {
        if ((fmpz_sgn(f + i)  == ad_sgn) || (fmpz_sgn(f + i)  == 0))
            continue;

        q2 = WORD_MAX;
        for (j = i + 1; j <= d; j++)
        {
            if (fmpz_sgn(f + j) != ad_sgn)
                continue;

            p = fmpz_bits(f + i) - fmpz_bits(f + j) - 1;
            q2 = FLINT_MIN(q2, p/(j-i) + 2);
        }

        q1 = FLINT_MAX(q1, q2);
    }

    if (q1 == WORD_MIN)
        return q1 = -1;

    return q1+1;
}

