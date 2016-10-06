/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

slong fmpz_poly_solve_root_lower_bound_2exp(const fmpz_poly_t F)
{
    slong q1, q2, p, i, j, a0_sgn;
    slong d, len;
    const fmpz * f;

    len = F->length;
    if (len == 0)
        return 0;

    d = len - 1;
    f = F->coeffs;

    a0_sgn = fmpz_sgn(f + 0);

    q1 = WORD_MIN;
    for (i = d; i > 0; i--)
    {
        if ((fmpz_sgn(f + i) == a0_sgn) || (fmpz_sgn(f + i)  == 0))
            continue;

        q2 = WORD_MAX;
        for (j = i - 1; j >= 0; j--)
        {
            if (fmpz_sgn(f + j) != a0_sgn)
                continue;

            p = fmpz_bits(f + i) - fmpz_bits(f + j) - 1;
            q2 = FLINT_MIN(q2, p/(i-j) + 2);
        }

        q1 = FLINT_MAX(q1, q2);
    }

    if (q1 == WORD_MIN)
        return q1 = -1;

    return -(q1+1);
}

