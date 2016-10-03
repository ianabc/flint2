/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

int fmpz_poly_solve_sgn_eval_at_c(const fmpz_poly_t P, const fmpz_t c)
{
    fmpz_t r;
    slong s;

    fmpz_poly_evaluate_fmpz(r, P, c);
    s = fmpz_sgn(r);
    fmpz_clear(r);

    return s;
}

