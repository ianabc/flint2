/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

double fmpz_bintvl_get_mid_d(fmpz_bintvl_srcptr I)
{
    mpq_t t;
    fmpq_t m;
    double r;

    mpq_init(t);
    fmpq_init(m);
    fmpz_bintvl_get_mid(m, I);   /* this can be optimised */
    fmpq_get_mpq(t, m);

    r = mpq_get_d(t);

    mpq_clear(t);
    fmpq_clear(m);
    return r;
}
