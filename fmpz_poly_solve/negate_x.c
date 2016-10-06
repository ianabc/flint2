/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"



void fmpz_poly_negate_x(fmpz_poly_t A)
{
    slong dg, j;

    dg = fmpz_poly_degree(A);
    
    for (j = 1; j <= dg; j += 2)
    {
        fmpz_neg(fmpz_poly_get_coeff_ptr(A, j), fmpz_poly_get_coeff_ptr(A, j));
    }
    
    if (fmpz_sgn(fmpz_poly_get_coeff_ptr(A, dg)) < 0)
    {
        fmpz_poly_neg(A, A);
    }
}


