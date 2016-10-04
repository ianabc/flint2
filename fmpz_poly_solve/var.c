/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"



/* computes the number of sign variations */

slong fmpz_poly_solve_var(const fmpz_poly_t f)
{
	slong i, j;
    slong v = 0;
    slong d = fmpz_poly_degree(f);

    j = 0;
    for ( i=1; i <= d; ++i) {
    	if ( fmpz_sgn(fmpz_poly_get_coeff_ptr(f, i)) * fmpz_sgn(fmpz_poly_get_coeff_ptr(f, j)) < 0 )
          {
    		++v;
    		j = i;
    	}
    }
    return v;
}
