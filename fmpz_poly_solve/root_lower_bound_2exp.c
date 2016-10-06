/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"


/* #define fmpz_bitsize(a) fmpz_sizeinbase((a), 2) */

slong fmpz_poly_solve_root_lower_bound_2exp(const fmpz_poly_t F)
{
    slong q1, q2, p, i, j, a0_sgn;
    slong d = fmpz_poly_degree(F);

    if (d <= 0) return 0;
    
    a0_sgn = fmpz_sgn(fmpz_poly_get_coeff_ptr(F, 0));

    q1 = WORD_MIN;
    for (i = d; i > 0; --i) {
	if ( (fmpz_sgn(fmpz_poly_get_coeff_ptr(F, i)) == a0_sgn) || 
             (fmpz_sgn(fmpz_poly_get_coeff_ptr(F, i))  == 0) ) continue;

    	q2 = WORD_MAX;
    	for (j = i-1; j >= 0; --j)
        {
			if ( fmpz_sgn(fmpz_poly_get_coeff_ptr(F, j)) != a0_sgn ) continue;
            p = fmpz_sizeinbase(fmpz_poly_get_coeff_ptr(F, i), 2) -
                fmpz_sizeinbase(fmpz_poly_get_coeff_ptr(F, j), 2) - 1;
            q2 = FLINT_MIN(q2, p/(i-j) +2);
    	}
    	q1 = FLINT_MAX(q1, q2);
    }
    if ( q1 == LONG_MIN )
    	return q1 = -1;
    return -(q1+1);
}


