/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

int fmpz_poly_solve_remove_content_2exp(fmpz_poly_t F)
{

    ulong cont, i, z;
    slong deg = fmpz_poly_degree(F);

    if (deg < 0) return 0;
	i = 0; while ( fmpz_sgn(fmpz_poly_get_coeff_ptr(F, i)) == 0 ) i++;
    cont = fmpz_val2(fmpz_poly_get_coeff_ptr(F, i));
    
	for( ; (i <= deg) && cont; i++ )
    {
		if ( fmpz_sgn(fmpz_poly_get_coeff_ptr(F, i)) != 0 )
        {
			z = fmpz_val2(fmpz_poly_get_coeff_ptr(F, i));
			if ( z < cont ) cont = z;
		}
	}
	if ( cont == 0 ) return 0;

	for ( i = 0; i <= deg; i++   )
    {
		fmpz_fdiv_q_2exp(fmpz_poly_get_coeff_ptr(F, i), fmpz_poly_get_coeff_ptr(F, i), cont);
    }
    
	return cont;
}

