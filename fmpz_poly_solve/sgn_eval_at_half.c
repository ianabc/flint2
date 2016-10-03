/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

/* Returns the sign of F(1/2) */
int fmpz_poly_solve_sgn_eval_at_half (const fmpz_poly_t P)
{
	slong j, p;
	int ret;
	fmpz_t x, y;
    
    slong deg = fmpz_poly_degree(P);
    
	fmpz_init(y);
	fmpz_init_set_ui(x, 0);

	p = deg + 1;
	for (j = 0; j <= deg; j++, p--)
    {
		fmpz_mul_2exp(y, fmpz_poly_get_coeff_ptr(P, j), p);
		fmpz_add(x, x, y);
	}

	fmpz_clear(y);
	ret = fmpz_sgn(x);
	fmpz_clear(x);
	return ret;
}
