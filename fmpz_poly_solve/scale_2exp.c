/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"


void __fmpz_poly_solve_scale_2exp_pos(fmpz_poly_t F, slong k)
{
	slong i, p;
    slong deg = fmpz_poly_degree(F);
	p = k;
	for (i = 1; i <= deg; i++, p += k) {
		fmpz_mul_2exp(fmpz_poly_get_coeff_ptr(F, i), fmpz_poly_get_coeff_ptr(F, i), p);
	}
	return;
}

void __fmpz_poly_solve_scale_2exp_neg(fmpz_poly_t F, slong k)
{
	slong i, p;
    slong deg = fmpz_poly_degree(F);
	p = deg * (-k);
   	for ( i = 0; i < deg ; i++, p += k ) {
        fmpz_mul_2exp(fmpz_poly_get_coeff_ptr(F, i), fmpz_poly_get_coeff_ptr(F, i), p);
	}
   	return;
}

int fmpz_poly_solve_scale_2exp(fmpz_poly_t F, slong k)
{
    if (fmpz_poly_degree(F) <= 0) return 0;
    (k > 0 ? __fmpz_poly_solve_scale_2exp_pos(F, k) : __fmpz_poly_solve_scale_2exp_neg(F, k));
	return fmpz_poly_solve_remove_content_2exp(F);
}

