/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"


int fmpz_poly_solve_sgn_eval_at_c_2exp(const fmpz_poly_t P, const fmpz_t c, slong k)
{
    fmpz_t r;
    fmpz_t y;
    slong j;
    
    slong d = fmpz_poly_degree(P);

    if (d < 0)   return 0;
    if (d == 0)  return fmpz_sgn(fmpz_poly_get_coeff_ptr(P, 0));
    
    fmpz_init_set_ui(r, 1);
    fmpz_init_set_ui(y, 1);

    /* fprintf(stderr, "c: ");mpz_out_str(stderr, 10, c); fprintf(stderr, " "); */
    /* debug("k %ld", k); */
    
    if (k <= 0)
    {
        j = fmpz_poly_solve_sgn_eval_at_c(P, c);
        fmpz_clear(r);
        fmpz_clear(y);
        return j;
    }

    fmpz_mul(r, fmpz_poly_get_coeff_ptr(P, d), c);
    for (j = d-1; j >= 1; j--)
    {    
        fmpz_mul_2exp(y, fmpz_poly_get_coeff_ptr(P, j), (d-j)*k);
        fmpz_add(r, r, y);
        fmpz_mul(r, r, c);
    }
    fmpz_mul_2exp(y, fmpz_poly_get_coeff_ptr(P, 0), d*k);
    fmpz_add(r, r, y);
    
    j = fmpz_sgn(r);
    /* fprintf(stderr, "r: ");mpz_out_str(stderr, 10, r); fprintf(stderr, " "); */
  
    fmpz_clear(r);
    fmpz_clear(y);
    
    return j;
}
