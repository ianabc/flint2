/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"




fmpz_bintvl_t*  fmpz_poly_solve_isol_vca_in_0_inf(const fmpz_poly_t A, slv_info_ptr info)
{
    fmpz_bintvl_t* roots = (fmpz_bintvl_t*) malloc(info->dg * sizeof(fmpz_bintvl_t)); 

    fmpz_poly_t F;
    fmpz_poly_init(F);
    fmpz_poly_set(F, A);


    /* Check if 0 is a root of the polynomial, and reduce its degree */
    fmpz_poly_solve_is_zero_a_root(F, roots, info); 

    long k = fmpz_poly_solve_root_upper_bound_2exp(F);
    info->bd = k;

    /* printf("\nP2:= "); fmpz_poly_print_pretty(F, "T"); printf(";\n\n");  */
    /* Put the roots in (0, 1) */

	fmpz_poly_solve_scale_2exp(F, k);
  
    /* printf("(inf) f: "); fmpz_poly_print(F); printf("\n\n"); */
    
    /* Isolate the real roots in (0, 1) */
    fmpz_poly_solve_isol_vca_in_0_1(F, roots, info); 
    
    
    fmpz_poly_clear(F);
    
    /* Compute the sign of the poly at the left endpoint of the isolating interval */
    fmpz_poly_solve_adjust_bintvl_signs(A, roots, info); 

    
    return roots;
}

