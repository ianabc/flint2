/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"


int fmpz_poly_solve_is_zero_a_root(fmpz_poly_t P,
                                   fmpz_bintvl_t* vec_bintvl, 
                                   slv_info_ptr info)
{
    slong j;
    fmpz_bintvl_t I;
    
    /* Check if 0 is a root */
    if ( fmpz_cmp_ui(fmpz_poly_get_coeff_ptr(P, 0), 0) == 0 ) {
        fmpz_bintvl_init(I);
        I->is_exact = 1;
 
        fmpz_bintvl_new_root(vec_bintvl, I, info);
        info->nb_roots++;

        /* the poly has a smaller degree */
        info->dg--; 
        for (j = 0; j <= info->dg; j++)
        {
            fmpz_set(fmpz_poly_get_coeff_ptr(P, j), fmpz_poly_get_coeff_ptr(P, j+1));
        }
        fmpz_bintvl_clear(I);

        return 1;
    }
    
    return 0;
}




