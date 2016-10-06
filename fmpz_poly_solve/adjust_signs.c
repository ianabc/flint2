/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void fmpz_poly_solve_adjust_bintvl_signs(const fmpz_poly_t P,
                                         fmpz_bintvl_t* roots,
                                         slv_info_srcptr info)
{
    int s;
    slong i;
    
    for (i = 0; i < info->nb_roots; i++)
    {   
        if ( roots[i]->is_exact ) { continue; }
        
        if ( roots[i]->k > 0 )
        {
            s = fmpz_poly_solve_sgn_eval_at_c_2exp(P, roots[i]->c, roots[i]->k);
        }
        else
        {
            s = fmpz_poly_solve_sgn_eval_at_c(P, roots[i]->c);
        }
        
        roots[i]->sgn_left = s;
    }
}

