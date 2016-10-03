/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void fmpz_bintvl_new_root(fmpz_bintvl_t * roots, fmpz_bintvl_srcptr I,
        slv_info_ptr info)
{
    int b;
    int j;
    /* REMOVE THIS EVENTUALLY */
    int k;
    fmpz_t c;    
    
    fmpz_init_set(c, I->c);
    k = I->k;
    j = info->nb_roots;  /** The current root */
    b = info->bd; 

    /* debug("#roots: %lu", info->nb_roots); */
    /* debug("\t sign: %d  \t b: %d  \t k: %d", sign, b, k); */

    roots[j]->sign = info->sign;
    fmpz_init(roots[j]->c);

    if (k <= b)
    {
        if (info->sign == -1)
        {
            fmpz_neg(roots[j]->c, c);
            fmpz_sub_ui(roots[j]->c, roots[j]->c, 1);
            fmpz_mul_2exp(roots[j]->c, roots[j]->c, b-k);
        }
        else
        {
            fmpz_mul_2exp(roots[j]->c, c, b-k); 
        }

        roots[j]->k = k - b;
        roots[j]->is_exact = I->is_exact;
    }
    else
    {
        if (info->sign == -1)
        {
            fmpz_neg(roots[j]->c, c);
            fmpz_sub_ui(roots[j]->c, roots[j]->c, 1);
        }
        else
        {
            fmpz_set(roots[j]->c, c);
        }

        roots[j]->k = k - b;
        roots[j]->is_exact = I->is_exact;
    }

    fmpz_clear(c);
    return;
}

