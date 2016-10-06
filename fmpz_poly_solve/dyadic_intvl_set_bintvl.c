/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void fmpz_dyadic_intvl_set_bintvl(fmpz_dyadic_intvl_t a, const fmpz_bintvl_t b)
{  
    a->k = b->k;

    if ( b->is_exact )
    {
        a->sgn_left = 0;
        if ( a->k < 0 ) a->k = 0;
    }
    else
    {
        a->sgn_left = b->sgn_left;
    }
    
    fmpz_set(a->c, b->c);
    
}

