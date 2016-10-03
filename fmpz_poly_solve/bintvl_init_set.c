/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void fmpz_bintvl_init_set(fmpz_bintvl_ptr a, fmpz_bintvl_srcptr b)
{
    a->is_exact = b->is_exact;
    a->k 		= b->k;
    a->sgn_left = b->sgn_left;
    a->sign 	= b->sign;
    fmpz_init_set(a->c, b->c);
}

