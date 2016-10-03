/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void fmpz_bintvl_init(fmpz_bintvl_ptr a)
{
    a->is_exact = 0;
    a->k 		= 0;
    a->sgn_left = 0;
    a->sign 	= 0;

    fmpz_init(a->c);
}

