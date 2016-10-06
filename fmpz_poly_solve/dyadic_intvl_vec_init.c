/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

fmpz_dyadic_intvl_struct *
fmpz_dyadic_intvl_vec_init(slong len)
{
    return  (fmpz_dyadic_intvl_struct*) flint_calloc(len, sizeof(fmpz_dyadic_intvl_struct));
}

