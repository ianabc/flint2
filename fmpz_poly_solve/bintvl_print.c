/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void fmpz_bintvl_print(fmpz_bintvl_srcptr I)
{
    flint_printf("I: ");
    fmpz_print(I->c);
    flint_printf("/2^%wd \n", I->k);
}

