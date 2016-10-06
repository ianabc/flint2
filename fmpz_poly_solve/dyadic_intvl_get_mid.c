/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void fmpz_dyadic_intvl_get_mid(fmpq_t m, const fmpz_dyadic_intvl_t I)
{
    fmpz_t t;
    fmpz_init_set(t, I->c);

    if (I->k >= 0)
    {
        fmpz_mul_2exp(t, t, 1);
        fmpz_add_ui(t, t, 1);
        fmpz_set(fmpq_numref(m), t);
        fmpz_set_ui(t, 1);
        fmpz_mul_2exp(t, t, I->k + 1);
        fmpz_set(fmpq_denref(m), t);
    }
    else
    {
        fmpz_set_ui(t, 1);
        fmpz_mul_2exp(t, t, -I->k - 1);
        fmpz_add(t, t, I->c);
        fmpz_set(fmpq_numref(m), t);
    }

    /* this could be optimised */
    fmpq_canonicalise(m);

    fmpz_clear(t);
}

