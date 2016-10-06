/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"


int fmpz_poly_fread2(FILE * file, fmpz_poly_t poly)
{
    int r;
    slong i, len;
    mpz_t t;

    mpz_init(t);

    /* We are reading the degree. Let length is len = dg + 1; */
    r = mpz_inp_str(t, file, 10);
    mpz_add_ui(t, t, 1);
    if (r == 0)
    {
        mpz_clear(t);
        return 0;
    }
    if (!mpz_fits_slong_p(t))
    {
        flint_printf("Exception (fmpz_poly_fread). Length does not fit into a slong.\n");
        abort();
    }
    len = flint_mpz_get_si(t);
    mpz_clear(t);

    fmpz_poly_fit_length(poly, len);

    for (i = 0; i < len; i++)
    {
        r = fmpz_fread(file, poly->coeffs + i);
        if (r <= 0)
            return r;
    }

    _fmpz_poly_set_length(poly, len);
    _fmpz_poly_normalise(poly);

    return 1;
}




