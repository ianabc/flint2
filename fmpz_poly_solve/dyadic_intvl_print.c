/*
  Copyright (C) 2016 Elias Tsigaridas

  This file is part of FLINT.

  FLINT is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"


void fmpz_dyadic_intvl_fprint(FILE *stream, const fmpz_dyadic_intvl_t I)
{
	fmpz_t c1, c2;
	slong k1 = 0;
	slong k2 = 0;

	fmpz_init(c1);
	fmpz_init(c2);

	fmpz_init_set(c1, I->c);

	if (I->k <= 0){
		/* mpz_out_str(stream, 10, I->c);  */
	} else {
		k1 = I->k;
	}

	if (I->k <= 0) {
		fmpz_set_ui(c2, 1);
		fmpz_mul_2exp(c2, c2, -I->k);
		fmpz_add(c2, I->c, c2);
	} else {
		fmpz_add_ui(c2, I->c, 1);
		k2 = I->k;
	}

    if (I->sgn_left == 0)
    {
        /* The root is a rational number */
        fprintf(stream, "%10.5f  ~  ", fmpz_dyadic_intvl_get_mid_d(I));
        fprintf(stream, "(");
        fmpz_fprint(stream, c1); fprintf(stream, "/2^%ld", k1);
        fprintf(stream, ", ");
        fmpz_fprint(stream, c1); fprintf(stream, "/2^%ld", k2);
        fprintf(stream, ")");
    }
    else
    {
        fprintf(stream, "%10.5f  ~  ", fmpz_dyadic_intvl_get_mid_d(I)); 
        fprintf(stream, "(");
        fmpz_fprint(stream, c1); fprintf(stream, "/2^%ld", k1);
        fprintf(stream, ", ");
        fmpz_fprint(stream, c2); fprintf(stream, "/2^%ld", k2);
        fprintf(stream, ")");
    }
    fmpz_clear(c1);
    fmpz_clear(c2);
    
	return ;
}


void fmpz_dyadic_intvl_print(const fmpz_dyadic_intvl_t I)
{
    fmpz_dyadic_intvl_fprint(stdout, I);
}

void fmpz_dyadic_intvl_vec_fprint(FILE* stream, fmpz_dyadic_intvl_struct* roots, slong nbr)
{
    slong i;
	fprintf(stream, "\n#roots = %lu \n", nbr);
	fprintf(stream, "[\n");
	for (i = 0; i < nbr; i++) {
		/* fprintf(stream, " %2d   ", roots[i]->sgn_left); */
		fmpz_dyadic_intvl_fprint(stream, roots+i);
		fprintf(stream, "\n");
	}
	fprintf(stream, "]\n");
}

