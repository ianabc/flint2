/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"


int fmpz_bisect_dyadic_intvl_bisect(const  fmpz_poly_t P,
                                    fmpz_dyadic_intvl_t I,
                                    slong t)
{
    int sgn_m;
    slong i;
    
    fmpz_t m;
    fmpz_init(m);

    /* fprintf(stderr, "c: ");mpz_out_str(stderr, 10, I->c); fprintf(stderr, " ");  */
    // fprintf(stderr, "2^%ld \n", I->k);
    for (i=0; i < t; i++) {
       /* fprintf(stderr, "bm: "); mpz_out_str(stderr, 10, I->c); fprintf(stderr, " "); */
       /* fprintf(stderr, " b  2^%ld \n", I->k); */
        if (I->k < 0) {
            fmpz_set_ui(m, 1);
            fmpz_mul_2exp(m, m, -I->k - 1);
            fmpz_add(m, m, I->c);
        } else {
            fmpz_mul_ui(m, I->c, 2);
            fmpz_add_ui(m, m, 1);
        }
        
        /* fprintf(stderr, "m: "); mpz_out_str(stderr, 10, m); fprintf(stderr, " "); */
        /* fprintf(stderr, "2^%ld \n", I->k); */
        sgn_m = fmpz_poly_solve_sgn_eval_at_c_2exp(P, m, I->k+1);
        
        
        /* printf("\n\t Signs:  %d %d  %d",  I->sgn_left, sgn_m,  -I->sgn_left); */
        if (sgn_m == 0) {
            fmpz_set(I->c, m);
            I->k++;
            I->sgn_left = 0;
            fmpz_clear(m);
            return 0;
        }

        if (sgn_m == I->sgn_left) {
            fmpz_set(I->c, m);
            I->k++;
        } else if (sgn_m == -I->sgn_left) {
            if (I->k >= 0) { fmpz_mul_ui(I->c, I->c, 2);}
            I->k++;
        } else {
            printf("We should never be here!  %d %d  %d",  I->sgn_left, sgn_m,  -I->sgn_left);
            return -1;
        }
           
    }
    fmpz_clear(m);
    return 1;
}



int fmpz_dyadic_intvl_refine_until( const fmpz_poly_t P, fmpz_dyadic_intvl_t I, slong t)
{
    // t should be negative
    // debug("here");
    /* check_debug(t >= 0, "The width should be small"); */
    long m;

    if (t >= I->k) {
        m = t - I->k;
        /* debug("m : %ld \t %ld", m, I->k); */
        return fmpz_bisect_dyadic_intvl_bisect(P, I, m);
    }
    return 1;
}



void fmpz_dyadic_intvl_vec_refine_until(const fmpz_poly_t P,
                                        fmpz_dyadic_intvl_struct* roots,
                                        long nbr,
                                        long t)
{
  long i;
	for (i = 0; i < nbr; i++) {
        fmpz_dyadic_intvl_refine_until(P, roots+i, t);
	}
}


