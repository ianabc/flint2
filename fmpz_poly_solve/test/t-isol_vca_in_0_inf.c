/*
  Copyright (C) 2016 Elias Tsigaridas

  This file is part of FLINT.

  FLINT is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_poly_solve.h"

int
main(void)
{
    int iter;
    FLINT_TEST_INIT(state);

    flint_printf("isol_vca_in_0_inf ... ");
    fflush(stdout);

    /* Check aliasing */
    for (iter = 0; iter < 10; iter++) /*1000 * flint_test_multiplier(); iter++)  */
    {
#if 0
        FILE *file_in = stdin;    
        
        char fname[] = "/Users/elias/ET/et-soft/data/rat-5.dat";
        file_in = fopen(fname, "r");
        
        fmpz_poly_t F;
        fmpz_poly_init(F);
        fmpz_poly_fread2(file_in, F);
        fmpz_poly_print_pretty(F, "T"); flint_printf("\n");

/* info */
        slv_info_t info;
        slv_info_init(info);
        
        info->dg = fmpz_poly_degree(F);
        /* flint_printf("dg = %wd \n", info->dg);   */
        if (info->dg <= 2)
        {
            fmpz_poly_clear(F);
            continue;
        }
        
        /* printf("\nf: "); fmpz_poly_print(f); printf("\n\n");  */
        fmpz_bintvl_t* roots = (fmpz_bintvl_t*) flint_malloc(info->dg * sizeof(fmpz_bintvl_t));
        
        /* Isolate the roots using VCA */
        roots = fmpz_poly_solve_isol_vca_in_0_inf(F, info); 
                
        fmpz_poly_solve_print_all_roots(stdout, roots, info->nb_roots); 
             
         
        slv_info_print(info); 
        

        fclose(file_in);
        
        fmpz_poly_clear(F);
    
        continue;
#endif
        
#if 1

        printf("=====================================================\n");
        fmpz_poly_t f;
        fmpz_t b;
        slong k, m;
        slong  i;
        fmpz_poly_t h;
                
        fmpz_init(b);
        
        fmpz_poly_init(f);
        fmpz_poly_set_coeff_si(f, 0, 1);
        fmpz_poly_init(h);
                
        m = WORD_MIN;
        for(i = 1; i <= 5; i++)
        {
            k = n_randint(state, 101);
            fmpz_poly_set_coeff_si(h, 1, k);
            k = n_randint(state, 102);
            fmpz_poly_set_coeff_si(h, 0, -k);
            fmpz_poly_mul(f, f, h);
        }
        fmpz_poly_solve_remove_content_2exp(f);
        
        fmpz_poly_derivative(h, f);
        fmpz_poly_gcd(h, f, h);
        ulong  dd;
        fmpz_poly_pseudo_div(h, &dd, f, h);
                
        /* printf("\nfac:= "); fmpz_poly_print(h); printf(";\n\n");  */
        /* printf("\nfac:= "); fmpz_poly_print_pretty(h, "T"); printf(";\n\n");  */
                
        
        
        slv_info_t info;
        slv_info_init(info);

        info->dg = fmpz_poly_degree(h);
        info->t_dg = fmpz_poly_degree(h);
        
        /* flint_printf("dg = %wd \n", info->dg);   */
        if (info->t_dg <= 2)
        {
                    
            fmpz_poly_clear(f);
            fmpz_poly_clear(h);
            continue;
        }

        /* printf("\nf: "); fmpz_poly_print(f); printf("\n\n");  */
        fmpz_bintvl_t* roots = (fmpz_bintvl_t*) flint_malloc(info->dg * sizeof(fmpz_bintvl_t));
        
        /* Isolate the roots using VCA */
        roots = fmpz_poly_solve_isol_vca_in_0_inf(h, info); 
                
        
        fmpz_poly_solve_print_all_roots(stdout, roots, info->nb_roots); 
         
        slv_info_print(info); 

        /* printf("\nf:= "); fmpz_poly_print(h);; printf(";\n\n");  */
        /* flint_printf("r: %wd \t %wd \n", info->nb_roots, info->dg); */
        
        if ( info->nb_roots > info->t_dg )
        {
            flint_printf("FAIL:\n");
            printf("\nf:= "); fmpz_poly_print_pretty(h, "T"); printf(";\n\n"); 
            /* fmpz_poly_print(f); printf("\n\n"); */
            printf("ERROR \n"); 
            abort();
        }
        
        flint_free(roots);
                
        fmpz_poly_clear(f);
        fmpz_poly_clear(h);
#endif

    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
