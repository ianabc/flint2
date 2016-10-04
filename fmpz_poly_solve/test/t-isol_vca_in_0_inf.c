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

    flint_printf("isol_vca_in_0_inf ....");
    fflush(stdout);

    /* Check aliasing */
    for (iter = 0; iter < 3; iter++) //1000 * flint_test_multiplier(); iter++) 
    {
        fmpz_poly_t f;
        
        fmpz_poly_init(f);
        fmpz_poly_randtest(f, state, n_randint(state, 50), 20);

        /* info */
        slv_info_t info;
        slv_info_init(info);
        
        info->dg = fmpz_poly_degree(f);
        /* flint_printf("dg = %wd \n", info->dg);   */
        if (info->dg <= 2)
        {
            fmpz_poly_clear(f);
            continue;
        }
        
        /* printf("\nf: "); fmpz_poly_print(f); printf("\n\n");  */
        fmpz_bintvl_t* roots = (fmpz_bintvl_t*) flint_malloc(info->dg * sizeof(fmpz_bintvl_t));
        
        /* Isolate the roots using VCA */
        roots = fmpz_poly_solve_isol_vca_in_0_inf(f, info); 
                
        
        /*print_roots_all(stdout, roots, info->nb_roots); */
         
        slv_info_print(info); 
        	
        /*
        if ( !fmpz_poly_equal(f, g) )
        {
            flint_printf("FAIL:\n");       
            fmpz_poly_print(f); printf("\n\n");
            printf("ERROR \n"); 
            abort();
        }
        */

        //flint_free(roots);

        fmpz_poly_clear(f);
    }

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
