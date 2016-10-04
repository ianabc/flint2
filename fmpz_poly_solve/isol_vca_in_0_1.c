/*
  Copyright (C) 2016 Elias Tsigaridas

  This file is part of FLINT.

  FLINT is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"



/* 
   Upper bound on the number of positive real roots using Descartes' rule of sign/
   Optimized procedure using various hacks from [Rouillier,Zimmermann:JCAM:2004]
*/ 

long fmpz_poly_solve_Descartes_test ( fmpz_poly_t P,
                                      fmpz_poly_t Q,
                                      long sigh, 
                                      long* status,
                                      slv_info_ptr info )
{  

    unsigned long V = 0;	/** number of sign variations. Initially 0 */
    long i, j, s, t;

    long deg = fmpz_poly_degree(P);

    /* printf("DT P: "); mpzx_pretty_print (stdout, P); */
    info->nb_nodes++;
    /* Prune the computation if all the coefficients are of the sign of P[deg]
       In that case any subsequent interval shall have the same property,
       we put *status = -1  and inform VCA_recursive.
    */
    j = deg; t = fmpz_sgn(fmpz_poly_get_coeff_ptr(P, deg));


    /* poly_print_list(stdout, P, deg); */
    /* Check if all the coefficients have the same sign */
    while ( (j >= 0) && ( (fmpz_sgn(fmpz_poly_get_coeff_ptr(P, j)) == t) || 
                          (fmpz_sgn(fmpz_poly_get_coeff_ptr(P, j)) == 0) ) )  {
        /* printf("%ld - %ld ", t, j); mpz_out_str(stdout, 10, P[j]); printf("\n");  */
        j--;
    }
    if ( j < 0 ) {
        info->nb_pos_hack_1++;
        *status = -1 ; /* ALL_COEFF_POS; */
        return V;
    }

    /* Copy P to Q */
    fmpz_poly_set(Q, P);
    // mpz_poly_set_unsafe(Q, P);

    /* Start computing  Q(1/(x+1)) with the naive algorithm */
    /* The 1st iteration gives the evaluation at x=1, i.e. P(1) */
    for ( j = 0; j <= deg-1; j++ ) { fmpz_add(fmpz_poly_get_coeff_ptr(Q, j+1), fmpz_poly_get_coeff_ptr(Q, j+1), fmpz_poly_get_coeff_ptr(Q, j)); }

    /* s is the sign of P(1) */
    s = fmpz_sgn(fmpz_poly_get_coeff_ptr(Q, deg));

    /* signs at 0  1/2  1 */
    /* s = sgn(P(1))
       All the roots are in (0, 1).
       The sign of P[deg] is the the sign of P(infty) and P(1).
       status = 1 if sgn(P(0)) = sgn(P(1/2)) <> sgn(P(1)) 
       Possibility of 2 roots in (0,1/2) and (1/2, 1) 
       status = 0 otherwise
    */
    *status = s && (s == fmpz_sgn(fmpz_poly_get_coeff_ptr(P, 0))) && (s == -sigh);
    /* fprintf(stderr, "sign(P[0]) = %d sign(P(1/2)) = %ld sign(P(1)) = %ld, flag = %ld\n", mpz_sgn(P->coeffs[0]), sigh, s, *status); */

    for ( i = 1; i <= deg-1; i++ ) {
        /* Prune the computation if all further coefficients are of the sign of Q[deg-i] */
        j = deg - i; t = s;
        /* Starting from j, find the first non-vanishing coeff of Q */
        while ( (j >= 0) && (t == 0) ) { t = fmpz_sgn(fmpz_poly_get_coeff_ptr(Q, j)); j--; }
        /* Check if all the remaining coeffs have the same sign,
           in order to exit early */
        while ( (j >= 0) && ( (fmpz_sgn(fmpz_poly_get_coeff_ptr(Q, j)) == t) ||
                              (fmpz_sgn(fmpz_poly_get_coeff_ptr(Q, j)) == 0) ) ) {  j--;  }
        if ( j < 0 ) {
            info->nb_pos_hack_2++;
            /* printf("THIS IS OPT 2 (%ld) \n", i); */
            /* fprintf(stderr, "Pruning at i=%ld, V = %lu\n", i, V); */
            return V;
        }

        /* Perform one more iteration of the naive algo to compute Q(1/(X+1)) */
        for ( j = 0; j <= deg - i - 1; j++ ) { fmpz_add(fmpz_poly_get_coeff_ptr(Q, j+1), fmpz_poly_get_coeff_ptr(Q, j+1), fmpz_poly_get_coeff_ptr(Q, j)); }

        if ( s == 0 ) {
            s = fmpz_sgn(fmpz_poly_get_coeff_ptr(Q, deg-i));
        } else {
            if ( s == -fmpz_sgn(fmpz_poly_get_coeff_ptr(Q, deg-i)) ) {
                if ( ((V == 1) && !*status) || (V == 2) ) {
                    return (V + 1);
                }
                V++; s = -s;
            }
        }
    }
    if ( s == -fmpz_sgn(fmpz_poly_get_coeff_ptr(Q, 0)) ) V++;

    return V;


}


/* Memory efficient variant. Perform all the computations with 2 polynomials.  */
void fmpz_poly_solve_isol_vca_in_0_1(fmpz_poly_t FF, 
                                     fmpz_bintvl_t* roots, 
                                     slv_info_ptr info)
{
    /* printf("FF: "); fmpz_poly_print(FF); printf("\n\n");  */
    long V;
    long k = 0;
    
    fmpz_t c;
    fmpz_t one;
    fmpz_one(one);

    long status = 0;
    long shalf= 1;
 
    long dg = fmpz_poly_degree(FF);
    info->dg = dg; 

    fmpz_lst_bintvl  queue;

    fmpz_bintvl_ptr I;
    //I = malloc(sizeof(slv_bintvl_t));

    fmpz_poly_t P;
    fmpz_poly_init2(P, dg+1);
    fmpz_poly_set(P, FF);

    fmpz_poly_t Q;
    fmpz_poly_init2(Q, dg+1);
    fmpz_poly_set(Q, P);


    V = fmpz_poly_solve_Descartes_test(P, Q, shalf, &status, info);
    if ( V == 0 ) {
        /* TODO: clear the memory */
        fmpz_poly_clear(P);
        fmpz_poly_clear(Q);

        return ;
    }
    if ( V == 1 ) {
		

        I = flint_malloc(sizeof(fmpz_bintvl_t));
        fmpz_bintvl_init(I);
        fmpz_bintvl_new_root(roots, I, info);
        info->nb_roots++;
        fmpz_bintvl_clear(I);
        flint_free(I);
        fmpz_poly_clear(P);
        fmpz_poly_clear(Q);

        return;
    }
    /* The poly has more than one sign variations.
     * Split and initialize the queue.  */
    FLINT_SLIST_INIT(&queue);
    info->max_depth = 1;

    fmpz_bintvl_ptr Il = flint_malloc(sizeof(fmpz_bintvl_t));
    fmpz_bintvl_init(Il);
    Il->k = 1;

    fmpz_bintvl_ptr Ir = flint_malloc(sizeof(fmpz_bintvl_t));
    fmpz_bintvl_init(Ir);
    Ir->k = 1;
    fmpz_set_ui(Ir->c, 1);

    FLINT_SLIST_PUSH(queue, Ir);
    FLINT_SLIST_PUSH(queue, Il);

    k = 0;
    fmpz_init_set_ui(c, 0);
    while ( (!FLINT_SLIST_EMPTY(&queue)) && info->max_depth <= 100000)  {
        FLINT_SLIST_POP(queue, I);
        info->max_depth = FLINT_MAX(info->max_depth, I->k);

        /* Construct the polynomial depending on the previous k */
        if ( k < I->k ) {
            fmpz_poly_solve_scale_2exp(P, -1);
            info->nb_homo++;
            k++;
        } else if ( k == I->k) {
            fmpz_poly_taylor_shift(P, P, one);
            // fmpz_poly_taylor_shift_by_1(P, P);
            info->nb_trans++;
        } else if ( k > I->k ) {
            fmpz_poly_taylor_shift(P, P, one);
            //fmpz_poly_taylor_shift_by_1(P, P);
            info->nb_trans++;
            fmpz_poly_solve_scale_2exp(P, k - (I->k));
            info->nb_homo++;
            k = I->k;
        }

        /*
          Compute the sign of P(1/2) ; 
          If # sign variations is 2  and  sign(P(0)) = sign(P(1)) = -sign(P(1/2)) 
          we have found two roots.
        */
        shalf = fmpz_poly_solve_sgn_eval_at_half(P);

        V = fmpz_poly_solve_Descartes_test(P, Q, shalf, &status, info);
        // printf("V: %ld \n", V);
        switch ( V ) {
        case 0:
            fmpz_bintvl_clear(I);
            flint_free(I);
            break;

        case 1:
            /* printf("Found 1 root!"); */
            fmpz_bintvl_new_root(roots, I, info);
            info->nb_roots++;

            fmpz_bintvl_clear(I);
            flint_free(I);
            break;

        case 2:
            if (status) {
                /* debug("status: %ld", status); */

                /* verb_info("Found 2 roots (1/2 hack) !"); */
                /* There is a root in (0, 1/2) and (1/2, 1) */
                fmpz_mul_2exp(I->c, I->c, 1);
                I->k++;
                fmpz_bintvl_new_root(roots, I, info);
                info->nb_roots++;

                fmpz_add_ui(I->c, I->c, 1);
                fmpz_bintvl_new_root(roots, I, info);
                info->nb_roots++;

                info->nb_half_opt++;

                fmpz_bintvl_clear(I);
                flint_free(I);
                break;
            }

        default:
            /* More than two roots. Split and continue */
            /* The left child */
            Il = flint_malloc(sizeof(fmpz_bintvl_t));
            fmpz_bintvl_init_set(Il, I);
            Il->k = I->k + 1;
            fmpz_mul_2exp(Il->c, Il->c, 1);

            /* The right child */
            Ir = flint_malloc(sizeof(fmpz_bintvl_t));
            fmpz_bintvl_init_set(Ir, I);
            Ir->k = I->k + 1;
            fmpz_set(Ir->c, Il->c);
            fmpz_add_ui(Ir->c, Ir->c, 1);

            
            FLINT_SLIST_PUSH(queue, Ir);
            FLINT_SLIST_PUSH(queue, Il);
        }
    }
    /* check_debug(FLINT_SLIST_EMPTY(&queue), "There are intervals that we DID NOT consider!"); */

    fmpz_clear(c);
    fmpz_clear(one);
    fmpz_poly_clear(P);
    fmpz_poly_clear(Q);
    return;

}



