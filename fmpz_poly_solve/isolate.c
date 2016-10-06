/*
  Copyright (C) 2016 Elias Tsigaridas

  This file is part of FLINT.

  FLINT is free software: you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License (LGPL) as published
  by the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"



fmpz_bintvl_t*
_fmpz_poly_solve_isolate_pos(const fmpz_poly_t F, slv_info_ptr info)
{
    fmpz_bintvl_t* roots;
 
    slong dg;

    dg = fmpz_poly_degree(F);
    info->dg = dg;
    info->t_dg = dg;
    
    roots =   (fmpz_bintvl_t*) flint_calloc(info->t_dg, sizeof(fmpz_bintvl_t));
    /* Isolate the roots using VCA */
    roots = fmpz_poly_solve_isol_vca_in_0_inf(F, info);

    /* fmpz_poly_solve_print_all_roots(stdout, roots, info->nb_roots); */
    return roots;
}


slong
fmpz_poly_solve_isolate(fmpz_dyadic_intvl_struct* vec_roots, slv_info_t info, const fmpz_poly_t P, int flag)
{
    slv_info_t t_info, p_info, n_info;
    
    fmpz_bintvl_t* roots;
    fmpz_bintvl_t* p_roots;
    fmpz_bintvl_t* n_roots;

    fmpz_poly_t A;
    slong dg, j;


    dg = fmpz_poly_degree(P);
    if ( dg <= 0 )
    {
        return 0;
    }

    if ( dg == 1)
    {
        printf("We  can handle degree 1 currently! \n");
    }
    
    fmpz_poly_init(A);
    fmpz_poly_set(A, P);

    switch ( flag )
    {
        
    case 1:
        /* Isolate the positive roots */
        slv_info_init(t_info);
        t_info->sign = 1;
       
        /* Isolate the roots using VCA */
        roots = _fmpz_poly_solve_isolate_pos(A, t_info);
        fmpz_poly_solve_adjust_bintvl_signs(P, roots, t_info);
        
        for(j = 0; j < t_info->nb_roots; j++)
        {
            fmpz_dyadic_intvl_set_bintvl(vec_roots+j, roots[j]);
        }

        if (info != NULL) slv_info_set(info, t_info);
    
        break;

    case -1:
        /* Isolate the positive roots */
        slv_info_init(t_info);
        t_info->sign = -1;

        fmpz_poly_negate_x(A);
        roots = _fmpz_poly_solve_isolate_pos(A, t_info);
        fmpz_poly_solve_adjust_bintvl_signs(P, roots, t_info);

        for (j = 0; j < t_info->nb_roots; j++)
        {
            fmpz_dyadic_intvl_set_bintvl(vec_roots+j, roots[j]);
        }

        if (info != NULL) slv_info_set(info, t_info);
        
        break;

    case 0:
    default:
        /* Isolate all the  roots */
        slv_info_init(p_info);
        p_info->dg = dg;
        p_info->t_dg = dg;
        p_info->sign = 1;

        /* Isolate the positive roots */
        p_roots = _fmpz_poly_solve_isolate_pos(A, p_info);
        fmpz_poly_solve_adjust_bintvl_signs(P, p_roots, p_info);

        slv_info_init(n_info);
        n_info->dg = dg;
        n_info->t_dg = dg;
        n_info->sign = -1;

        
        /* Isolate the negative roots */

        fmpz_poly_negate_x(A);
        n_roots = _fmpz_poly_solve_isolate_pos(A, n_info);
        fmpz_poly_solve_adjust_bintvl_signs(P, n_roots, n_info);
        
        for (j = 0; j < n_info->nb_roots; j++)
        {
            fmpz_dyadic_intvl_set_bintvl(vec_roots+j, n_roots[j]);
        }
        for (j = 0; j < p_info->nb_roots; j++)
        {
            fmpz_dyadic_intvl_set_bintvl(vec_roots + j + n_info->nb_roots, p_roots[j]);                
        }
 
        if (info != NULL) slv_info_merge(info, n_info, p_info);
    }        

    if (info == NULL) printf("RRR \n");
    fmpz_poly_clear(A); 
    return info->nb_roots;
}



