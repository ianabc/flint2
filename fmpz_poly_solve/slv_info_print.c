/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void slv_info_print(slv_info_srcptr info)
{
    printf("/*------------------------------------------------------------*/\n");
    printf("  Statistics: \n");

    printf( "  #degree : %3lu \n", info->dg);
    printf( "  #bound  : %3lu \n", info->bd);
    
    printf( "  #roots  : %3lu \n", info->nb_roots);
    printf( "  #nodes  : %3lu \n", info->nb_nodes);
    printf( "  #depth  : %3lu \n", info->max_depth);

    printf( "  #trans  : %3lu \n", info->nb_trans);
    printf( "  #homo   : %3lu \n\n", info->nb_homo);

    printf( "  #pos_h_1  : %3lu \n", info->nb_pos_hack_1);
    printf( "  #pos_h_2  : %3lu \n", info->nb_pos_hack_2);
    printf( "  #half_h   : %3lu \n", info->nb_half_opt);
    
    printf("/*------------------------------------------------------------*/\n");
}

