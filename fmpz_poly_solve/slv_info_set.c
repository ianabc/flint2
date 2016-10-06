/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
     (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void slv_info_set(slv_info_t info, const slv_info_t info1)
{
    info->dg = info1->dg;
    info->t_dg = info1->dg;
    info->bd = info1->bd;
    info->nb_roots = info1->nb_roots;

    info->nb_nodes = info1->nb_nodes;
    info->max_depth = info1->max_depth;
    info->nb_homo = info1->nb_homo;
    info->nb_trans = info1->nb_trans;
    info->nb_half_opt = info1->nb_half_opt;
    info->nb_pos_hack_1 = info1->nb_pos_hack_1;
    info->nb_pos_hack_2 = info1->nb_pos_hack_2;

	return;
}

