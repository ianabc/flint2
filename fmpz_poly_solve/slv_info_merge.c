/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
     (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void slv_info_merge(slv_info_t info, const slv_info_t info1, const slv_info_t info2)
{
    info->dg = FLINT_MAX(info1->dg, info2->dg);
    info->t_dg = info1->dg;
    info->bd = FLINT_MAX(info1->bd, info2->bd);
    info->nb_roots = info1->nb_roots + info2->nb_roots;

    info->nb_nodes = info1->nb_nodes + info2->nb_nodes;
    info->max_depth = FLINT_MAX(info1->max_depth, info2->max_depth);
    info->nb_homo = info1->nb_homo + info2->nb_homo;
    info->nb_trans = info1->nb_trans + info2->nb_trans;
    info->nb_half_opt = info1->nb_half_opt + info2->nb_half_opt;
    info->nb_pos_hack_1 = info1->nb_pos_hack_1 + info2->nb_pos_hack_1;
    info->nb_pos_hack_2 = info1->nb_pos_hack_2 + info2->nb_pos_hack_2;

	return;
}

