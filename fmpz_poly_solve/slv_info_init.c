/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_solve.h"

void slv_info_init(slv_info_ptr info)
{
	info->max_depth = 0;
	info->nb_nodes = 0;
	info->nb_homo = 0;
	info->nb_trans = 0;
	info->nb_half_opt = 0;
    info->nb_pos_hack_1 = 0;
    info->nb_pos_hack_2 = 0;
    
	info->t_dg = 0;
    info->dg = 0;
	info->bd = -1;
    info->nb_roots = 0;

	return;
}

