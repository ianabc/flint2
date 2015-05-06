/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2009 William Hart
    Copyright (C) 2015 Tommy Hofmann

******************************************************************************/

#define NMOD_MAT_INLINES_C

#define ulong ulongxx /* interferes with system includes */
#include <stdlib.h>
#undef ulong
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "nmod_mat.h"

mp_limb_t _nmod_mat_get_entry(nmod_mat_t mat, slong i, slong j)
{
  mp_limb_t x;
  x = nmod_mat_entry(mat, i, j);
  return x;
}

void _nmod_mat_set_entry(nmod_mat_t mat, slong i, slong j, mp_limb_t x)
{
  nmod_mat_entry(mat, i, j) = x;
}
