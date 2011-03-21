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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <mpir.h>
#include <mpfr.h>
#include "fmpr.h"


void fixed_exp_mpfr(mp_ptr y, mp_srcptr x, mp_size_t limbs, long tol_bits)
{
    mpfr_t t;
    mpfr_init2(t, limbs * FLINT_BITS + 2);

    fixed_get_mpfr(t, x, limbs, limbs);
    mpfr_exp(t, t, GMP_RNDD);
    fixed_set_mpfr(y, t, limbs);

    mpfr_clear(t);
}
