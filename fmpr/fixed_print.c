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


void fixed_print(mp_srcptr x, long limbs, long frac_limbs)
{
    mpz_t y;
    mpfr_t f;
    mpz_init(y);
    mpz_import(y, limbs, -1, sizeof(mp_limb_t), 0, 0, x);
    mpfr_init2(f, limbs*FLINT_BITS);
    mpfr_set_z_2exp(f, y, -frac_limbs*FLINT_BITS, GMP_RNDN);
    mpfr_printf("%.75Rf\n", f);
    mpfr_clear(f);
    mpz_clear(y);
}
