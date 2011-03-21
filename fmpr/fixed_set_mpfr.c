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


mp_size_t
fixed_set_mpfr(mp_ptr y, mpfr_srcptr t, mp_size_t prec)
{
    mp_size_t n;
    mpz_t z;
    mpfr_t u;

    mpz_init(z);
    mpfr_init2(u, mpfr_get_prec(t));

    mpfr_mul_2exp(u, t, prec * FLINT_BITS, GMP_RNDD);
    mpfr_get_z(z, u, GMP_RNDD);
    n = z->_mp_size;

    mpn_copyi(y, z->_mp_d, n);
    mpn_zero(y + n, prec - n);

    mpz_clear(z);
    mpfr_clear(u);

    return n;
}
