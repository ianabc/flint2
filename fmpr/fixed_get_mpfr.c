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


void
fixed_get_mpfr(mpfr_t dest, mp_srcptr src, mp_size_t size, mp_size_t prec)
{
    long i;
    mpz_t z;

    for (i = size; i > 0 && src[i-1] == 0; i--);

    z->_mp_d = (mp_ptr) src;
    z->_mp_alloc = size;
    z->_mp_size = i;

    mpfr_set_z_2exp(dest, z, -FLINT_BITS * prec, GMP_RNDD);
}
