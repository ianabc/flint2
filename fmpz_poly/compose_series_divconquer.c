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

    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/*
    See fmpz_poly/compose_divconquer.c
 */

void
_fmpz_poly_compose_series_divconquer(fmpz * res, const fmpz * poly1, long len1, 
                                                 const fmpz * poly2, long len2, 
                                                 long N)
{
    long i, j, k, n;
    long *hlen, alloc, powlen, v_alloc;
    fmpz *v, **h, *pow, *temp;

    if (len1 <= 2 || len2 == 1 || N == 1)
    {
        _fmpz_poly_compose_series_horner(res, poly1, len1, poly2, len2, N);
        return;
    }

    /* Initialisation */

    hlen = (long *) flint_malloc(((len1 + 1) / 2) * sizeof(long));

    for (k = 1; (2 << k) < len1; k++) ;

    hlen[0] = hlen[1] = FLINT_MIN(N, ((1 << k) - 1) * (len2 - 1) + 1);
    for (i = k - 1; i > 0; i--)
    {
        long hi = (len1 + (1 << i) - 1) / (1 << i);
        long t  = FLINT_MIN(N, ((1 << i) - 1) * (len2 - 1) + 1);
        for (n = (hi + 1) / 2; n < hi; n++)
            hlen[n] = t;
    }
    powlen = FLINT_MIN(N, (1 << k) * (len2 - 1) + 1);

    alloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc += hlen[i];

    v_alloc = alloc + 2 * powlen;
    v = _fmpz_vec_init(v_alloc);
    h = (fmpz **) flint_malloc(((len1 + 1) / 2) * sizeof(fmpz *));
    h[0] = v;
    for (i = 0; i < (len1 - 1) / 2; i++)
    {
        h[i + 1] = h[i] + hlen[i];
        hlen[i]  = 0;
    }
    hlen[(len1 - 1) / 2] = 0;
    pow  = v + alloc;
    temp = pow + powlen;

    /* Let's start the actual work */
    
    for (i = 0, j = 0; i < len1 / 2; i++, j += 2)
    {
        if (poly1[j + 1] != 0L || 1)
        {
            _fmpz_vec_scalar_mul_fmpz(h[i], poly2, len2, poly1 + j + 1);
            fmpz_add(h[i], h[i], poly1 + j);
            hlen[i] = len2;
        }
        else if (poly1[j] != 0L)
        {
            fmpz_set(h[i], poly1 + j);
            hlen[i] = 1;
        }
    }

    if ((len1 & 1L))
    {
        if (poly1[j] != 0)
        {
            fmpz_set(h[i], poly1 + j);
            hlen[i] = 1;
        }
    }

    powlen = FLINT_MIN(N, 2 * len2 - 1);
    _fmpz_poly_mullow(pow, poly2, len2, poly2, len2, powlen);


    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        if (hlen[1] > 0)
        {
            long templen = FLINT_MIN(N, powlen + hlen[1] - 1);

            _fmpz_poly_mullow(temp, pow, powlen, h[1], hlen[1], templen);
            _fmpz_poly_add(h[0], temp, templen, h[0], hlen[0]);
            hlen[0] = FLINT_MAX(hlen[0], templen);
        }
        
        for (i = 1; i < n / 2; i++)
        {
            if (hlen[2*i + 1] > 0)
            {
                hlen[i] = FLINT_MIN(N, hlen[2*i + 1] + powlen - 1);
                _fmpz_poly_mullow(h[i], pow, powlen, h[2*i + 1], hlen[2*i + 1],
                    hlen[i]);
            }
            else
            {
                hlen[i] = 0;
            }

            _fmpz_poly_add(h[i], h[i], hlen[i], h[2*i], hlen[2*i]);
            hlen[i] = FLINT_MAX(hlen[i], hlen[2*i]);
        }

        if ((n & 1L))
        {
            hlen[i] = FLINT_MIN(N, hlen[2*i]);
            _fmpz_vec_set(h[i], h[2*i], hlen[i]);
        }
        
        _fmpz_poly_mullow(temp, pow, powlen, pow, powlen, 
                          FLINT_MIN(N, 2 * powlen - 1));
        powlen = FLINT_MIN(N, 2 * powlen - 1);
        {
            fmpz * t = pow;
            pow      = temp;
            temp     = t;
        }
    }

    _fmpz_poly_mullow(res, pow, powlen, h[1], hlen[1], 
                      FLINT_MIN(N, powlen + hlen[1] - 1));
    _fmpz_vec_add(res, res, h[0], hlen[0]);

    _fmpz_vec_clear(v, v_alloc);
    flint_free(h);
    flint_free(hlen);
}

void 
fmpz_poly_compose_series_divconquer(fmpz_poly_t res, 
    const fmpz_poly_t poly1, const fmpz_poly_t poly2, long N)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long lenr;

    if (len2 != 0 && !fmpz_is_zero(poly2->coeffs))
    {
        printf("exception: fmpz_poly_compose_series_divconquer: inner polynomial "
                "must have zero constant term\n");
        abort();
    }

    if (len1 == 0 || N == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        fmpz_poly_set_fmpz(res, poly1->coeffs);
        return;
    }
    
    lenr = FLINT_MIN(N, (len1 - 1) * (len2 - 1) + 1);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if (res != poly1 && res != poly2)
    {
        fmpz_poly_fit_length(res, lenr);
        _fmpz_poly_compose_series_divconquer(res->coeffs, 
            poly1->coeffs, len1, poly2->coeffs, len2, lenr);
        _fmpz_poly_set_length(res, lenr);
        _fmpz_poly_normalise(res);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, lenr);
        _fmpz_poly_compose_series_divconquer(t->coeffs, 
            poly1->coeffs, len1, poly2->coeffs, len2, lenr);
        _fmpz_poly_set_length(t, lenr);
        _fmpz_poly_normalise(t);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
