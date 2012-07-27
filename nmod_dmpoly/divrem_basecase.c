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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "nmod_dmpoly.h"


/* TODO: aliasing? */

void __nmod_dmpoly_divrem_basecase(arr_ptr Q, arr_ptr R, arr_srcptr A,
    long lenA, arr_srcptr B, long lenB, int vars, nmod_t mod)
{
    if (vars == 1)
    {
        _nmod_poly_divrem((mp_ptr) Q, (mp_ptr) R, (mp_srcptr) A, lenA, (mp_srcptr) B, lenB, mod);
    }
    else
    {
        arr_struct tmp;
        arr_srcptr leadB = B + (lenB - 1);
        long i, iQ, iR;

        _nmod_dmpoly_init(&tmp, vars - 1);

        if (R != A)
        {
            for (i = 0; i < lenA; i++)
                _nmod_dmpoly_set(R + i, A + i, vars - 1);
        }

        for (iQ = lenA - lenB, iR = lenA - 1; iQ >= 0; iQ--, iR--)
        {
            if (((arr_ptr) R + iR)->length < leadB->length)
            {
                _nmod_dmpoly_zero(Q + iQ, vars - 1);
            }
            else
            {
                /* TODO: this could use div, not writing out the remainder */
                _nmod_dmpoly_divrem_basecase(Q + iQ, &tmp, R + iR, leadB, vars - 1, mod);

                for (i = 0; i < lenB; i++)
                {
                    /* R[iQ + i] -= B[i] * Q[iQ] */
                    _nmod_dmpoly_mul(&tmp, B + i, Q + iQ, vars - 1, mod);
                    _nmod_dmpoly_sub(R + iQ + i, R + iQ + i, &tmp, vars - 1, mod);
                }
            }
        }

        _nmod_dmpoly_clear(&tmp, vars - 1);
    }
}

void
_nmod_dmpoly_divrem_basecase(arr_ptr Q, arr_ptr R,
                                        arr_srcptr A, arr_srcptr B, int vars, nmod_t mod)
{
    long lenA, lenB, lenQ, lenR;

    lenA = A->length;
    lenB = B->length;

    if (lenA < lenB)
    {
        _nmod_dmpoly_set(R, A, vars);
        _nmod_dmpoly_zero(Q, vars);
        return;
    }

    if (lenA == 0)
    {
        _nmod_dmpoly_zero(Q, vars);
        _nmod_dmpoly_zero(R, vars);
        return;
    }

    if (vars == 1)
    {
        lenQ = lenA - lenB + 1;
        lenR = lenB - 1;
    }
    else
    {
        lenQ = lenA - lenB + 1;
        lenR = lenA;
    }

    _nmod_dmpoly_fit_length(Q, lenQ, vars);
    _nmod_dmpoly_fit_length(R, lenR, vars);

    __nmod_dmpoly_divrem_basecase(Q->coeffs, R->coeffs, A->coeffs, lenA,
        B->coeffs, lenB, vars, mod);

    Q->length = lenQ;
    R->length = lenR;
    _nmod_dmpoly_normalise(Q, vars);
    _nmod_dmpoly_normalise(R, vars);
}

void
nmod_dmpoly_divrem_basecase(nmod_dmpoly_t Q, nmod_dmpoly_t R, const nmod_dmpoly_t A,
                                    const nmod_dmpoly_t B)
{
    if (B->arr.length == 0)
    {
        printf("Exception: nmod_dmpoly_divrem: divide by zero\n");
        abort();
    }

    _nmod_dmpoly_divrem_basecase(&Q->arr, &R->arr, &A->arr, &B->arr, A->vars, A->mod);
}
