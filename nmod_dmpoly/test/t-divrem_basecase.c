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

#include <stdio.h>
#include <stdlib.h>
#include "nmod_dmpoly.h"

int main()
{
    flint_rand_t state;
    long iter;

    printf("divrem_basecase....");
    fflush(stdout);

    flint_randinit(state);

    /* Test exact division */
    for (iter = 0; iter < 10000; iter++)
    {
        nmod_dmpoly_t A, B, C, Q, R;
        mp_limb_t mod;
        long len;
        int vars;

        mod = n_randtest_prime(state, 0);
        len = 1 + n_randint(state, 10);
        vars = 1 + n_randint(state, 3);

        nmod_dmpoly_init(A, vars, mod);
        nmod_dmpoly_init(B, vars, mod);
        nmod_dmpoly_init(C, vars, mod);
        nmod_dmpoly_init(Q, vars, mod);
        nmod_dmpoly_init(R, vars, mod);

        do { nmod_dmpoly_randtest(B, state, len); }
        while (nmod_dmpoly_is_zero(B));

        nmod_dmpoly_randtest(C, state, len);
        nmod_dmpoly_mul(A, B, C);
        nmod_dmpoly_divrem_basecase(Q, R, A, B);

        if (!nmod_dmpoly_equal(Q, C) || !nmod_dmpoly_is_zero(R))
        {
            printf("FAIL\n");
            printf("A:\n"); nmod_dmpoly_print(A); printf("\n\n");
            printf("B:\n"); nmod_dmpoly_print(B); printf("\n\n");
            printf("C:\n"); nmod_dmpoly_print(C); printf("\n\n");
            printf("Q:\n"); nmod_dmpoly_print(Q); printf("\n\n");
            printf("R:\n"); nmod_dmpoly_print(R); printf("\n\n");
            printf("mod: %lu\n", mod);
            abort();
        }

        nmod_dmpoly_clear(A);
        nmod_dmpoly_clear(B);
        nmod_dmpoly_clear(C);
        nmod_dmpoly_clear(Q);
        nmod_dmpoly_clear(R);
    }

    /* Test A = Q*B + R */
    for (iter = 0; iter < 10000; iter++)
    {
        nmod_dmpoly_t A, B, C, Q, R;
        mp_limb_t mod;
        long len;
        int vars;

        mod = n_randtest_prime(state, 0);
        len = 1 + n_randint(state, 10);
        vars = 1 + n_randint(state, 3);

        nmod_dmpoly_init(A, vars, mod);
        nmod_dmpoly_init(B, vars, mod);
        nmod_dmpoly_init(C, vars, mod);
        nmod_dmpoly_init(Q, vars, mod);
        nmod_dmpoly_init(R, vars, mod);

        nmod_dmpoly_randtest(A, state, len);

        do { nmod_dmpoly_randtest(B, state, len); }
        while (nmod_dmpoly_is_zero(B));

        nmod_dmpoly_divrem_basecase(Q, R, A, B);

        nmod_dmpoly_mul(C, Q, B);
        nmod_dmpoly_add(C, C, R);

        if (!nmod_dmpoly_equal(A, C))
        {
            printf("FAIL\n");
            printf("A:\n"); nmod_dmpoly_print(A); printf("\n\n");
            printf("B:\n"); nmod_dmpoly_print(B); printf("\n\n");
            printf("Q:\n"); nmod_dmpoly_print(Q); printf("\n\n");
            printf("R:\n"); nmod_dmpoly_print(R); printf("\n\n");
            printf("Q*B + R:\n"); nmod_dmpoly_print(C); printf("\n\n");
            printf("mod: %lu\n", mod);
            abort();
        }

        nmod_dmpoly_clear(A);
        nmod_dmpoly_clear(B);
        nmod_dmpoly_clear(C);
        nmod_dmpoly_clear(Q);
        nmod_dmpoly_clear(R);
    }

    flint_randclear(state);
    printf("PASS\n");
    return 0;
}
