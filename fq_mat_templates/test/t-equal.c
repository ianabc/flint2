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

    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("equal....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        TEMPLATE(T, ctx_t) ctx;
        TEMPLATE(T, mat_t) A, B, C, D, E;
        TEMPLATE(T, t) x;
        slong m, n, j;

        TEMPLATE(T, ctx_randtest)(ctx, state);

        TEMPLATE(T, init)(x);

        m = n_randint(state, 20);
        n = n_randint(state, 20);

        TEMPLATE(T, mat_init)(A, m, n, ctx);
        TEMPLATE(T, mat_init)(B, m, n, ctx);
        TEMPLATE(T, mat_init)(C, m, n, ctx);
        TEMPLATE(T, mat_init)(D, m+1, n, ctx);
        TEMPLATE(T, mat_init)(E, m, n+1, ctx);

        if (TEMPLATE(T, mat_equal)(A, D) || TEMPLATE(T, mat_equal)(A, E))
        {
            printf("FAIL: different dimensions should not be equal\n");
            abort();
        }

        TEMPLATE(T, mat_randtest)(A, state, ctx);
        TEMPLATE(T, mat_set)(B, A);

        if (!TEMPLATE(T, mat_equal)(A, B))
        {
            printf("FAIL: copied matrices should be equal\n");
            abort();
        }

        if (m && n)
        {
            j = n_randint(state, m*n);
            TEMPLATE(T, one)(x, ctx);
            TEMPLATE(T, add)(A->entries + j, A->entries + j, x, ctx);

            if (TEMPLATE(T, mat_equal)(A, B))
            {
                printf("FAIL: modified matrices should not be equal\n");
                abort();
            }
        }

        TEMPLATE(T, clear)(x);
        TEMPLATE(T, mat_clear)(A);
        TEMPLATE(T, mat_clear)(B);
        TEMPLATE(T, mat_clear)(C);
        TEMPLATE(T, mat_clear)(D);
        TEMPLATE(T, mat_clear)(E);

        TEMPLATE(T, ctx_clear)(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return 0;
}


#endif