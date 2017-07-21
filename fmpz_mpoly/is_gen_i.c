/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

int fmpz_mpoly_is_gen_i(const fmpz_mpoly_t poly,
                                           slong k, const fmpz_mpoly_ctx_t ctx)
{
   int deg, rev;
   slong i;
   ulong * exps;

   TMP_INIT;

   if (poly->length != 1)
      return 0;

   if (!fmpz_is_one(poly->coeffs + 0))
      return 0;

   degrev_from_ord(deg, rev, ctx->ord);

   if (deg)
   {
      if ((poly->exps[0] >> (FLINT_BITS - poly->bits)) != 1)
         return 0;
   }

   TMP_START;

   exps = (ulong *) TMP_ALLOC((ctx->n - deg)*sizeof(ulong));

   fmpz_mpoly_get_monomial(exps, poly, 0, ctx);

   for (i = 0; i < k; i++)
   {
      if (exps[i] != 0)
         return 0;
   }

   if (exps[i++] != 1)
      return 0;

   for ( ; i < ctx->n - deg; i++)
   {
      if (exps[i] != 0)
         return 0;
   }   

   TMP_END;

   return 1;
}
