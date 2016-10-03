/*
    Copyright (C) 2016 Elias Tsigaridas

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_POLY_SOLVE_H
#define FMPZ_POLY_SOLVE_H

#ifdef FMPZ_POLY_SOLVE_INLINES_C
#define FMPZ_POLY_SOLVE_INLINE FLINT_DLL
#else
#define FMPZ_POLY_SOLVE_INLINE static __inline__
#endif

#undef ulong
#define ulong ulongxx /* interferes with system includes */
#include <stdio.h>
#undef ulong
#include <gmp.h>
#define ulong mp_limb_t

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

   FLINT_DLL int fmpz_poly_solve_sgn_eval_at_half(const fmpz_poly_t P);
   FLINT_DLL int fmpz_poly_solve_sgn_eval_at_c(const fmpz_poly_t P, const fmpz_t c);
   
   FLINT_DLL int fmpz_poly_solve_sgn_eval_at_c_2exp(const fmpz_poly_t P,
                                                    const fmpz_t c,
                                                    long k);


   
#if 0
   FLINT_DLL void fmpz_poly_solve_init(fmpz_poly_solve_t fac);

FLINT_DLL void fmpz_poly_solve_init2(fmpz_poly_solve_t fac, slong alloc);

FLINT_DLL void fmpz_poly_solve_realloc(fmpz_poly_solve_t fac, slong alloc);

FLINT_DLL void fmpz_poly_solve_fit_length(fmpz_poly_solve_t fac, slong len);

FLINT_DLL void fmpz_poly_solve_clear(fmpz_poly_solve_t fac);

FLINT_DLL void fmpz_poly_solve_set(fmpz_poly_solve_t res, const fmpz_poly_solve_t fac);

FLINT_DLL void fmpz_poly_solve_insert(fmpz_poly_solve_t fac, 
                             const fmpz_poly_t p, slong exp);

FLINT_DLL void fmpz_poly_solve_concat(fmpz_poly_solve_t res, 
                             const fmpz_poly_solve_t fac);

FLINT_DLL void fmpz_poly_solve_print(const fmpz_poly_solve_t fac);

FLINT_DLL void fmpz_poly_solve_zassenhaus_recombination(fmpz_poly_solve_t final_fac, 
	const fmpz_poly_solve_t lifted_fac, 
    const fmpz_poly_t F, const fmpz_t P, slong exp);
    
FLINT_DLL void fmpz_poly_solve_squarefree(fmpz_poly_solve_t fac, const fmpz_poly_t F);

FLINT_DLL void _fmpz_poly_solve_zassenhaus(fmpz_poly_solve_t final_fac, 
								  slong exp, const fmpz_poly_t f, slong cutoff);

FLINT_DLL void fmpz_poly_solve_zassenhaus(fmpz_poly_solve_t fac, const fmpz_poly_t G);

#endif /* if 0 */
   
#ifdef __cplusplus
}
#endif

#endif

