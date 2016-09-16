/*
    Copyright (C) 2016 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZ_MPOLY_H
#define FMPZ_MPOLY_H

#ifdef FMPZ_MPOLY_INLINES_C
#define FMPZ_MPOLY_INLINE FLINT_DLL
#else
#define FMPZ_MPOLY_INLINE static __inline__
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

#ifdef __cplusplus
 extern "C" {
#endif

/*  Type definitions *********************************************************/

typedef enum {
   ORD_LEX, ORD_REVLEX, ORD_DEGLEX, ORD_DEGREVLEX
} ordering_t;

typedef struct
{
   slong n;        /* number of elements in exponent vector (including deg) */
   ordering_t ord; /* polynomial ordering */
} fmpz_mpoly_ctx_struct;

typedef fmpz_mpoly_ctx_struct fmpz_mpoly_ctx_t[1];

typedef struct
{
   fmpz * coeffs; /* alloc fmpzs */
   ulong * exps;  
   slong alloc;
   slong length;
   slong bits;     /* number of bits per exponent */
} fmpz_mpoly_struct;

typedef fmpz_mpoly_struct fmpz_mpoly_t[1];

typedef struct fmpz_mpoly_heap_s
{
   ulong exp;
   slong i;
   slong j;
   struct fmpz_mpoly_heap_s * next;
} fmpz_mpoly_heap_s;

/* Macros ********************************************************************/

#define degrev_from_ord(deg, rev, ord)                                \
   (deg) = (rev) = 0;                                                 \
   do {                                                               \
      switch (ord)                                                    \
      {                                                               \
      case ORD_LEX:                                                   \
         break;                                                       \
      case ORD_DEGLEX:                                                \
         (deg) = 1; break;                                            \
      case ORD_REVLEX:                                                \
         (rev) = 1; break;                                            \
      case ORD_DEGREVLEX:                                             \
         (deg) = (rev) = 1; break;                                    \
      default:                                                        \
         flint_throw(FLINT_ERROR, "Invalid ordering in fmpz_mpoly");  \
      }                                                               \
   } while (0)

/* Context object ************************************************************/

FLINT_DLL void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, 
                                            slong nvars, const ordering_t ord);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx)
{
   /* nothing to be done at the moment */
}

/* Monomials *****************************************************************/

#if FLINT64

void _fmpz_mpoly_unpack_monomials_8to64(ulong * exps1, const ulong * exps2, 
                                                           slong n, slong len);

void _fmpz_mpoly_unpack_monomials_16to64(ulong * exps1, const ulong * exps2, 
                                                           slong n, slong len);

void _fmpz_mpoly_unpack_monomials_32to64(ulong * exps1, const ulong * exps2, 
                                                           slong n, slong len);

#endif

void _fmpz_mpoly_unpack_monomials_8to32(ulong * exps1, const ulong * exps2, 
                                                           slong n, slong len);

void _fmpz_mpoly_unpack_monomials_16to32(ulong * exps1, const ulong * exps2, 
                                                           slong n, slong len);

void _fmpz_mpoly_unpack_monomials_8to16(ulong * exps1, const ulong * exps2, 
                                                           slong n, slong len);

ulong * _fmpz_mpoly_unpack_monomials(slong bits1, const ulong * exps2, 
                                              slong bits2, slong n, slong len);

/*  Memory management ********************************************************/

FLINT_DLL void fmpz_mpoly_init(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_init2(fmpz_mpoly_t poly, slong alloc, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_realloc(fmpz_mpoly_t poly, slong alloc, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_fit_length(fmpz_mpoly_t poly, slong len, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_clear(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_normalise(fmpz_mpoly_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void _fmpz_mpoly_set_length(fmpz_mpoly_t poly, slong newlen, 
                                                   const fmpz_mpoly_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        slong i;
        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i); 
    }
    poly->length = newlen;
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_truncate(fmpz_mpoly_t poly, slong newlen, 
                                                   const fmpz_mpoly_ctx_t ctx)
{
    if (poly->length > newlen)
    {
        slong i;

        for (i = newlen; i < poly->length; i++)
            _fmpz_demote(poly->coeffs + i);

        poly->length = newlen;

        _fmpz_mpoly_normalise(poly, ctx);
    }  
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_fit_bits(fmpz_mpoly_t poly,
                                        slong bits, const fmpz_mpoly_ctx_t ctx)
{
   if (bits > poly->bits)
   {
      slong N = (bits*ctx->n - 1)/FLINT_BITS + 1;

      poly->exps = flint_realloc(poly->exps, N*poly->alloc*sizeof(ulong));
      poly->bits = bits;
   }   
}

/*  Basic manipulation *******************************************************/

FLINT_DLL ulong * _fmpz_mpoly_max_degrees1(const ulong * exps, slong len, 
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL ulong * fmpz_mpoly_max_degrees(const fmpz_mpoly_t poly, 
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_gen1(fmpz * poly, ulong * exps, slong i,
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_ui(fmpz_mpoly_t poly,
                                          ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_si(fmpz_mpoly_t poly,
                                          slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_fmpz(fmpz_mpoly_t poly,
                                   const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_equal_ui(const fmpz_mpoly_t poly,
                                          ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_equal_si(const fmpz_mpoly_t poly,
                                          slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_equal_fmpz(const fmpz_mpoly_t poly,
                                   const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
void fmpz_mpoly_swap(fmpz_mpoly_t poly1, 
                                fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
{
   fmpz * t1;
   ulong * t2;
   slong t3;

   t1 = poly1->coeffs;
   poly1->coeffs = poly2->coeffs;
   poly2->coeffs = t1;

   t2 = poly1->exps;
   poly1->exps = poly2->exps;
   poly2->exps = t2;

   t3 = poly1->length;
   poly1->length = poly2->length;
   poly2->length = t3;

   t3 = poly1->alloc;
   poly1->alloc = poly2->alloc;
   poly2->alloc = t3;

   t3 = poly1->bits;
   poly1->bits = poly2->bits;
   poly2->bits = t3;
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_zero(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
   _fmpz_mpoly_set_length(poly, 0, ctx);
}

FMPZ_MPOLY_INLINE
void fmpz_mpoly_one(fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set_ui(poly, UWORD(1), ctx);
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_is_zero(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
   return poly->length == 0;
}

FMPZ_MPOLY_INLINE
int fmpz_mpoly_is_one(const fmpz_mpoly_t poly, const fmpz_mpoly_ctx_t ctx)
{
   return fmpz_mpoly_equal_ui(poly, 1, ctx);
}

FLINT_DLL int fmpz_mpoly_is_gen(const fmpz_mpoly_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int fmpz_mpoly_is_gen_i(const fmpz_mpoly_t poly,
                                          slong k, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_coeff_fmpz(fmpz_t x,
                 const fmpz_mpoly_t poly, slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL ulong fmpz_mpoly_get_coeff_ui(const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong fmpz_mpoly_get_coeff_si(const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_fmpz(fmpz_mpoly_t poly, 
                          slong n, const fmpz_t x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_ui(fmpz_mpoly_t poly,
                                 slong n, ulong x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_coeff_si(fmpz_mpoly_t poly,
                                 slong n, slong x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_get_monomial(ulong * exps, const fmpz_mpoly_t poly, 
                                          slong n, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_set_monomial(fmpz_mpoly_t poly, 
                      slong n, const ulong * exps, const fmpz_mpoly_ctx_t ctx);

FMPZ_MPOLY_INLINE
int _fmpz_mpoly_monomial_is_zero(const ulong * exps, slong m)
{
   slong i;

   for (i = 0; i < m; i++)
   {
      if (exps[i] != 0)
         return 0;
   }

   return 1;
}

#define fmpz_mpoly_get_coeff_ptr(poly, n, ctx) \
    ((n) < (poly)->length ? (poly)->coeffs + (n) : NULL)

#define fmpz_mpoly_get_monomial_ptr(poly, n, ctx) \
    ((n) < (poly)->length ? (poly)->exps + (n)*(((ctx)->n - 1)/(FLINT_BITS/(poly)->bits) + 1) : NULL)

FLINT_DLL void _fmpz_mpoly_renormalise(fmpz_mpoly_t poly,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Set and negate */

FLINT_DLL void _fmpz_mpoly_set(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong m);

FLINT_DLL void fmpz_mpoly_set(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                   const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_neg(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong m);

FLINT_DLL void fmpz_mpoly_neg(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Comparison */

FLINT_DLL int _fmpz_mpoly_equal(fmpz * poly1, ulong * exps1,
                    const fmpz * poly2, const ulong * exps2, slong n, slong m);

FLINT_DLL int fmpz_mpoly_equal(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                   const fmpz_mpoly_ctx_t ctx);

/* Arithmetic ****************************************************************/

FLINT_DLL void fmpz_mpoly_add_ui(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_add_si(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_add_fmpz(fmpz_mpoly_t poly1,
         const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sub_ui(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sub_si(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void fmpz_mpoly_sub_fmpz(fmpz_mpoly_t poly1,
         const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_add1(fmpz * poly1, ulong * exps1,
                const fmpz * poly2, const ulong * exps2, slong len2,
                const fmpz * poly3, const ulong * exps3, slong len3,
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL void fmpz_mpoly_add(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                         const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL slong _fmpz_mpoly_sub1(fmpz * poly1, ulong * exps1,
                const fmpz * poly2, const ulong * exps2, slong len2,
                const fmpz * poly3, const ulong * exps3, slong len3,
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL void fmpz_mpoly_sub(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                         const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_mul_fmpz(fmpz * poly1, ulong * exps1,
 const fmpz * poly2, const ulong * exps2, slong len2, slong N, const fmpz_t c);

FLINT_DLL void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t poly1,
         const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_mul_si(fmpz * poly1, ulong * exps1,
        const fmpz * poly2, const ulong * exps2, slong len2, slong N, slong c);

FLINT_DLL void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_mul_ui(fmpz * poly1, ulong * exps1,
        const fmpz * poly2, const ulong * exps2, slong len2, slong N, ulong c);

FLINT_DLL void fmpz_mpoly_scalar_mul_ui(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_divexact_fmpz(fmpz * poly1, ulong * exps1,
 const fmpz * poly2, const ulong * exps2, slong len2, slong N, const fmpz_t c);

FLINT_DLL void fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_t poly1,
         const fmpz_mpoly_t poly2, const fmpz_t c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_divexact_si(fmpz * poly1, ulong * exps1,
        const fmpz * poly2, const ulong * exps2, slong len2, slong N, slong c);

FLINT_DLL void fmpz_mpoly_scalar_divexact_si(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, slong c, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL void _fmpz_mpoly_scalar_divexact_ui(fmpz * poly1, ulong * exps1,
        const fmpz * poly2, const ulong * exps2, slong len2, slong N, ulong c);

FLINT_DLL void fmpz_mpoly_scalar_divexact_ui(fmpz_mpoly_t poly1,
                const fmpz_mpoly_t poly2, ulong c, const fmpz_mpoly_ctx_t ctx);

/* Input/output **************************************************************/

FLINT_DLL void _exp_get_degrees1(ulong * expvec, ulong v, slong bits,
                                                    slong n, int deg, int rev);

FLINT_DLL char * _fmpz_mpoly_get_str_pretty1(const fmpz * poly,
                          const ulong * exps, slong len, const char ** x, 
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t poly,
                                  const char ** x, const fmpz_mpoly_ctx_t ctx);

FLINT_DLL int _fmpz_mpoly_fprint_pretty1(FILE * file, const fmpz * poly, 
                           const ulong * exps, slong len, const char ** x,
                                        slong bits, slong n, int deg, int rev);

FLINT_DLL int fmpz_mpoly_fprint_pretty(FILE * file, 
         const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
