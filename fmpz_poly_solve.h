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

/* --- Linked list macros adapted from the BSD queue.h --- */

/*
 * Copyright (c) 1991, 1993
 *	The Regents of the University of California.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *	@(#)queue.h	8.5 (Berkeley) 8/20/94
 */

#ifdef QUEUE_MACRO_DEBUG
#define FLINT_Q_INVALIDATE(a) (a) = ((void *)-1)
#else
#define FLINT_Q_INVALIDATE(a)
#endif

/*
 * Singly-linked List definitions.
 */
#define FLINT_SLIST_HEAD(name, type) \
struct name { \
    struct type *slh_first; /* first element */ \
}

#define	FLINT_SLIST_HEAD_INITIALIZER(head) \
    { NULL }

#define FLINT_SLIST_ENTRY(type) \
struct { \
    struct type *sle_next; /* next element */ \
}

/*
 * Singly-linked List access methods.
 */
#define	FLINT_SLIST_FIRST(head) ((head)->slh_first)
#define	FLINT_SLIST_END(head) NULL
#define	FLINT_SLIST_EMPTY(head) (FLINT_SLIST_FIRST(head) == FLINT_SLIST_END(head))
#define	FLINT_SLIST_NEXT(elm, field) ((elm)->field.sle_next)

#define	FLINT_SLIST_FOREACH(var, head, field) \
    for((var) = FLINT_SLIST_FIRST(head); \
        (var) != FLINT_SLIST_END(head); \
        (var) = FLINT_SLIST_NEXT(var, field))

#define	FLINT_SLIST_FOREACH_PREVPTR(var, varp, head, field) \
    for ((varp) = &FLINT_SLIST_FIRST((head)); \
        ((var) = *(varp)) != FLINT_SLIST_END(head); \
        (varp) = &FLINT_SLIST_NEXT((var), field))

/*
 * Singly-linked List functions.
 */
#define	FLINT_SLIST_INIT(head) { \
    FLINT_SLIST_FIRST(head) = FLINT_SLIST_END(head); \
}

#define	FLINT_SLIST_INSERT_AFTER(slistelm, elm, field) do { \
    (elm)->field.sle_next = (slistelm)->field.sle_next; \
    (slistelm)->field.sle_next = (elm); \
} while (0)

#define	FLINT_SLIST_INSERT_HEAD(head, elm, field) do { \
    (elm)->field.sle_next = (head)->slh_first; \
    (head)->slh_first = (elm); \
} while (0)

#define	FLINT_SLIST_REMOVE_NEXT(head, elm, field) do { \
    (elm)->field.sle_next = (elm)->field.sle_next->field.sle_next; \
} while (0)

#define	FLINT_SLIST_REMOVE_HEAD(head, field) do { \
    (head)->slh_first = (head)->slh_first->field.sle_next; \
} while (0)

#define FLINT_SLIST_REMOVE(head, elm, type, field) do { \
    if ((head)->slh_first == (elm)) { \
        FLINT_SLIST_REMOVE_HEAD((head), field); \
    } else { \
        struct type *curelm = (head)->slh_first; \
                                    \
        while (curelm->field.sle_next != (elm)) \
            curelm = curelm->field.sle_next; \
        curelm->field.sle_next = \
            curelm->field.sle_next->field.sle_next; \
        FLINT_Q_INVALIDATE((elm)->field.sle_next); \
    }                               \
} while (0)

/* --- End of linked list macros adapted from the BSD queue.h --- */

struct slv_info
{
    slong 	max_depth;
    slong 	nb_nodes;
    slong	nb_homo;
    slong 	nb_trans;
    slong 	nb_half_opt;
    slong   nb_pos_hack_1;
    slong   nb_pos_hack_2;

    slong 	t_dg;       /* The degree of the input */
    slong 	dg;         /* The current degree. It might not be equal
                           to the t_dg if we find a rational root
                           in the process of isolation  */

    int     sign;       /* -1 if we solve for negative roots */
    slong   bd;         /* The root bound is 2^{bd} */
    slong   nb_roots;   /* The number of roots */
};

typedef struct slv_info slv_info_t[1];
typedef struct slv_info * slv_info_ptr;
typedef const struct slv_info * slv_info_srcptr;

FLINT_DLL void slv_info_init(slv_info_ptr info);
FLINT_DLL void slv_info_print(slv_info_srcptr info);

/* Struct for keeping an interval of a subdivision algorithm. */
struct fmpz_bintvl
{
    fmpz_t c; /* The interval is \f$(c/2^k, (c+1)/2^k )\f$ */
    slong k;
    unsigned int is_exact;
    int sgn_left;      /* sign on the left endpoint */
    unsigned int sign; /* are we solving in a negative interval? */
    FLINT_SLIST_ENTRY(fmpz_bintvl) pnext; /* pointer to the next entry in the list. */
};

/* List of Intervals */
FLINT_SLIST_HEAD(fmpz_lst_intvl_t, fmpz_bintvl) ;
typedef struct fmpz_lst_intvl_t fmpz_lst_bintvl;

/* Abbreviation for basic list operations */
#define FLINT_SLIST_PUSH(queue, I) FLINT_SLIST_INSERT_HEAD(&queue, I, pnext);

#define FLINT_SLIST_POP(queue, I)  \
    I = FLINT_SLIST_FIRST(&queue); \
    FLINT_SLIST_REMOVE_HEAD(&queue, pnext); \

typedef struct fmpz_bintvl          fmpz_bintvl_t[1];
typedef struct fmpz_bintvl *        fmpz_bintvl_ptr;
typedef const struct fmpz_bintvl *  fmpz_bintvl_srcptr;

FLINT_DLL void fmpz_bintvl_print(fmpz_bintvl_srcptr I);
FLINT_DLL void fmpz_bintvl_init(fmpz_bintvl_ptr a);
FLINT_DLL void fmpz_bintvl_init_set(fmpz_bintvl_ptr a, const fmpz_bintvl_srcptr b);
FLINT_DLL void fmpz_bintvl_clear(fmpz_bintvl_ptr a);

FLINT_DLL void fmpz_bintvl_get_mid(fmpq_t m, const fmpz_bintvl_srcptr I);
FLINT_DLL double fmpz_bintvl_get_mid_d(const fmpz_bintvl_srcptr I);

FLINT_DLL void fmpz_bintvl_new_root(fmpz_bintvl_t* vec_bintvl,
                     fmpz_bintvl_srcptr I,
                     slv_info_ptr info);

FLINT_DLL int fmpz_is_zero_a_root(fmpz_poly_t P,
                        fmpz_bintvl_t* vec_bintvl,
                        slv_info_ptr info);

   


   FLINT_DLL int fmpz_poly_solve_sgn_eval_at_half(const fmpz_poly_t P);
   FLINT_DLL int fmpz_poly_solve_sgn_eval_at_c(const fmpz_poly_t P, const fmpz_t c);
   
   FLINT_DLL int fmpz_poly_solve_sgn_eval_at_c_2exp(const fmpz_poly_t P,
                                                    const fmpz_t c,
                                                    slong k);
   
   FLINT_DLL slong fmpz_poly_solve_remove_content_2exp(fmpz_poly_t F);

   FLINT_DLL slong fmpz_poly_solve_scale_2exp(fmpz_poly_t F, slong k);

   
   FLINT_DLL slong fmpz_poly_solve_root_upper_bound_2exp(const fmpz_poly_t F);
   FLINT_DLL slong fmpz_poly_solve_root_lower_bound_2exp(const fmpz_poly_t F);

   FLINT_DLL slong fmpz_poly_solve_var(const fmpz_poly_t f);

   FLINT_DLL
   fmpz_bintvl_t* fmpz_poly_solve_isol_vca_in_0_inf(const fmpz_poly_t A, slv_info_ptr info);

     FLINT_DLL
     void fmpz_poly_solve_isol_vca_in_0_1(fmpz_poly_t FF, 
                                          fmpz_bintvl_t* roots, 
                                          slv_info_ptr info);
     

   
#ifdef __cplusplus
}
#endif

#endif

