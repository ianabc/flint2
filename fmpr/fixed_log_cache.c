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

static const mp_limb_t atanh_denom_1[] = {
    6332659870762850625UL
};

static const mp_limb_t atanh_denom_2[] = {
    885821206052908125UL,  9807130727058790125UL
};

static const mp_limb_t atanh_denom_3[] = {
    23481740411754625UL,  168059789309637225UL,  606997343490162625UL
};

static const mp_limb_t atanh_denom_4[] = {
    5555477684381625UL,  31974528653557875UL,  96229744645159125UL,
    232032717002861775UL
};

static const mp_limb_t atanh_coeffs_1[] = {
    6332659870762850625UL, 2110886623587616875UL, 1266531974152570125UL,
    904665695823264375UL, 703628874529205625UL, 575696351887531875UL,
    487127682366373125UL, 422177324717523375UL, 372509404162520625UL,
    333297887934886875UL, 301555231941088125UL, 275333037859254375UL,
    253306394830514025UL, 234542958176401875UL, 218367581750443125UL,
    204279350669769375UL, 191898783962510625UL,
};

static const mp_limb_t atanh_coeffs_2[] = {
    885821206052908125UL, 3269043575686263375UL, 177164241210581625UL,
    1401018675294112875UL, 98424578450323125UL, 891557338823526375UL,
    68140092773300625UL, 653808715137252675UL, 52107129767818125UL,
    516164775108357375UL, 42181962192995625UL, 426396988132990875UL,
    35432848242116325UL, 363227063965140375UL, 30545558829410625UL,
    316359055711573875UL, 26843066850088125UL, 280203735058822575UL,
    23941113677105625UL, 251464890437404875UL, 21605395269583125UL,
    228072807606018375UL, 19684915690064625UL, 208662355894867875UL,
    18077983796998125UL, 192296680922721375UL, 16713607661375625UL,
    178311467764705275UL,
};

static const mp_limb_t atanh_coeffs_3[] = {
    23481740411754625UL, 56019929769879075UL, 121399468698032525UL,
    3354534344536375UL, 18673309923293025UL, 55181576680923875UL,
    1806287723981125UL, 11203985953975815UL, 35705726087656625UL,
    1235881074302875UL, 8002847109982725UL, 26391188847398375UL,
    939269616470185UL, 6224436641097675UL, 20930942878971125UL,
    757475497153375UL, 5092720888170825UL, 17342781242576075UL,
    634641632750125UL, 4309225366913775UL, 14804813255857625UL,
    546086986319875UL, 3734661984658605UL, 12914837095535375UL,
    479219192076625UL, 3295289986463475UL, 11452780065852125UL,
    426940734759175UL, 2948417356309425UL, 10288090567629875UL,
    384946564127125UL, 2667615703327575UL, 9338420669079425UL,
    350473737488875UL, 2435649120429525UL, 8549258359016375UL,
};

static const mp_limb_t atanh_coeffs_4[] = {
    5555477684381625UL, 10658176217852625UL, 19245948929031825UL,
    33147531000408825UL, 617275298264625UL, 2906775332141625UL,
    7402288049627625UL, 15468847800190785UL, 326792804963625UL,
    1682869929134625UL, 4582368792626625UL, 10088379000124425UL,
    222219107375265UL, 1184241801983625UL, 3318267056729625UL,
    7484926354931025UL, 168347808617625UL, 913557961530225UL,
    2600803909328625UL, 5949556846227225UL, 135499455716625UL,
    743593689617625UL, 2138438769892425UL, 4936866319209825UL,
    113377095599625UL, 626951542226625UL, 1815655559342625UL,
    4218776672779305UL, 97464520778625UL, 541941163619625UL,
    1577536797461625UL, 3683059000045425UL, 85468887452025UL,
    477231770948625UL, 1394633980364625UL, 3268066436660025UL,
    76102434032625UL, 426327048714105UL, 1249736943443625UL,
    2937123000036225UL, 68586144251625UL, 385235284982625UL,
    1132114642884225UL, 2667042724170825UL,
};

static const mp_limb_t * atanh_coeffs[] =
{
    NULL,
    atanh_coeffs_1,
    atanh_coeffs_2,
    atanh_coeffs_3,
    atanh_coeffs_4
};

static const mp_limb_t * atanh_denom[] =
{
    NULL,
    atanh_denom_1,
    atanh_denom_2,
    atanh_denom_3,
    atanh_denom_4
};

#define ATANH_MAXTERMS_1 (sizeof(atanh_coeffs_1) / sizeof(mp_limb_t))
#define ATANH_MAXTERMS_2 (sizeof(atanh_coeffs_2) / sizeof(mp_limb_t))
#define ATANH_MAXTERMS_3 (sizeof(atanh_coeffs_3) / sizeof(mp_limb_t))
#define ATANH_MAXTERMS_4 (sizeof(atanh_coeffs_4) / sizeof(mp_limb_t))


static __inline__ long
log_needed_terms(long reduced, long tol_bits)
{
    return tol_bits / (2*reduced) + 1;
}


void fixed_log_cache(mp_ptr y, mp_srcptr x, mp_size_t limbs, long tol_bits)
{
    mp_limb_t top;
    long terms;
    int sums;
    int i1;

    if (log_cache_initialised == 0)
        log_cache_init();

    /*
    log(1+a+b) = log(1+a) + 2*atanh(b/(2+b+2a)) = log(1+a) + 2*s*F(s^2)

    where s = b/(2+b+2a), F(x) = 1 + x/3 + x^2/5 + x^4/7 + ...
    */
    top = x[limbs - 1];
    i1 = top >> (FLINT_BITS - LOG_CACHE_BITS);
    terms = log_needed_terms(LOG_CACHE_BITS, tol_bits);

    if      (terms <= ATANH_MAXTERMS_1) sums = 1;
    else if (terms <= ATANH_MAXTERMS_2) sums = 2;
    else if (terms <= ATANH_MAXTERMS_3) sums = 3;
    else if (terms <= ATANH_MAXTERMS_4) sums = 4;

    if (terms <= ATANH_MAXTERMS_4)
    {
        mp_limb_t a[LOG_CACHE_PREC_LIMBS*2 + 2];
        mp_limb_t b[LOG_CACHE_PREC_LIMBS*2 + 2];
        mp_limb_t s[LOG_CACHE_PREC_LIMBS*2 + 2];
        mp_limb_t t[LOG_CACHE_PREC_LIMBS*2 + 2];

        /* b = x - a, shifted for division */
        mpn_zero(b, limbs);
        mpn_copyi((b + limbs), x, limbs - 1);
        (b + limbs)[limbs - 1] = (top << LOG_CACHE_BITS) >> LOG_CACHE_BITS;

        /* t = b + 2 + 2a */
        mpn_copyi(t, b + limbs, limbs);
        t[limbs - 1] |= (((mp_limb_t) i1) << (FLINT_BITS - LOG_CACHE_BITS + 1));
        t[limbs] = 2 + (i1 >> (LOG_CACHE_BITS - 1));

        /* s = b / t */
        mpn_tdiv_q(s, b, 2*limbs, t, limbs + 1);
        mpn_sqr_basecase(t, s, limbs);

        fixed_eval_series_1(a, atanh_coeffs[sums], atanh_denom[sums],
            sums != 1, t + limbs, limbs, terms, FLINT_MAX(sums,2));

        mpn_mul_basecase(t, a, limbs + 1, s, limbs);
        mpn_lshift(t + limbs, t + limbs, limbs, 1);

        mpn_add_n(y, t + limbs, log_cache[i1]
                                    + (LOG_CACHE_PREC_LIMBS - limbs), limbs);
        return;
    }

    printf("fixed_log: not implemented: %ld terms\n", terms);
    abort();
}
