#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <assert.h>
#include <mpfr.h>
#include <limits.h>

#include "stringify.h"
#include "subr.h"
#include "mpfr_pi_generic.h"


/*
 * Compute PI using MPFR abitrary precision floating point library to N digits,
 * using Srinivasa Ramanujan's formula from 1910.
 *
 * Copyright (C) Fio Cattaneo <fio@cattaneo.us>, All Rights Reserved.
 *
 * This code is distributed under dual BSD/GPLv2 open source license.
 *
 */

/*
 * MPFR arbitrary precision floating point library docs:
 *
 * http://cs.swan.ac.uk/~csoliver/ok-sat-library/internet_html/doc/doc/Mpfr/3.0.0/mpfr.html/index.html#Top
 */

/*
 * Srinivasa Ramanujan 1910 formula, optimized version.
 *
 * More info on Ramanujan PI formulas:
 * https://en.wikipedia.org/wiki/Srinivasa_Ramanujan
 * https://en.wikipedia.org/wiki/Approximations_of_%CF%80
 * https://en.wikipedia.org/wiki/Ramanujan%E2%80%93Sato_series
 *
 *
 * Standard Formula
 * =================================================================
 *
 * 1/PI = CMULT * SUM(k, 0..infinity) TERM(k)
 * 
 * CMULT = (2 * sqrt(2)) / 9801			      # constant
 *                                                    # 9801 = 99^2
 *
 * TERM(k) = [ (4 * k)! * (1103 + 26390 * k) ] /      # dividend
 *           [ ((k!) ^ 4) * (396 ^ (4 * k)) ]         # divisor
 *                                                    # 396 = 99 * 4
 * 
 * Optimized Formula
 * =================================================================
 *
 * rewrite equations with some sub-expressions calculated based on the previous iteration.
 *
 * observe that K! = factorial(k) is the same as:
 *      factorial(0) = 1
 *      factorial(k) = factorial(k - 1) * k
 *      
 * =================================================================
 *
 * thus, rewrite
 *
 * TERM(k) = [ (4 * k)! * (1103 + 26390 * k) ] /      # dividend
 *           [ ((k!) ^ 4) * (396 ^ (4 * k)) ]         # divisor
 *                                                    # 396 = 99 * 4
 *
 * as
 *
 * TERM(k) = [ (4k)! * (1103 + 26390 * k) ] /         # dividend
 *           [ (FACT(k) ^ 4) * (396 ^ (4k)) ]         # divisor
 *                                                    # 396 = 99 * 4
 * where
 *
 * next(k) = k+1
 * next(4k) = 4k+4
 *
 * FACT(k == 0): 1                                                   
 * FACT(k != 0): FACT(k - 1) * k
 *
 * =================================================================
 *
 * similarly, rewrite
 *
 * TERM(k) = [ (4k)! * (1103 + 26390 * k) ] /         # dividend
 *           [ (FACT(k) ^ 4) * (396 ^ (4k)) ]         # divisor
 *                                                    # 396 = 99 * 4
 *
 * as
 *
 * TERM(k) = [ FACT4(4k) * (1103 + 26390 * k) ] /     # dividend
 *           [ (FACT(k) ^ 4) * (396 ^ (4k)) ]         # divisor
 *                                                    # 396 = 99 * 4
 *
 * where
 *
 * prev(k) = k-1
 * prev(4k) = 4k-4
 * next(k) = k+1
 * next(4k) = 4k+4
 *
 * FACT(k == 0): 1                                                   
 * FACT(k != 0): FACT(k - 1) * k
 *
 * FACT4(4k == 0): 1
 * FACT4(4k != 0): FACT(4k-4)) * 4k * (4k-1) * (4k-2) * (4k-3)
 *
 */

static const char *pi_impl_ramanujan_1910_opt_get_name(void)
{
	return "Ramanujan 1910 Formula (optimized)";
}

static void pi_impl_ramanujan_1910_opt_deinitialize(struct mpfr_pi_impl *impl);
static const mpfr_t * pi_impl_ramanujan_1910_opt_compute_next_term(struct mpfr_pi_impl *impl, unsigned long *out_k);
static const void pi_impl_ramanujan_1910_opt_get_value(const mpfr_t *partial_sum, mpfr_t *result_out);

/* actual implementation struct for this algorithm */
struct __mpfr_pi_impl {
	/* generic part */
	struct mpfr_pi_impl g;
	/* private part */
	unsigned long start_k; /* starting k, just for record */
	unsigned long curr_k; /* current k */
	unsigned long curr_4k; /* current 4k */
	/* various state variables needed */
	/*
	 * temp variables reused at each iteration
	 */
	mpfr_t curr_fact_k;
	mpfr_t curr_fact_4k;
	mpfr_t term_dividend;
	mpfr_t term_divisor;
	mpfr_t term;
	mpfr_t t0;
	/* sum of all current terms */
	mpfr_t term_sum;
	/* cmult constant, computed at initialization */
	mpfr_t cmult;
	/* actual pi, computed on demand or every now and then */
	mpfr_t pi;
};

#if 0
#define DIGITS_TO_K(d)	(((d) / 8L) + 1L)		/* number of iterations to get "d" digits */
/* arg reused, must pass l-value */
#define K_TO_DIGITS(k)	((k) >= SLACK_K ? ((k) * 8L) - SLACK_K : 0L)
#endif

#define __SAFE_LONG_MAX		(LONG_MAX / 2L)
#define __SAFE_ULONG_MAX	(ULONG_MAX / 2UL)

struct mpfr_pi_impl *pi_impl_ramanujan_1910_opt_initialize(unsigned long digits, const unsigned long k_start)
{
	struct __mpfr_pi_impl *__impl = malloc(sizeof (struct __mpfr_pi_impl));
	assert(__impl != NULL);
	__impl->g.f_impl_get_name = pi_impl_ramanujan_1910_opt_get_name;
	__impl->g.f_initialize = pi_impl_ramanujan_1910_opt_initialize;
	__impl->g.f_deinitialize = pi_impl_ramanujan_1910_opt_deinitialize;
	__impl->g.f_pi_compute_next_term = pi_impl_ramanujan_1910_opt_compute_next_term;
	__impl->g.f_pi_get_value = pi_impl_ramanujan_1910_opt_get_value;

	/*
	 * make number of iterations for requested digits will not make the divisor in
	 * the term sum overflow, not for (4 * k), not for the expression below:
	 *
	 * mpfr_mul_ui(__impl->term_dividend, __impl->curr_fact_4k, (1103UL + 26390UL * k), CFG_MPFR_RND);
	 */
	assert(k_start < (__SAFE_ULONG_MAX / 4UL));
	assert(k_start < (__SAFE_ULONG_MAX / 26390UL));
	__impl->start_k = k_start;
	__impl->curr_k = __impl->start_k;
	__impl->curr_4k = __impl->start_k * 4UL;
	/* various state variables needed */
	/*
	 * FIXME: need to dynamically calculate desired precision by using "digits" input arg
	 */
	mpfr_init2(__impl->curr_fact_k, CFG_MPFR_PREC);
	mpfr_init2(__impl->curr_fact_4k, CFG_MPFR_PREC);
	mpfr_init2(__impl->term_dividend, CFG_MPFR_PREC);
	mpfr_init2(__impl->term_divisor, CFG_MPFR_PREC);
	mpfr_init2(__impl->term, CFG_MPFR_PREC);
	mpfr_init2(__impl->term_sum, CFG_MPFR_PREC);
	mpfr_init2(__impl->cmult, CFG_MPFR_PREC);
	mpfr_init2(__impl->pi, CFG_MPFR_PREC);
	mpfr_init2(__impl->t0, CFG_MPFR_PREC);

	/* set term_sum */
	mpfr_set_ui(__impl->term_sum, 0UL, CFG_MPFR_RND);

	/* set FACT(k) */
	mpfr_fac_ui(__impl->curr_fact_k, __impl->start_k, CFG_MPFR_RND);
	/* set FACT4(k) = FACT4(k) = FACT(4 * k) */
	mpfr_fac_ui(__impl->curr_fact_4k, (__impl->start_k * 4UL), CFG_MPFR_RND);

	return (struct mpfr_pi_impl *)__impl;
}

static void pi_impl_ramanujan_1910_opt_deinitialize(struct mpfr_pi_impl *impl)
{
	struct __mpfr_pi_impl *__impl = (struct __mpfr_pi_impl *)impl;

	mpfr_clear(__impl->curr_fact_k);
	mpfr_clear(__impl->curr_fact_4k);
	mpfr_clear(__impl->term_dividend);
	mpfr_clear(__impl->term_divisor);
	mpfr_clear(__impl->term);
	mpfr_clear(__impl->term_sum);
	mpfr_clear(__impl->cmult);
	mpfr_clear(__impl->pi);
	mpfr_clear(__impl->t0);
	free(__impl);
}

static const mpfr_t * pi_impl_ramanujan_1910_opt_compute_next_term(struct mpfr_pi_impl *impl, unsigned long *out_k)
{
	struct __mpfr_pi_impl *__impl = (struct __mpfr_pi_impl *)impl;
	const unsigned long k = __impl->curr_k;
	const unsigned long _4k = __impl->curr_4k;
	int ret;

	/*
	 * calculate dividend
	 *
	 * [ FACT4(4k) * (1103 + 26390 * k) ]
	 *
	 */

	/*
	 * this being 64-bits, we can safely compute the expression with normal arithmetic
	 *
	 * (1) mpfr_set_ui(__impl->t0, (1103 + 26390 * k), CFG_MPFR_RND);
	 *     mpfr_mul(__impl->term_dividend, __impl->curr_fact_4k, __impl->t0, CFG_MPFR_RND);
	 *
	 * (2) mpfr_mul_ui(__impl->term_dividend, __impl->curr_fact_4k, (1103 + 26390 * k), CFG_MPFR_RND);
	 */
#if 0
	mpfr_set_ui(__impl->t0, 26390UL, CFG_MPFR_RND);
	mpfr_mul_ui(__impl->t0, __impl->t0, k, CFG_MPFR_RND);
	mpfr_add_ui(__impl->t0, __impl->t0, 1103UL, CFG_MPFR_RND);
	/* t0 has (1103 + 26390 * k) */
	mpfr_mul(__impl->term_dividend, __impl->curr_fact_4k, __impl->t0, CFG_MPFR_RND);
#endif
	mpfr_mul_ui(__impl->term_dividend, __impl->curr_fact_4k, (1103UL + 26390UL * k), CFG_MPFR_RND);
	/* term_dividend has 4k! * (1103 + 26390 * k) */

	/* term_dividend calculated */
	//printf("make_pi: term_dividend(%d) = ", k);
	//mpfr_out_str(stdout, 10, 0, term_dividend, CFG_MPFR_RND);
	//printf("\n");

	/*
	 * calculate divisor
	 *
	 * [ (FACT(k) ^ 4) * (396 ^ (4 * k)) ]
	 *
	 */
	mpfr_pow_ui(__impl->term_divisor, __impl->curr_fact_k, 4UL, CFG_MPFR_RND);
	/* term_divisor now has ((k!) ^ 4) */
	mpfr_set_ui(__impl->t0, 396UL, CFG_MPFR_RND);
	mpfr_pow_ui(__impl->t0, __impl->t0, _4k, CFG_MPFR_RND);
	/* t0 has (396 ^ (4 * k)) */
	mpfr_mul(__impl->term_divisor, __impl->term_divisor, __impl->t0, CFG_MPFR_RND);
	/* term_divisor calculated */

	//printf("make_pi: term_divisor(%d) = ", k);
	//mpfr_out_str(stdout, 10, 0, term_divisor, CFG_MPFR_RND);
	//printf("\n");

	/*
	 * calculate term
	 */
	mpfr_div(__impl->term, __impl->term_dividend, __impl->term_divisor, CFG_MPFR_RND);
	//printf("make_pi: term(%d) = ", k);
	//mpfr_out_str(stdout, 10, 0, term, CFG_MPFR_RND);
	//printf("\n");

	/* term_dividend holds term */

	/*
	 * calculate term_sum
	 */
	mpfr_add(__impl->term_sum, __impl->term_sum, __impl->term, CFG_MPFR_RND);
	//printf("make_pi: term_sum(%d) = ", k);
	//mpfr_out_str(stdout, 10, 0, term_sum, CFG_MPFR_RND);
	//printf("\n");


	/*
	 * calculate out values and retval.
	 */
	*out_k = __impl->curr_k;


	/*
	 * setup for next iteration.
	 */
	/*
	 * next K
	 */
	__impl->curr_k += 1;
	/*
	 * next 4*K
	 */
	__impl->curr_4k += 4;
	/*
	 * compute FACT(k) as:
	 * 		FACT(k) = FACT(k - 1) * k
	 */
	mpfr_mul_ui(__impl->curr_fact_k, __impl->curr_fact_k, __impl->curr_k, CFG_MPFR_RND);
	/*
	 * compute FACT4(k) as:
	 *              FACT4(4k) = FACT(PREV(4k)) * 4k * (4k - 1) * (4k - 2) * (4k - 3)
	 */
	mpfr_mul_ui(__impl->curr_fact_4k, __impl->curr_fact_4k, __impl->curr_4k, CFG_MPFR_RND);
	mpfr_mul_ui(__impl->curr_fact_4k, __impl->curr_fact_4k, (__impl->curr_4k - 1UL), CFG_MPFR_RND);
	mpfr_mul_ui(__impl->curr_fact_4k, __impl->curr_fact_4k, (__impl->curr_4k - 2UL), CFG_MPFR_RND);
	mpfr_mul_ui(__impl->curr_fact_4k, __impl->curr_fact_4k, (__impl->curr_4k - 3UL), CFG_MPFR_RND);

	return &__impl->term_sum;
}

static const void pi_impl_ramanujan_1910_opt_get_value(const mpfr_t *partial_sum, mpfr_t *result_out)
{
	int ret;

	/*
	 * CMULT = (2 * sqrt(2)) / 9801			      # constant
	 *                                                    # 9801 = 99^2
	 */
	/*
	 * calculate CMULT constant directly into result
	 */
	mpfr_sqrt_ui((*result_out), 2UL, CFG_MPFR_RND);
	mpfr_mul_ui((*result_out), (*result_out), 2UL, CFG_MPFR_RND);
	mpfr_div_ui((*result_out), (*result_out), 9801UL, CFG_MPFR_RND);
	// printf("make_pi: cmult = ");
	// mpfr_out_str(stdout, 10, 0, cmult, CFG_MPFR_RND);
	// printf("\n");

	/*
	 * make sure to use max_k in digits estimation,
	 * as k is above k just to add some slack space.
	 */

	/*
	 * calculate 1 / PI = cmult * term_sum
	 */
	mpfr_mul((*result_out), (*result_out), (*partial_sum), CFG_MPFR_RND);
	//printf("1/pi(%d) = ", k);
	//mpfr_out_str(stdout, 10, 0, pi, CFG_MPFR_RND);
	//printf("\n");
	/*
	 * calculate PI from 1 / PI
	 */
	mpfr_ui_div((*result_out), 1UL, (*result_out), CFG_MPFR_RND);
	//printf("pi(%d - %d) = %s\n", k, k * 8, get_pi_value(pi, k*8));
	//mpfr_out_str(stdout, 10, 0, pi, CFG_MPFR_RND);
	//printf("\n");
}
