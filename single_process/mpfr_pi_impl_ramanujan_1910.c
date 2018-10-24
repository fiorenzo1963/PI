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
 * Srinivasa Ramanujan 1910 formula:
 *
 * More info on Ramanujan PI formulas:
 * https://en.wikipedia.org/wiki/Srinivasa_Ramanujan
 * https://en.wikipedia.org/wiki/Approximations_of_%CF%80
 * https://en.wikipedia.org/wiki/Ramanujan%E2%80%93Sato_series
 *
 * 1/PI = CMULT * SUM(k, 0..infinity) TERM(k)
 * 
 * CMULT = (2 * sqrt(2)) / 9801			      # constant
 *                                                    # 9801 = 99^2
 *
 * TERM(k) = [ (4 * k)! * (1103 + 26390 * k) ] /      # dividend
 *           [ ((k!) ^ 4) * (396 ^ (4 * k)) ]         # divisor
 *                                                    # 396 = 99 * 4
 */

static const char *pi_impl_ramanujan_1910_get_name(void)
{
	return "Ramanujan 1910 Formula";
}

static void pi_impl_ramanujan_1910_deinitialize(struct mpfr_pi_impl *impl);
static int pi_impl_ramanujan_1910_compute_next_term(struct mpfr_pi_impl *impl, unsigned long *out_k, long *digits_out);
static mpfr_t *pi_impl_ramanujan_1910_get_value(struct mpfr_pi_impl *impl, long *digits_out);

/* actual implementation struct for this algorithm */
struct __mpfr_pi_impl {
	/* generic part */
	struct mpfr_pi_impl g;
	/* private part */
	unsigned long curr_k; /* current iteration -- compute Ki must be called with this iteration number. compute Ki will increment k */
	long curr_digits;
	long desired_digits; /* desired digits */
	unsigned long max_k; /* max_k to reach desired digits */
	/* various state variables needed */
	/*
	 * temp variables reused at each iteration
	 */
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

#define DIGITS_TO_K(d)	(((d) / 8) + 1)		/* number of iterations to get "d" digits */
#define SLACK_K		DIGITS_TO_K(8)		/* slack factor added to the above just to be sure */
/* arg reused, must pass l-value */
#define K_TO_DIGITS(k)	((k) >= SLACK_K ? ((k) * 8) - SLACK_K : 0L)

#define __SAFE_LONG_MAX		(LONG_MAX / 2L)
#define __SAFE_ULONG_MAX	(ULONG_MAX / 2UL)

struct mpfr_pi_impl *pi_impl_ramanujan_1910_initialize(const long digits, unsigned long *out_max_k)
{
	struct __mpfr_pi_impl *__impl = malloc(sizeof (struct __mpfr_pi_impl));
	assert(__impl != NULL);
	__impl->g.f_impl_get_name = pi_impl_ramanujan_1910_get_name;
	__impl->g.f_initialize = pi_impl_ramanujan_1910_initialize;
	__impl->g.f_deinitialize = pi_impl_ramanujan_1910_deinitialize;
	__impl->g.f_pi_compute_next_term = pi_impl_ramanujan_1910_compute_next_term;
	__impl->g.f_pi_get_value = pi_impl_ramanujan_1910_get_value;

	__impl->curr_k = 0UL;
	__impl->curr_digits = 0L;
	__impl->desired_digits = digits;
	printf("pi_impl_ramanujan_1910_initialize: desired digits = %ld\n", __impl->desired_digits);
	/* iterations needed */
	__impl->max_k = DIGITS_TO_K(digits) + SLACK_K;
	assert(digits < __SAFE_LONG_MAX);
	assert(__impl->max_k < __SAFE_ULONG_MAX);
	/* algorithm computes 4k directly with unsigned longs */
	assert(__impl->max_k < __SAFE_ULONG_MAX / 4);
	printf("pi_impl_ramanujan_1910_initialize: max_k = %lu\n", __impl->max_k);
	/* various state variables needed */
	mpfr_init2(__impl->term_dividend, CFG_MPFR_PREC);
	mpfr_init2(__impl->term_divisor, CFG_MPFR_PREC);
	mpfr_init2(__impl->term, CFG_MPFR_PREC);
	mpfr_init2(__impl->term_sum, CFG_MPFR_PREC);
	mpfr_init2(__impl->cmult, CFG_MPFR_PREC);
	mpfr_init2(__impl->pi, CFG_MPFR_PREC);
	mpfr_init2(__impl->t0, CFG_MPFR_PREC);

	/*
	 * CMULT = (2 * sqrt(2)) / 9801			      # constant
	 *                                                    # 9801 = 99^2
	 */
	/* calculate cmult constant */
	mpfr_sqrt_ui(__impl->cmult, 2UL, CFG_MPFR_RND);
	mpfr_mul_ui(__impl->cmult, __impl->cmult, 2UL, CFG_MPFR_RND);
	mpfr_div_ui(__impl->cmult, __impl->cmult, 9801UL, CFG_MPFR_RND);
	// printf("make_pi: cmult = ");
	// mpfr_out_str(stdout, 10, 0, cmult, CFG_MPFR_RND);
	// printf("\n");

	/* set term_sum */
	mpfr_set_ui(__impl->term_sum, 0, CFG_MPFR_RND);

	*out_max_k = __impl->max_k;
	return (struct mpfr_pi_impl *)__impl;
}

static void pi_impl_ramanujan_1910_deinitialize(struct mpfr_pi_impl *impl)
{
	struct __mpfr_pi_impl *__impl = (struct __mpfr_pi_impl *)impl;
	
	mpfr_clear(__impl->term_dividend);
	mpfr_clear(__impl->term_divisor);
	mpfr_clear(__impl->term);
	mpfr_clear(__impl->term_sum);
	mpfr_clear(__impl->cmult);
	mpfr_clear(__impl->pi);
	mpfr_clear(__impl->t0);
	free(__impl);
}

static int pi_impl_ramanujan_1910_compute_next_term(struct mpfr_pi_impl *impl, unsigned long *out_k, long *digits_out)
{
	struct __mpfr_pi_impl *__impl = (struct __mpfr_pi_impl *)impl;
	const unsigned long k = __impl->curr_k;
	int ret;

	/*
	 * calculate dividend
	 *
	 * [ (4 * k)! * (1103 + 26390 * k) ]
	 *
	 */

	mpfr_fac_ui(__impl->term_dividend, (4UL * k), CFG_MPFR_RND);
	/* term_dividend now has (4*k)! */
	mpfr_set_ui(__impl->t0, 26390UL, CFG_MPFR_RND);
	mpfr_mul_ui(__impl->t0, __impl->t0, k, CFG_MPFR_RND);
	mpfr_add_ui(__impl->t0, __impl->t0, 1103UL, CFG_MPFR_RND);
	/* t0 has (1103 + 26390 * k) */
	mpfr_mul(__impl->term_dividend, __impl->term_dividend, __impl->t0, CFG_MPFR_RND);
	/* term_dividend calculated */
	//printf("make_pi: term_dividend(%d) = ", k);
	//mpfr_out_str(stdout, 10, 0, term_dividend, CFG_MPFR_RND);
	//printf("\n");

	/*
	 * calculate divisor
	 *
	 * [ ((k!) ^ 4) * (396 ^ (4 * k)) ]
	 *
	 */
	mpfr_fac_ui(__impl->term_divisor, k, CFG_MPFR_RND);
	mpfr_pow_ui(__impl->term_divisor, __impl->term_divisor, 4UL, CFG_MPFR_RND);
	/* term_divisor now has ((k!) ^ 4) */
	mpfr_set_ui(__impl->t0, 396UL, CFG_MPFR_RND);
	mpfr_pow_ui(__impl->t0, __impl->t0, (4UL * k), CFG_MPFR_RND);
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
	*digits_out = K_TO_DIGITS(__impl->curr_k);
	ret = (__impl->curr_k >= __impl->max_k) ? 1 : 0;
	/*
	 * setup for next iteration.
	 */
	__impl->curr_k++;

	return ret;
}

static mpfr_t *pi_impl_ramanujan_1910_get_value(struct mpfr_pi_impl *impl, long *digits_out)
{
	struct __mpfr_pi_impl *__impl = (struct __mpfr_pi_impl *)impl;
	int ret;

	/*
	 * use (curr_k - 1) as when we call this, curr_k has not been done yet
	 */

	// printf("pi_impl_ramanujan_1910_get_value: ret=%d, curr_k-1=%ld, max_k=%ld\n", ret, (__impl->curr_k - 1),__impl->max_k);
	if (__impl->curr_k == 0UL || K_TO_DIGITS(__impl->curr_k - 1) == 0) {
		*digits_out = 0L;
		return NULL;
	}

	/*
	 * make sure to use max_k in digits estimation,
	 * as k is above k just to add some slack space.
	 */

	/*
	 * calculate 1 / PI = cmult * term_sum
	 */
	mpfr_mul(__impl->pi, __impl->cmult, __impl->term_sum, CFG_MPFR_RND);
	//printf("1/pi(%d) = ", k);
	//mpfr_out_str(stdout, 10, 0, pi, CFG_MPFR_RND);
	//printf("\n");
	/*
	 * calculate PI from 1 / PI
	 */
	mpfr_ui_div(__impl->pi, 1UL, __impl->pi, CFG_MPFR_RND);
	//printf("pi(%d - %d) = %s\n", k, k * 8, get_pi_value(pi, k*8));
	//mpfr_out_str(stdout, 10, 0, pi, CFG_MPFR_RND);
	//printf("\n");

	*digits_out = K_TO_DIGITS(__impl->curr_k - 1);

	return &__impl->pi;
}
