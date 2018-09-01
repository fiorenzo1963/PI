#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <assert.h>
#include <mpfr.h>

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
 * TERM(k) = [ (4 * k)! * (1103 + 26390 * k) ] /      # divend
 *           [ ((k!) ^ 4) * (396 ^ (4 * k)) ]         # divisor
 *                                                    # 396 = 99 * 4
 */


#if 0
struct mpfr_pi_impl {
	/*
	 * return name of the implementation.
	 */
	const char * (*f_impl_get_name)(void);
	/*
	 * initialize an implementation and return its struct
	 * the init function will add extra state variables to this struct, so do not make any
	 * size assumptions on it.
	 * sets iteration K value to 0 and other implementation specific constants
	 * set out_iterations to rhe number of estimated iterations, if known, otherwise 0.
	 * if the implementation optimizes this calculation and keeps intermediate state, initialize this state.
	 */
	struct mpfr_pi_impl * (*f_initialize)(const long digits, long *out_iterations);
	/*
	 * free an implementation struct.
	 */
	void (*f_deinitialize)(struct mpfr_pi_impl *impl);
	/*
	 * compute K(i) term of the series.
	 * the implementation is free to optimize this calculation and keep intermediate state.
	 * return 0 if more iterations are needed, 1 if the desired precision has been reached.
	 * sets computed k in k_out.
	 * set digits_out if the digits are known, otherwise sets to 0
	 * (may set digits_out every now and then, so it's not guaranteed to be updated at every iteration).
	 */
	int (*f_pi_compute_next_term)(struct mpfr_pi_impl *impl, long *ki_out, long *digits_out);
	/*
	 * this can be called anytime, return NULL if k is still 0.
	 * it computes PI based on the current values.
	 * if called after f_pi_computer_term returns 1, it's guaranteed to have at least
	 * the desired number of digits of precision.
	 * caller must not modify the returned value in any way.
	 */
	 mpfr_t * (*f_pi_get_value)(struct mpfr_pi_impl *impl, long *digits_out);
};
#endif

const char *pi_impl_ramananujan_1910_get_name(void)
{
	return "Ramananujan 1910 Formula";
}

void pi_impl_ramananujan_1910_deinitialize(struct mpfr_pi_impl *impl);
int pi_impl_ramananujan_1910_compute_next_term(struct mpfr_pi_impl *impl, long *ki_out, long *digits_out);
mpfr_t *pi_impl_ramananujan_1910_get_value(struct mpfr_pi_impl *impl, long *digits_out);

/* actual implementation struct for this algorithm */
struct __mpfr_pi_impl {
	/* generic part */
	struct mpfr_pi_impl g;
	/* private part */
	long curr_k; /* current iteration -- compute Ki must be called with this iteration number. compute Ki will increment k */
	long curr_digits;
	long desired_digits; /* desired digits */
	long max_k; /* max_k to reach desired digits */
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

#define DIGITS_TO_K(d)	(((d) / 8) + 1)		/* estimate number of iterations to get "d" digits */
#define SLACK_K		DIGITS_TO_K(8)		/* slack factor added to the above just to be sure */
/* arg reused, must pass l-value */
#define K_TO_DIGITS(k)	((k) >= SLACK_K ? ((k) * 8) - SLACK_K : 0L)

struct mpfr_pi_impl *pi_impl_ramananujan_1910_initialize(const long digits, long *out_iterations)
{
	struct __mpfr_pi_impl *__impl = malloc(sizeof (struct __mpfr_pi_impl));
	assert(__impl != NULL);
	__impl->g.f_impl_get_name = pi_impl_ramananujan_1910_get_name;
	__impl->g.f_initialize = pi_impl_ramananujan_1910_initialize;
	__impl->g.f_deinitialize = pi_impl_ramananujan_1910_deinitialize;
	__impl->g.f_pi_compute_next_term = pi_impl_ramananujan_1910_compute_next_term;
	__impl->g.f_pi_get_value = pi_impl_ramananujan_1910_get_value;

	__impl->curr_k = 0L;
	__impl->curr_digits = 0L;
	__impl->desired_digits = digits;
	printf("pi_impl_ramananujan_1910_initialize: desired digits = %ld\n", __impl->desired_digits);
	/* iterations needed */
	__impl->max_k = DIGITS_TO_K(digits) + SLACK_K;
	printf("pi_impl_ramananujan_1910_initialize: estimated iterations = %ld\n", __impl->max_k);
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
	mpfr_sqrt_ui(__impl->cmult, 2, CFG_MPFR_RND);
	mpfr_mul_ui(__impl->cmult, __impl->cmult, 2, CFG_MPFR_RND);
	mpfr_div_ui(__impl->cmult, __impl->cmult, 9801, CFG_MPFR_RND);
	// printf("make_pi: cmult = ");
	// mpfr_out_str(stdout, 10, 0, cmult, CFG_MPFR_RND);
	// printf("\n");

	/* set term_sum */
	mpfr_set_ui(__impl->term_sum, 0, CFG_MPFR_RND);

	*out_iterations = __impl->max_k;
	return (struct mpfr_pi_impl *)__impl;
}

void pi_impl_ramananujan_1910_deinitialize(struct mpfr_pi_impl *impl)
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

int pi_impl_ramananujan_1910_compute_next_term(struct mpfr_pi_impl *impl, long *ki_out, long *digits_out)
{
	struct __mpfr_pi_impl *__impl = (struct __mpfr_pi_impl *)impl;
	long k = __impl->curr_k;
	int ret;

	/*
	 * calculate dividend
	 *
	 * [ (4 * k)! * (1103 + 26390 * k) ]
	 *
	 */

	mpfr_fac_ui(__impl->term_dividend, (4 * k), CFG_MPFR_RND);
	/* term_dividend now has (4*k)! */
	mpfr_set_ui(__impl->t0, 26390, CFG_MPFR_RND);
	mpfr_mul_ui(__impl->t0, __impl->t0, k, CFG_MPFR_RND);
	mpfr_add_ui(__impl->t0, __impl->t0, 1103, CFG_MPFR_RND);
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
	mpfr_pow_ui(__impl->term_divisor, __impl->term_divisor, 4, CFG_MPFR_RND);
	/* term_divisor now has ((k!) ^ 4) */
	mpfr_set_ui(__impl->t0, 396, CFG_MPFR_RND);
	mpfr_pow_ui(__impl->t0, __impl->t0, (4 * k), CFG_MPFR_RND);
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
	*ki_out = __impl->curr_k;
	*digits_out = K_TO_DIGITS(__impl->curr_k);
	ret = (__impl->curr_k >= __impl->max_k) ? 1 : 0;
	/*
	 * setup for next iteration.
	 */
	__impl->curr_k++;

	return ret;
}

mpfr_t *pi_impl_ramananujan_1910_get_value(struct mpfr_pi_impl *impl, long *digits_out)
{
	struct __mpfr_pi_impl *__impl = (struct __mpfr_pi_impl *)impl;
	int ret;

	/* curr_k is the next iteration, so use (curr_k - 1) */
	ret = ((__impl->curr_k - 1) >= __impl->max_k) ? 1 : 0;
	// printf("pi_impl_ramananujan_1910_get_value: ret=%d, curr_k-1=%ld, max_k=%ld\n", ret, (__impl->curr_k - 1),__impl->max_k);
	if (ret == 0) {
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
	mpfr_ui_div(__impl->pi, 1, __impl->pi, CFG_MPFR_RND);
	//printf("pi(%d - %d) = %s\n", k, k * 8, get_pi_value(pi, k*8));
	//mpfr_out_str(stdout, 10, 0, pi, CFG_MPFR_RND);
	//printf("\n");

	*digits_out = K_TO_DIGITS(__impl->curr_k - 1);

	return &__impl->pi;
}
