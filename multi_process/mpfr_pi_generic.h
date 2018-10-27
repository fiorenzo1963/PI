#ifndef _MPFR_PI_GENERIC_
#define _MPRF_PI_GENERIC_


#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>
#include <assert.h>
#include <mpfr.h>

#include "stringify.h"

/* config variables, can be changed */
#define CFG_MPFR_PREC		10000000
#define CFG_MPFR_RND		MPFR_RNDD
/* derived from above */
#define CFG_MPFR_PREC_STR	__stringify(CFG_MPFR_PREC)

#define MAX_DIGITS_PREC		((unsigned long)(((double)CFG_MPFR_PREC / 4) * 0.95))

static char *mpfr_t_to_str(mpfr_t *value, long chars)
{
	char *buf;
	assert(value != NULL);
	assert(chars > 0);
	// printf("get_float_to_str(prec=%d, digits=%d)\n", CFG_MPFR_PREC, chars);
	buf = malloc(chars + 100);
	assert(buf != NULL);
	mpfr_snprintf(buf, chars, "%." CFG_MPFR_PREC_STR "R*f", CFG_MPFR_RND, *value);
	return buf;
}

static void free_mpfr_str(char *buf)
{
	free(buf);
}

struct mpfr_pi_impl {
	/*
	 * return name of the implementation.
	 */
	const char * (*f_impl_get_name)(void);
	/*
	 * initialize an implementation and return its struct
	 * the init function will add extra state variables to this struct, so do not make any
	 * size assumptions on it.
	 * sets initial iteration K value to "start_k" and other implementation specific constants.
	 * if the implementation optimizes this calculation and keeps intermediate state, initialize this state.
	 * digits is the approximate maximum number of digits the calculation wants to achieve in total.
	 */
	struct mpfr_pi_impl * (*f_initialize)(const unsigned long digits, const unsigned long k_start);
	/*
	 * free an implementation struct.
	 */
	void (*f_deinitialize)(struct mpfr_pi_impl *impl);
	/*
	 * compute K(i) term of the series.
	 * the implementation is free to optimize this calculation and keep intermediate state.
	 * sets computed k in k_out, which means k_out gets incremented by 1 on each call.
	 * if this is called N times, the k values partial sum will be from [k_start .. k_start + N].
	 * (may set digits_out every now and then, so it's not guaranteed to be updated at every iteration).
	 *
	 * returns partial sum of [k_start .. k_start + N].
	 */
	const mpfr_t * (*f_pi_compute_next_term)(struct mpfr_pi_impl *impl, unsigned long *k_out);
	/*
	 * this can be called anytime, 
	 * caller must add all the partial sums to have to have a full sequence from 0 to M into
	 * partial_sum, and will get an approximate value for PI in return into result_out.
	 * Ramanujan's algorithm only yields exactly 8 digits per iteration up to 100K digits or so,
	 * at which point is begins be less than 8 (something like 98% off 8 dights per iteration after it).
	 */
	 void (*f_pi_get_value)(const mpfr_t *partial_sum, mpfr_t *result_out);
};

#endif
