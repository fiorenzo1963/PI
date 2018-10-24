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
	 * sets iteration K value to 0 and other implementation specific constants
	 * set out_max_k to rhe number of iterations (i.e., max K value), if known, otherwise 0.
	 * if the implementation optimizes this calculation and keeps intermediate state, initialize this state.
	 */
	struct mpfr_pi_impl * (*f_initialize)(const long digits, unsigned long *out_max_k);
	/*
	 * free an implementation struct.
	 */
	void (*f_deinitialize)(struct mpfr_pi_impl *impl);
	/*
	 * compute K(i) term of the series.
	 * the implementation is free to optimize this calculation and keep intermediate state.
	 * return 0 if more iterations are needed, 1 if the desired precision has been reached.
	 * sets computed k in k_out.
	 * set digits_out if the digits are known, otherwise sets to 0.
	 * (may set digits_out every now and then, so it's not guaranteed to be updated at every iteration).
	 */
	int (*f_pi_compute_next_term)(struct mpfr_pi_impl *impl, unsigned long *k_out, long *digits_out);
	/*
	 * this can be called anytime, return NULL if estimated digits is still 0.
	 * it computes PI based on the current values.
	 * if called after f_pi_computer_term returns 1, it's guaranteed to have at least
	 * the desired number of digits of precision.
	 */
	 mpfr_t * (*f_pi_get_value)(struct mpfr_pi_impl *impl, long *digits_out);
};

#endif
