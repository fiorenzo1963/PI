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
#define CFG_MPFR_PREC		1000000
#define CFG_MPFR_RND		MPFR_RNDD
/* derived from above */
#define CFG_MPFR_PREC_STR	__stringify(CFG_MPFR_PREC)

static char *mpfr_t_to_str(mpfr_t *value, int digits)
{
	char *buf;
	assert(value != NULL);
	assert(digits > 0);
	// printf("get_float_to_str(prec=%d, digits=%d)\n", CFG_MPFR_PREC, digits);
	buf = malloc(digits + 10);
	assert(buf != NULL);
	mpfr_snprintf(buf, digits, "%." CFG_MPFR_PREC_STR "R*f", CFG_MPFR_RND, *value);
	return buf;
}

static void free_mpfr_str(char *buf)
{
	free(buf);
}

struct mpfr_pi_impl {
	/*
	 * initialize an implementation and return its struct
	 * the init function will add extra state variables to this struct, so do not make any
	 * size assumptions on it.
	 */
	struct mpfr_pi_impl * (*f_initialize)(void);
	/* free an implementation struct */
	void (*f_deinitialize)(struct mpfr_pi_impl *impl);
	/*
	 * start PI calculation. sets K value to 0 and other implementation specific constants
	 * returns number of estimated iterations, if known, otherwise returns 0.
	 * if the implementation optimizes this calculation and keeps intermediate state, initialize this state.
	 */
	void (*f_start_pi_compute)(struct mpfr_pi_impl *impl, const int digits, int *estimated_iterations);
	/*
	 * compute K(i) term of the series.
	 * the implementation is free to optimize this calculation and keep intermediate state.
	 * return 0 if more iterations are needed, 1 if the desired precision has been reached.
	 */
	int (*f_pi_compute_term)(struct mpfr_pi_impl *impl, const int Ki);
	/*
	 * this can be called anytime. it computes PI based on the current values.
	 * if called after f_pi_computer_term returns 1, it's guaranteed to have at least
	 * the desired number of digits of precision.
	 */
	int (*f_pi_get_current_pi)(struct mpfr_pi_impl *impl, const int Ki);
};

#endif
