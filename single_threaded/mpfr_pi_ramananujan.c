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
 * 1/PI = CMULT * SUM(k, 0..infinity) TERM(k)
 * 
 * CMULT = (2 * sqrt(2)) / 9801			      # constant
 *                                                    # 9801 = 99^2
 *
 * TERM(k) = [ (4 * k)! * (1103 + 26390 * k) ] /      # divend
 *           [ ((k!) ^ 4) * (396 ^ (4 * k)) ]         # divisor
 *                                                    # 396 = 99 * 4
 */
#define DIGITS_TO_K(d)	((d) / 8)
#define SLACK_K		DIGITS_TO_K(128)

#define CHARACTERS_PER_LINE	100

#define MPFR_PREC_USED		1000000
#define MPFR_PREC_USED_STR	__stringify(MPFR_PREC_USED)

char *get_float_to_str(mpfr_t *pi, int digits)
{
	char *buf;
	assert(pi != NULL);
	assert(digits > 0);
	// printf("get_float_to_str(prec=%d, digits=%d)\n", MPFR_PREC_USED, digits);
	buf = malloc(digits + 10);
	assert(buf != NULL);
	mpfr_snprintf(buf, digits, "%." MPFR_PREC_USED_STR "R*f", MPFR_RNDD, *pi);
	return buf;
}

void free_float_str(char *buf)
{
	free(buf);
}

void print_pi(const char *pi_string)
{
	size_t i, sz = strlen(pi_string);
	for (i = 0; i < sz; ) {
		char buf[CHARACTERS_PER_LINE + 1];
		int cc = CHARACTERS_PER_LINE;
		if ((sz - i) < cc)
			cc = sz - i;
		strncpy(buf, &pi_string[i], cc);
		buf[cc] = '\0';
		printf("%s\n", buf);
		i += cc;
	}
	assert(i == sz);
}

void make_pi(int digits)
{
	int k, max_k;
	int kt0, tk1;
	uint64_t time0, time1, time2;
        uint64_t tss3, tss4;
	char datebuf[128];
	char offsetbuf[128];

	mpfr_t term_dividend;
	mpfr_t term_divisor;
	mpfr_t term;
	mpfr_t term_sum;
	mpfr_t cmult;
	mpfr_t pi;
	mpfr_t t0;

	char *s;

	mpfr_init2(term_dividend, MPFR_PREC_USED);
	mpfr_init2(term_divisor, MPFR_PREC_USED);
	mpfr_init2(term, MPFR_PREC_USED);
	mpfr_init2(term_sum, MPFR_PREC_USED);
	mpfr_init2(cmult, MPFR_PREC_USED);
	mpfr_init2(pi, MPFR_PREC_USED);
	mpfr_init2(t0, MPFR_PREC_USED);

	/* calculate cmult constant */
	mpfr_sqrt_ui(cmult, 2, MPFR_RNDD);
	mpfr_mul_ui(cmult, cmult, 2, MPFR_RNDD);
	mpfr_div_ui(cmult, cmult, 9801, MPFR_RNDD);
	// printf("make_pi: cmult = ");
	// mpfr_out_str(stdout, 10, 0, cmult, MPFR_RNDD);
	// printf("\n");

	/* set term_sum */
	mpfr_set_ui(term_sum, 0, MPFR_RNDD);

	time0 = gettimestamp_nsecs();
	ts_to_date_str(datebuf, sizeof (datebuf), time0);
	ts_to_offset_str(offsetbuf, sizeof (offsetbuf), time0 - time0);
	max_k = DIGITS_TO_K(digits);

	printf("%s: %s: make_pi, digits = %d, max k = %d\n", datebuf, offsetbuf, digits, max_k);

	/*
	 * the extra iterations are not technically necessary, but just to be safe .....
	 */
	tss3 = gettimestamp_nsecs();
	for (k = 0; k <= max_k + SLACK_K; k++) {

		// printf("make_pi: k = %d (estimated digits = %d)\n", k, (k+1)*8);


		/*
		 * calculate dividend
		 *
		 * [ (4 * k)! * (1103 + 26390 * k) ]
		 *
		 */

		mpfr_fac_ui(term_dividend, (4 * k), MPFR_RNDD);
		/* term_dividend now has (4*k)! */
		mpfr_set_ui(t0, 26390, MPFR_RNDD);
		mpfr_mul_ui(t0, t0, k, MPFR_RNDD);
		mpfr_add_ui(t0, t0, 1103, MPFR_RNDD);
		/* t0 has (1103 + 26390 * k) */
		mpfr_mul(term_dividend, term_dividend, t0, MPFR_RNDD);
		/* term_dividend calculated */
		//printf("make_pi: term_dividend(%d) = ", k);
		//mpfr_out_str(stdout, 10, 0, term_dividend, MPFR_RNDD);
		//printf("\n");

		/*
		 * calculate divisor
		 *
		 * [ ((k!) ^ 4) * (396 ^ (4 * k)) ]
		 *
		 */
		mpfr_fac_ui(term_divisor, k, MPFR_RNDD);
		mpfr_pow_ui(term_divisor, term_divisor, 4, MPFR_RNDD);
		/* term_divisor now has ((k!) ^ 4) */
		mpfr_set_ui(t0, 396, MPFR_RNDD);
		mpfr_pow_ui(t0, t0, (4 * k), MPFR_RNDD);
		/* t0 has (396 ^ (4 * k)) */
		mpfr_mul(term_divisor, term_divisor, t0, MPFR_RNDD);
		/* term_divisor calculated */
		//printf("make_pi: term_divisor(%d) = ", k);
		//mpfr_out_str(stdout, 10, 0, term_divisor, MPFR_RNDD);
		//printf("\n");

		/*
		 * calculate term
		 */
		mpfr_div(term, term_dividend, term_divisor, MPFR_RNDD);
		//printf("make_pi: term(%d) = ", k);
		//mpfr_out_str(stdout, 10, 0, term, MPFR_RNDD);
		//printf("\n");

		/*
		 * calculate term_sum
		 */
		mpfr_add(term_sum, term_sum, term, MPFR_RNDD);
		//printf("make_pi: term_sum(%d) = ", k);
		//mpfr_out_str(stdout, 10, 0, term_sum, MPFR_RNDD);
		//printf("\n");

		tss4 = gettimestamp_nsecs();
		if (ts_secs_portion(tss4 - tss3) >= 10) {
			ts_to_date_str(datebuf, sizeof (datebuf), gettimestamp_nsecs());
			ts_to_offset_str(offsetbuf, sizeof (offsetbuf), tss4 - tss3);
			printf("%s: %s: k = %d, max_k = %d\n", datebuf, offsetbuf, k, max_k);
			tss3 = tss4;
		}
	}

	tss4 = gettimestamp_nsecs();
	ts_to_date_str(datebuf, sizeof (datebuf), gettimestamp_nsecs());
	ts_to_offset_str(offsetbuf, sizeof (offsetbuf), tss4 - tss3);
	printf("%s: %s: k = %d, max_k = %d\n", datebuf, offsetbuf, k, max_k);
	tss3 = tss4;

	/*
	 * make sure to use max_k in digits estimation,
	 * as k is above k just to add some slack space.
	 */

	/*
	 * calculate 1 / PI = cmult * term_sum
	 */
	mpfr_mul(pi, cmult, term_sum, MPFR_RNDD);
	//printf("1/pi(%d) = ", k);
	//mpfr_out_str(stdout, 10, 0, pi, MPFR_RNDD);
	//printf("\n");
	/*
	 * calculate PI from 1 / PI
	 */
	mpfr_ui_div(pi, 1, pi, MPFR_RNDD);
	//printf("pi(%d - %d) = %s\n", k, k * 8, get_pi_value(pi, k*8));
	//mpfr_out_str(stdout, 10, 0, pi, MPFR_RNDD);
	//printf("\n");

	time1 = gettimestamp_nsecs();

	/*
	 * print PI.
	 * conversion from internal binary representation to decimal takes a long time.
	 */
	s = get_float_to_str(&pi, digits + 1); /* the +1 accounts for the decimal point */

	time2 = gettimestamp_nsecs();
	ts_to_date_str(datebuf, sizeof (datebuf), time2);
	ts_to_offset_str(offsetbuf, sizeof (offsetbuf), time2 - tss3);
	printf("%s: %s: (finalization and conversion base 10)\n", datebuf, offsetbuf);
	printf("pi(k = %d, d = %d):\n", max_k, digits);
	printf("\n");
	print_pi(s);

	free_float_str(s);

	mpfr_clear(term_dividend);
	mpfr_clear(term_divisor);
	mpfr_clear(term);
	mpfr_clear(term_sum);
	mpfr_clear(cmult);
	mpfr_clear(pi);
	mpfr_clear(t0);
	mpfr_free_cache();
}

int main(int argc, char **argv)
{
	int digits;

	printf("MPFR library: %-12s\nMPFR header:  %s (based on %d.%d.%d)\n",
	       mpfr_get_version (), MPFR_VERSION_STRING, MPFR_VERSION_MAJOR,
	       MPFR_VERSION_MINOR, MPFR_VERSION_PATCHLEVEL);
	printf("MPFR_PREC_MAX = %ld\n", MPFR_PREC_MAX);
	printf("\n");
	printf("MPFR_PREC_USED = %d\n", MPFR_PREC_USED);
	printf("mpfr_custom_get_size(MPFR_PREC_USED) = %ld\n", mpfr_custom_get_size(MPFR_PREC_USED));
	printf("approximate decimals for MPFR_PREC_USED = %ld (upper bound)\n", (long)((float)MPFR_PREC_USED / 3.5));
	printf("\n");

	if (argc != 2)
		digits = 1000;
	else
		digits = (int)strtoul(argv[1], NULL, 0);
	if (digits < 8) {
		printf("invalid %d parameter for digits (minimum value is 8)\n", digits);
		exit(1);
	}
	if ((long)digits >= (long)((float)MPFR_PREC_USED / 3.5)) {
		printf("this build does not support %d digits\n", digits);
		exit(1);
	}

	/* round up to the next multiple of 8 */
	digits += 07;
	digits &= ~07;

	printf("calculating pi to %d digits\n", digits);

	make_pi(digits);

	return 0;
}
