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
 * using various PI expansion formulas. See each implementation for details.
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

#define CHARACTERS_PER_LINE	100

#define RAMANUJAN_1910_OPT_ALGO	0
extern struct mpfr_pi_impl *pi_impl_ramanujan_1910_opt_initialize(unsigned long digits, const unsigned long k_start);

void writeout_pi(FILE *fd, const char *pi_string)
{
	size_t i, sz = strlen(pi_string);
	for (i = 0; i < sz; ) {
		char buf[CHARACTERS_PER_LINE + 1];
		int cc = CHARACTERS_PER_LINE;
		if ((sz - i) < cc)
			cc = sz - i;
		strncpy(buf, &pi_string[i], cc);
		buf[cc] = '\0';
		fprintf(fd, "%s\n", buf);
		i += cc;
	}
	fclose(fd);
	assert(i == sz);
}

void make_pi(unsigned long digits, const char *pi_algo_name, const int pi_algo)
{
/*
	FILE *fd;
*/
	unsigned long last_k, curr_k = 0UL;
	struct mpfr_pi_impl *impl = NULL;
	mpfr_t term_sum, term_sum_plusplus, pi_value;
	/*
	 * timers stuff
	 */
	uint64_t time0, time1, time2;
        uint64_t tss3, tss4;
	char datebuf[128];
	char offsetbuf[128];
	char filename[256];
	unsigned long cc;

#if 0
	/*
	 * open results file right away, we don't want to compute for hour only to find out that
	 * this fails.
	 */
	snprintf(filename, sizeof (filename), "FPI_%ld_%s.txt", digits, algorithm);
	fd = fopen(filename, "w");
	assert(fd != NULL);
#endif

	printf("make_pi: digits = %lu, algorithm: %s\n", digits, pi_algo_name);

	mpfr_init2(term_sum, CFG_MPFR_PREC);
	mpfr_init2(term_sum_plusplus, CFG_MPFR_PREC);
	mpfr_init2(pi_value, CFG_MPFR_PREC);

	time0 = gettimestamp_nsecs();
	last_k = 0;

	ts_to_date_str(datebuf, sizeof (datebuf), time0);
	ts_to_offset_str(offsetbuf, sizeof (offsetbuf), time0 - time0);

	printf("%s: %s: make_pi, digits = %ld\n", datebuf, offsetbuf, digits);

	/*
	 * the extra iterations are not technically necessary, but just to be safe .....
	 */
	mpfr_set_ui(term_sum, 0UL, CFG_MPFR_PREC);
	last_k = curr_k;
	for (cc = 0; ; cc++) {

		unsigned long k_out;
		mpfr_t *curr_term_sum;

		if (impl == NULL) {
			assert(pi_algo == RAMANUJAN_1910_OPT_ALGO);
			printf("::::\n");
			impl = pi_impl_ramanujan_1910_opt_initialize(digits, curr_k);
			printf(":::::\n");
		}
		curr_term_sum = (mpfr_t *)(*impl->f_pi_compute_next_term)(impl, &k_out);
		assert(curr_k == k_out);
		curr_k++;

		// printf("ret=%d, curr_k=%lu, digits_out=%ld\n", ret, curr_k, curr_digits);

		tss4 = gettimestamp_nsecs();
		if (ts_secs_portion(tss4 - tss3) >= 60) {
			char *pi_value_s, *pi_value_s_plusplus;
			unsigned long c = 0;

			printf(":\n");

			mpfr_set(term_sum_plusplus, term_sum, CFG_MPFR_PREC);

			mpfr_add(term_sum, term_sum, (*curr_term_sum), CFG_MPFR_PREC);
			(*impl->f_pi_get_value)(&term_sum, &pi_value);
			pi_value_s = mpfr_t_to_str(&pi_value, digits + 1); /* the +1 accounts for the decimal point */

			printf(":1\n");

			curr_term_sum = (mpfr_t *)(*impl->f_pi_compute_next_term)(impl, &k_out);
			mpfr_add(term_sum_plusplus, term_sum_plusplus, (*curr_term_sum), CFG_MPFR_PREC);
			(*impl->f_pi_get_value)(&term_sum_plusplus, &pi_value);
			pi_value_s_plusplus = mpfr_t_to_str(&pi_value, digits + 1); /* the +1 accounts for the decimal point */

			printf(":2\n");

			ts_to_date_str(datebuf, sizeof (datebuf), gettimestamp_nsecs());
			ts_to_offset_str(offsetbuf, sizeof (offsetbuf), tss4 - time0);
			printf("%s: %s: k = %lu, k_delta = %lu\n", datebuf, offsetbuf, curr_k, (curr_k - last_k));
			tss3 = tss4;
			last_k = curr_k;

			for (c = 0; ; c++) {
				if (pi_value_s[c] == '\0' || pi_value_s_plusplus[c] == '\0') {
					printf("   pi_digits = %lu\n", c - 2);
					printf("   %s\n", pi_value_s);
					break;
				}
				if (pi_value_s[c] != pi_value_s_plusplus[c]) {
					printf("   pi_digits = %lu\n", c - 2);
					pi_value_s[c] = '\0';
					printf("   %s\n", pi_value_s);
					break;
				}
			}

			free_mpfr_str(pi_value_s);
			free_mpfr_str(pi_value_s_plusplus);
			(*impl->f_deinitialize)(impl);

			printf(":3\n");

			impl = NULL;

			printf("::\n");
		}
	}

	mpfr_clear(term_sum);
	mpfr_clear(term_sum_plusplus);
	mpfr_clear(pi_value);
}

int main(int argc, char **argv)
{
	unsigned long digits;

	setbuf(stdout, NULL);
	setbuf(stderr, NULL);
	setbuf(stdin, NULL);
	fclose(stdin);

	printf("MPFR library: %-12s\nMPFR header:  %s (based on %d.%d.%d)\n",
	       mpfr_get_version (), MPFR_VERSION_STRING, MPFR_VERSION_MAJOR,
	       MPFR_VERSION_MINOR, MPFR_VERSION_PATCHLEVEL);
	printf("MPFR_PREC_MAX = %ld\n", MPFR_PREC_MAX);
	printf("\n");
	printf("CFG_MPFR_PREC = %d\n", CFG_MPFR_PREC);
	printf("mpfr_custom_get_size(CFG_MPFR_PREC) = %ld\n", mpfr_custom_get_size(CFG_MPFR_PREC));
	printf("approximate decimals for CFG_MPFR_PREC = %ld (upper bound)\n", (long)((long)CFG_MPFR_PREC / 4));
	printf("\n");

	if (argc != 2) {
		printf("mpfr_pi: usage: mpfr_pi digits\n");
		exit(1);
	}

	digits = strtoul(argv[1], NULL, 0);
	if (digits <= 0) {
		printf("invalid %lu parameter for digits\n", digits);
		exit(1);
	}
	if (digits >= MAX_DIGITS_PREC) {
		printf("this build does not support %lu digits (max is %lu)\n", digits, MAX_DIGITS_PREC);
		exit(1);
	}

	printf("calculating pi to %lu digits (max digits is %lu) using %s algorithm\n", digits, MAX_DIGITS_PREC, argv[2]);

	make_pi(digits, "ramanujan_1910_opt", RAMANUJAN_1910_OPT_ALGO);

	return 0;
}
