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

extern struct mpfr_pi_impl *pi_impl_ramanujan_1910_initialize(const long digits, long *out_iterations);
extern struct mpfr_pi_impl *pi_impl_ramanujan_1910_opt_initialize(const long digits, long *out_iterations);

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

void make_pi(long digits, const char *algorithm)
{
	FILE *fd;
	unsigned long last_k, max_k;
	struct mpfr_pi_impl *impl = NULL;
	mpfr_t *pi_value;
	char *pi_value_s;
	long pi_value_digits;
	/*
	 * timers stuff
	 */
	uint64_t time0, time1, time2;
        uint64_t tss3, tss4;
	char datebuf[128];
	char offsetbuf[128];
	char filename[256];

	/*
	 * only implementation available for now
	 */
	if (strcmp(algorithm, "ramanujan_1910") == 0) {
		impl = pi_impl_ramanujan_1910_initialize(digits, &max_k);
		assert(impl != NULL);
	}
	if (strcmp(algorithm, "ramanujan_1910_opt") == 0) {
		impl = pi_impl_ramanujan_1910_opt_initialize(digits, &max_k);
		assert(impl != NULL);
	}
	if (impl == NULL) {
		printf("make_pi: unknon algorithm %s\n", algorithm);
		printf("make_pi: supported algorithms:\n");
		printf("                ramanujan_1910\n");
		printf("                ramanujan_1910_opt\n");
		exit(3);
	}

	/*
	 * open results file right away, we don't want to compute for hour only to find out that
	 * this fails.
	 */
	snprintf(filename, sizeof (filename), "FPI_%ld_%s.txt", digits, algorithm);
	fd = fopen(filename, "w");
	assert(fd != NULL);

	printf("make_pi: algorithm: %s\n", (*impl->f_impl_get_name)());

	time0 = gettimestamp_nsecs();
	last_k = 0;

	ts_to_date_str(datebuf, sizeof (datebuf), time0);
	ts_to_offset_str(offsetbuf, sizeof (offsetbuf), time0 - time0);

	printf("%s: %s: make_pi, digits = %ld, max_k = %lu\n", datebuf, offsetbuf, digits, max_k);

	/*
	 * the extra iterations are not technically necessary, but just to be safe .....
	 */
	for (;;) {

		unsigned long curr_k;
		long curr_digits;
		int ret = (*impl->f_pi_compute_next_term)(impl, &curr_k, &curr_digits);

		// printf("ret=%d, curr_k=%lu, digits_out=%ld\n", ret, curr_k, curr_digits);

		tss4 = gettimestamp_nsecs();
		if (ts_secs_portion(tss4 - tss3) >= 10) {
			ts_to_date_str(datebuf, sizeof (datebuf), gettimestamp_nsecs());
			ts_to_offset_str(offsetbuf, sizeof (offsetbuf), tss4 - time0);
			printf("%s: %s: k = %lu, k_delta = %lu, max_k = %lu\n", datebuf, offsetbuf, curr_k, (curr_k - last_k), max_k);
			tss3 = tss4;
			last_k = curr_k;
		}

		if (ret) {
			last_k = curr_k;
			break;
		}
	}

	pi_value = (*impl->f_pi_get_value)(impl, &pi_value_digits);
	assert(pi_value != NULL);

	tss4 = gettimestamp_nsecs();
	ts_to_date_str(datebuf, sizeof (datebuf), gettimestamp_nsecs());
	ts_to_offset_str(offsetbuf, sizeof (offsetbuf), tss4 - time0);
	printf("%s: %s: k = %lu, max_k = %lu, digits = %ld\n", datebuf, offsetbuf, last_k, max_k, pi_value_digits);
	tss3 = tss4;

	time1 = gettimestamp_nsecs();

	/*
	 * print PI.
	 * conversion from internal binary representation to decimal takes a long time.
	 */
	pi_value_s = mpfr_t_to_str(pi_value, digits + 1); /* the +1 accounts for the decimal point */

	time2 = gettimestamp_nsecs();
	ts_to_date_str(datebuf, sizeof (datebuf), time2);
	ts_to_offset_str(offsetbuf, sizeof (offsetbuf), time2 - tss3);
	printf("%s: %s: (finalization and conversion base 10)\n", datebuf, offsetbuf);
	writeout_pi(fd, pi_value_s);

	free_mpfr_str(pi_value_s);

	(*impl->f_deinitialize)(impl);

	printf("%s: %s: all done, output in %s\n", datebuf, offsetbuf, filename);
}

int main(int argc, char **argv)
{
	long digits;

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

	if (argc != 3) {
		printf("mpfr_pi: usage: mpfr_pi digits algorithm\n");
		exit(1);
	}

	digits = strtoul(argv[1], NULL, 0);
	if (digits <= 0) {
		printf("invalid %ld parameter for digits\n", digits);
		exit(1);
	}
	if ((long)digits >= (long)((float)CFG_MPFR_PREC / 3.5) - 100L) {
		printf("this build does not support %ld digits (max is %ld)\n", digits, (long)((float)CFG_MPFR_PREC / 3.5) - 100L);
		exit(1);
	}

	printf("calculating pi to %ld digits using %s algorithm\n", digits, argv[2]);

	make_pi(digits, argv[2]);

	return 0;
}
