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

double gettimestamp_secs(void)
{
	struct timespec ts;
	double r;
	clock_gettime(CLOCK_REALTIME, &ts);
	r = (double)ts.tv_sec;
	r += (double)ts.tv_nsec / (1000.0 * 1000.0 * 1000.0);
	return r;
}
