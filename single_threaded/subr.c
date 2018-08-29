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

uint64_t gettimestamp_nsecs(void)
{
	struct timespec ts;
	uint64_t r;
	clock_gettime(CLOCK_REALTIME, &ts);
	r = ts.tv_sec * (1000LL * 1000LL * 1000LL);
	r += ts.tv_nsec;
	return r;
}

void ts_to_offset_str(char *buf, size_t sz, uint64_t ts)
{
	unsigned long hours = ts_secs_portion(ts) / 3600LL;
	unsigned long mins = (ts_secs_portion(ts) - (hours * 3600LL)) / 60LL;
	unsigned long secs = (ts_secs_portion(ts) - (hours * 3600LL)) % 60LL;
	snprintf(buf, sz, "%lu:%02lu:%02lu.%06lu", hours, mins, secs, (unsigned long)ts_secs_usecs_portion(ts));
}

void ts_to_date_str(char *buf, size_t sz, uint64_t ts)
{
	struct tm tm;
	size_t _sz;
	time_t t = (time_t)ts_secs_portion(ts);
	localtime_r(&t, &tm);
	_sz = strftime(buf, sz, "%F:%T", &tm);
	if (_sz < sz)
		snprintf(&buf[_sz], sz - _sz, ".%06llu", (long long unsigned int)ts_secs_usecs_portion(ts));
}
