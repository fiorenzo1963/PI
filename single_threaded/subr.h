#ifndef _SUBR_H_
#define _SUBR_H_

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>

extern uint64_t gettimestamp_nsecs(void);

static inline uint64_t nsecs_to_secs(uint64_t ns)
{
	return ns / (1000LL * 1000LL * 1000LL);
}
static inline uint64_t ts_secs_portion(uint64_t ts)
{
	return nsecs_to_secs(ts);
}
static inline uint64_t ts_secs_nsecs_portion(uint64_t ts)
{
	return ts % (1000LL * 1000LL * 1000LL);
}
static inline uint64_t ts_secs_usecs_portion(uint64_t ts)
{
	return ts_secs_nsecs_portion(ts) / 1000LL;
}
static inline uint64_t ts_secs_msecs_portion(uint64_t ts)
{
	return ts_secs_usecs_portion(ts) / 1000LL;
}

extern void ts_to_offset_str(char *buf, size_t sz, uint64_t ts);
extern void ts_to_date_str(char *buf, size_t sz, uint64_t ts);

#endif
