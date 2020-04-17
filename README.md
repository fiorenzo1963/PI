# PI
* Compute PI with arbitrary precision using MPFR arbitrary precision floating point library, using various algorithms.
* At the moment only Ramanuja's 1910 algorithm is used.
* The precomputed PI digits in here are taken from publicly available sources, and used to compare algorithm accuracy.
* Currently only serialized computation is supported, see single_process/ directory.

# Build
* make sure libgmp and libmpfr are installed. if not, install the packages. if the package is not available you can download the source code and build them.
* go to sinle_process, and use make to build. by default it will link with static libraries, so the executable can be copied to another machine

# Run

* ./mpfr_pi <number_of_desired_digits> <algorithm>
* example:
```
	./mpfr_pi 1000 ramanujan_1910_opt
```
* sample output:
```
[fcattane@linux-oel77 single_process]$ ./mpfr_pi 100 ramanujan_1910_opt
MPFR library: 4.0.1-p13   
MPFR header:  4.0.1-p13 (based on 4.0.1)
MPFR_PREC_MAX = 9223372036854775551

CFG_MPFR_PREC = 10000000
mpfr_custom_get_size(CFG_MPFR_PREC) = 1250000
approximate decimals for CFG_MPFR_PREC = 2500000 (upper bound)

calculating pi to 100 digits using ramanujan_1910_opt algorithm
pi_impl_ramanujan_1910_opt_initialize: desired digits = 100
pi_impl_ramanujan_1910_opt_initialize: max_k = 16
make_pi: algorithm: Ramanujan 1910 Formula (optimized)
2020-04-16:20:24:57.552808: 0:00:00.000000: make_pi, digits = 100, max_k = 16
2020-04-16:20:24:58.176462: 0:00:00.623653: k = 0, k_delta = 0, max_k = 16
2020-04-16:20:25:01.031660: 0:00:03.478848: k = 16, max_k = 16, digits = 125
2020-04-16:20:25:03.110989: 0:00:02.079331: (finalization and conversion base 10)
2020-04-16:20:25:03.110989: 0:00:02.079331: all done, output in FPI_100_ramanujan_1910_opt.txt
[fcattane@linux-oel77 single_process]$ 
```
* output is places in the file shown in the output timing statistics. in the above example:
```
[fcattane@linux-oel77 single_process]$ cat FPI_100_ramanujan_1910_opt.txt 
3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706
[fcattane@linux-oel77 single_process]$ 
```

# BUGS
After a certain number of requested digits, precision might be a litte less. For instance when specifying 100000 digits, you may actually a few hundred digits less precision. This will be fixed soon.
