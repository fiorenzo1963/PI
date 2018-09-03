#!/bin/bash
DIGITS="\
	1000 \
	10000 \
	100000 \
	200000 \
	500000 \
	1000000 \
	2000000 \
"
for d in $DIGITS
do
	echo $0: `date`: begin: ./mpfr_pi $d ramananujan_1910_opt
	./mpfr_pi $d ramananujan_1910_opt
	echo $0: `date`: end: ./mpfr_pi $d ramananujan_1910_opt
done
