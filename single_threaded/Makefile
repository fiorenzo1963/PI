FILES_H := mpfr_pi_generic.h stringify.h subr.h
FILES_C_IMPL := mpfr_pi_impl_ramanujan_1910.c mpfr_pi_impl_ramanujan_1910_opt.c
FILES_C := mpfr_pi.c subr.c

mpfr_pi: $(FILES_H) $(FILES_C) $(FILES_C_IMPL)
	cc -o mpfr_pi $(FILES_C) $(FILES_C_IMPL) -lmpfr -lgmp

clean:
	rm -f mpfr_pi *.o core *.log *.out FPI*txt
