FILES_H := mpfr_pi_generic.h stringify.h subr.h
FILES_C_IMPL := mpfr_pi_impl_ramanujan_1910.c mpfr_pi_impl_ramanujan_1910_opt.c
FILES_C := mpfr_pi.c subr.c
OPT := -O3
LOCAL_H := -I/usr/local/include
LOCAL_LIB_PATH := /usr/local/lib
LOCAL_LIB_LD := -L$(LOCAL_LIB_PATH)

#mpfr_pi: $(FILES_H) $(FILES_C) $(FILES_C_IMPL)
#	cc $(OPT) -o mpfr_pi $(FILES_C) $(FILES_C_IMPL) -lmpfr -lgmp

mpfr_pi: $(FILES_H) $(FILES_C) $(FILES_C_IMPL)
	cc $(OPT) $(LOCAL_H) -o mpfr_pi $(FILES_C) $(FILES_C_IMPL) $(LOCAL_LIB_PATH)/libmpfr.a $(LOCAL_LIB_PATH)/libgmp.a

clean:
	rm -f mpfr_pi mpfr_pi.x *.o core *.log *.out FPI*txt
