#INCLUDES = syspara.h

bin_PROGRAMS = 1D_mORD_LP
1D_mORD_LP_CFLAGS = -O2 -axsse4.2 -mkl -parallel -Winline -inline-level=2 -ip -opt-mem-layout-trans=3 -pg -par-report2

1D_mORD_LP_INCLUDE = -I./include
1D_mORD_LP_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

1D_mORD_LP_SOURCES = ./src/ord.c \
	./src/exptable.c ./src/linear_mat.c ./src/pardiso_mat.c\
	./src/eular.c ./src/input_para.c ./src/system.c ./src/data_out.c\
	./src/mem.c ./src/var_consts.c ./src/syspara.h

