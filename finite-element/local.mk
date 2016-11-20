hpgmg-fe-y.c += $(call thisdir, \
	fefas.c \
	fefas-test.c \
	fmg.c \
	grid.c \
	memusage.c \
	sampler.c \
	tensor.c \
	tensor-fma.c \
	tensor-avx512.c \
	tensor-qpx.c \
	)

include $(call incsubdirs,op)
