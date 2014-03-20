fefas-y.c += $(call thisdir, \
	fefas.c \
	fefas-test.c \
	fmg.c \
	grid.c \
	tensor.c \
	)

include $(call incsubdirs,op)
