hpgmg-fv-y.c += $(call thisdir, \
	level.c \
	operators.7pt.c \
	mg.c \
	solvers.c \
	hpgmg.c \
	)

hpgmg-fv-$(CONFIG_TIMER_MPI).c += $(call thisdir,timer.mpi.c)
hpgmg-fv-$(CONFIG_TIMER_OMP).c += $(call thisdir,timer.omp.c)

