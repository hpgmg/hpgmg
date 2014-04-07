hpgmg-fv-y.c += $(call thisdir, \
	level.c \
	operators.7pt.c \
	mg.c \
	solvers.c \
	hpgmg.c \
	)

hpgmg-fv-$(CONFIG_TIMER_X86).c += $(call thisdir,timer.x86.c)
hpgmg-fv-$(CONFIG_TIMER_BGQ).c += $(call thisdir,timer.bgq.c)
hpgmg-fv-$(CONFIG_TIMER_OMP).c += $(call thisdir,timer.omp.c)
hpgmg-fv-$(CONFIG_TIMER_MPI).c += $(call thisdir,timer.mpi.c)
hpgmg-fv-$(CONFIG_TIMER_SPARC).c += $(call thisdir,timer.sparc.c)
