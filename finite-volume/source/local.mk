hpgmg-fv-y.c += $(call thisdir, \
	level.c \
	operators.7pt.c \
	mg.c \
	solvers.c \
	hpgmg.c \
	)

hpgmg-fv-$(CONFIG_X86).c += $(call thisdir,timer.x86.c)
hpgmg-fv-$(CONFIG_BGQ).c += $(call thisdir,timer.bgq.c)

