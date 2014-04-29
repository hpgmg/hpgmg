op-impls.c := $(wildcard $(call thisdir,op-*.c))
genregister := $(call thisdir,genregister.py)
register.c := $(OBJDIR)/register.c

hpgmg-fe-y.c += $(call thisdir, \
	op.c \
	) $(op-impls.c) $(register.c)

$(register.c) : $(genregister) $(op-impls.c) | $$(@D)/.DIR
	$(PYTHON) $(genregister) $@ $(op-impls.c)

HPGMG_FE_DIR := $(call thisdir,..)
$(OBJDIR)/register.o : HPGMG_CPPFLAGS += -I$(HPGMG_FE_DIR)
