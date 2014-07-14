# This file is meant to be included from $(HPGMG_ARCH)/Makefile

.SECONDEXPANSION:		# to expand $$(@D)/.DIR
.SUFFIXES:	                # Clear .SUFFIXES because we don't use implicit rules
.DELETE_ON_ERROR:               # Delete likely-corrupt target file if rule fails

OBJDIR ?= obj
LIBDIR ?= lib
BINDIR ?= bin
INCDIR ?= include

# function to prefix directory that contains most recently-parsed
# makefile (current) if that derictory is not ./
thisdir = $(addprefix $(dir $(lastword $(MAKEFILE_LIST))),$(1))
incsubdirs = $(addsuffix /local.mk,$(call thisdir,$(1)))
srctoobj = $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(filter-out $(OBJDIR)/%,$(1)))

hpgmg-fe-y.c :=
hpgmg-fv-y.c :=

all : $(if $(CONFIG_FE),hpgmg-fe) $(if $(CONFIG_FV),hpgmg-fv)

# Recursively include files for all targets
include $(SRCDIR)/local.mk

#### Rules ####
ifeq ($(V),)
  quiet_HELP := "Use \"$(MAKE) V=1\" to see the verbose compile lines.\n"
  quiet = @printf $(quiet_HELP)$(eval quiet_HELP:=)"  %10s %s\n" "$1$2" "$@"; $($1)
else ifeq ($(V),0)		# Same, but do not print any help
  quiet = @printf "  %10s %s\n" "$1$2" "$@"; $($1)
else				# Show the full command line
  quiet = $($1)
endif

CONFIG_XLCOMPILER := $(if $(findstring IBM XL,$(shell $(HPGMG_CC) -qversion 2>/dev/null || true)),y,)

%.$(AR_LIB_SUFFIX) : | $$(@D)/.DIR
	$(call quiet,AR) $(AR_FLAGS) $@ $^
	$(call quiet,RANLIB) $@

# gcc/gfortran style dependency flags; these are set in petscvariables starting with petsc-3.5
C_DEPFLAGS ?= $(if $(CONFIG_XLCOMPILER),-qmakedep=gcc,-MMD -MP)

# GCC-style syntax for C99.  Use "make C99FLAGS=-qlanglvl=extc99" or similar
# on systems that use different syntax to specify C99.
C99FLAGS := $(if $(findstring c99,$(PCC_FLAGS) $(HPGMG_CFLAGS) $(CFLAGS)),,$(if $(CONFIG_XLCOMPILER),-qlanglvl=extc99,-std=c99))

HPGMG_COMPILE.c = $(call quiet,CC) -c $(C99FLAGS) $(HPGMG_CPPFLAGS) $(CPPFLAGS) $(HPGMG_CFLAGS) $(CFLAGS) $(C_DEPFLAGS)
HPGMG_LINK = $(call quiet,CCLD) $(HPGMG_CFLAGS) $(CFLAGS) $(HPGMG_LDFLAGS) $(LDFLAGS) -o $@
CC = $(HPGMG_CC)
CCLD = $(if $(PCC_LINKER),$(PCC_LINKER),$(HPGMG_CC))

hpgmg-fe = $(BINDIR)/hpgmg-fe
hpgmg-fe : $(hpgmg-fe)
hpgmg-fe-y.o := $(patsubst %.c,%.o,$(filter $(OBJDIR)/%,$(hpgmg-fe-y.c))) $(call srctoobj,$(hpgmg-fe-y.c))
$(hpgmg-fe) : $(hpgmg-fe-y.o) | $$(@D)/.DIR
	$(HPGMG_LINK) $^ $(HPGMG_LDLIBS) $(LDLIBS) $(PETSC_SNES_LIB)

hpgmg-fv = $(BINDIR)/hpgmg-fv
hpgmg-fv : $(hpgmg-fv)
hpgmg-fv-y.o := $(call srctoobj,$(hpgmg-fv-y.c))
$(hpgmg-fv-y.o) : CPPFLAGS += $(CONFIG_FV_CPPFLAGS)
$(hpgmg-fv) : $(hpgmg-fv-y.o) | $$(@D)/.DIR
	$(HPGMG_LINK) $^ $(HPGMG_LDLIBS) $(LDLIBS) -lm

$(OBJDIR)/%.o: $(OBJDIR)/%.c
	$(HPGMG_COMPILE.c) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $$(@D)/.DIR
	$(HPGMG_COMPILE.c) $< -o $@

test: test-fe
prove: prove-fe

test-fe : $(hpgmg-fe)
	make -C "$(SRCDIR)/finite-element/test" PETSC_DIR="$(PETSC_DIR)" PETSC_ARCH="$(PETSC_ARCH)" HPGMG_BINDIR="$(abspath $(BINDIR))" all

prove-fe : $(hpgmg-fe)
	make -C "$(SRCDIR)/finite-element/test" PETSC_DIR="$(PETSC_DIR)" PETSC_ARCH="$(PETSC_ARCH)" HPGMG_BINDIR="$(abspath $(BINDIR))" prove

%/.DIR :
	@mkdir -p $(@D)
	@touch $@

.PRECIOUS: %/.DIR

.PHONY: all clean print hpgmg-fe hpgmg-fv test test-fe

clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(BINDIR)

# make print VAR=the-variable
print:
	@echo $($(VAR))

srcs.c := $(hpgmg-fe-y.c) $(hpgmg-fv-y.c)
srcs.o := $(call srctoobj,$(srcs.c))
srcs.d := $(srcs.o:%.o=%.d)
# Tell make that srcs.d are all up to date.  Without this, the include
# below has quadratic complexity, taking more than one second for a
# do-nothing build of PETSc (much worse for larger projects)
$(srcs.d) : ;

-include $(srcs.d)
