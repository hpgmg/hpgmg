# This file is meant to be included from $(HPGMG_ARCH)/Makefile

.SECONDEXPANSION:		# to expand $$(@D)/.DIR
.SUFFIXES:	                # Clear .SUFFIXES because we don't use implicit rules
.DELETE_ON_ERROR:               # Delete likely-corrupt target file if rule fails

OBJDIR ?= obj
LIBDIR ?= lib
BINDIR ?= bin

# function to prefix directory that contains most recently-parsed
# makefile (current) if that derictory is not ./
thisdir = $(addprefix $(dir $(lastword $(MAKEFILE_LIST))),$(1))
incsubdirs = $(addsuffix /local.mk,$(call thisdir,$(1)))

fefas-y.c :=

all : fefas

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
ifeq ($(PETSC_LANGUAGE),CXXONLY)
  cc_name := CXX
else
  cc_name := CC
endif

# Detect HBM performance counting on ALCF Blue Gene/Q machines
CONFIG_HBM := $(shell test -f /soft/perftools/hpctw/lib/libmpihpm.a -a -f /bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a && echo y || true)
ifeq ($(CONFIG_HBM),y)
  CCPPFLAGS += -DCONFIG_HBM
  FEFAS_LDLIBS += /soft/perftools/hpctw/lib/libmpihpm.a /bgsys/drivers/ppcfloor/bgpm/lib/libbgpm.a
endif

CONFIG_XLCOMPILER := $(if $(findstring IBM XL,$(shell $(CC) -qversion 2>/dev/null || true)),y,)

%.$(AR_LIB_SUFFIX) : | $$(@D)/.DIR
	$(call quiet,AR) $(AR_FLAGS) $@ $^
	$(call quiet,RANLIB) $@

# gcc/gfortran style dependency flags; these are set in petscvariables starting with petsc-3.5
C_DEPFLAGS ?= -MMD -MP

# GCC-style syntax for C99.  Use "make C99FLAGS=-qlanglvl=extc99" or similar
# on systems that use different syntax to specify C99.
C99FLAGS := $(if $(findstring c99,$(PCC_FLAGS) $(CFLAGS)),,$(if $(CONFIG_XLCOMPILER),-qlanglvl=extc99,-std=c99))

FEFAS_COMPILE.c = $(call quiet,$(cc_name)) -c $(C99FLAGS) $(PCC_FLAGS) $(CCPPFLAGS) $(CFLAGS) $(C_DEPFLAGS)

fefas = $(BINDIR)/fefas
fefas : $(fefas)
fefas-y.o := $(patsubst %.c,%.o,$(filter $(OBJDIR)/%,$(fefas-y.c))) $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(filter-out $(OBJDIR)/%,$(fefas-y.c)))
$(BINDIR)/fefas : $(fefas-y.o) | $$(@D)/.DIR
	$(call quiet,CLINKER) -o $@ $^ $(LDLIBS) $(FEFAS_LDLIBS) $(PETSC_SNES_LIB) $(LIBZ_LIB)

$(OBJDIR)/%.o: $(OBJDIR)/%.c
	$(FEFAS_COMPILE.c) $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $$(@D)/.DIR
	$(FEFAS_COMPILE.c) $< -o $@

test: test-fe

test-fe : $(fefas)
	make -C "$(SRCDIR)/finite-element/test" PETSC_DIR="$(PETSC_DIR)" PETSC_ARCH="$(PETSC_ARCH)" FEFAS_BINDIR="$(abspath $(BINDIR))" all

%/.DIR :
	@mkdir -p $(@D)
	@touch $@

.PRECIOUS: %/.DIR

.PHONY: all clean print fefas test test-fe

clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(BINDIR)

# make print VAR=the-variable
print:
	@echo $($(VAR))

srcs.c := $(fefas-y.c)
srcs.o := $(srcs.c:%.c=$(OBJDIR)/%.o)
srcs.d := $(srcs.o:%.o=%.d)
# Tell make that srcs.d are all up to date.  Without this, the include
# below has quadratic complexity, taking more than one second for a
# do-nothing build of PETSc (much worse for larger projects)
$(srcs.d) : ;

-include $(srcs.d)
