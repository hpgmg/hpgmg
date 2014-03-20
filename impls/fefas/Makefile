.SECONDEXPANSION:		# to expand $$(@D)/.DIR
.SUFFIXES:	                # Clear .SUFFIXES because we don't use implicit rules
.DELETE_ON_ERROR:               # Delete likely-corrupt target file if rule fails

include $(PETSC_DIR)/conf/variables

OBJDIR ?= $(PETSC_ARCH)/obj
LIBDIR ?= $(PETSC_ARCH)/lib
BINDIR ?= $(PETSC_ARCH)/bin

# function to prefix directory that contains most recently-parsed
# makefile (current) if that derictory is not ./
thisdir = $(subst ./,,$(addprefix $(dir $(lastword $(MAKEFILE_LIST))),$(1)))
incsubdirs = $(addsuffix /local.mk,$(call thisdir,$(1)))

libfefas-y.c :=
fefas-y.c :=

all : fefas

# Recursively include files for all targets
include local.mk

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

libfefas = $(LIBDIR)/libfefas.$(AR_LIB_SUFFIX)
libfefas : $(libfefas)
$(libfefas) : $(libfefas-y.c:%.c=$(OBJDIR)/%.o)

%.$(AR_LIB_SUFFIX) : | $$(@D)/.DIR
	$(call quiet,AR) $(AR_FLAGS) $@ $^
	$(call quiet,RANLIB) $@

# gcc/gfortran style dependency flags; these are set in petscvariables starting with petsc-3.5
C_DEPFLAGS ?= -MMD -MP

C99FLAGS := $(if $(filter c99,$(PCC_FLAGS)),,-std=c99)
FEFAS_COMPILE.c = $(call quiet,$(cc_name)) -c $(C99FLAGS) $(PCC_FLAGS) $(CCPPFLAGS) $(CFLAGS) $(C_DEPFLAGS)

fefas = $(BINDIR)/fefas
fefas : $(fefas)
fefas-y.o := $(patsubst %.c,%.o,$(filter $(OBJDIR)/%,$(fefas-y.c))) $(patsubst %.c,$(OBJDIR)/%.o,$(filter-out $(OBJDIR)/%,$(fefas-y.c)))
$(BINDIR)/fefas : $(fefas-y.o) $(libfefas) | $$(@D)/.DIR
	$(call quiet,CLINKER) -o $@ $^ $(PETSC_SNES_LIB) $(LIBZ_LIB)

$(OBJDIR)/%.o: $(OBJDIR)/%.c
	$(FEFAS_COMPILE.c) -I. $(abspath $<) -o $@

$(OBJDIR)/%.o: %.c | $$(@D)/.DIR
	$(FEFAS_COMPILE.c) $(abspath $<) -o $@

test : $(fefas)
	make -C test PETSC_ARCH=${PETSC_ARCH} all

%/.DIR :
	@mkdir -p $(@D)
	@touch $@

.PRECIOUS: %/.DIR

.PHONY: all clean print fefas libfefas test

clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(BINDIR)

# make print VAR=the-variable
print:
	@echo $($(VAR))

srcs.c := $(libfefas-y.c) $(fefas-y.c)
srcs.o := $(srcs.c:%.c=$(OBJDIR)/%.o)
srcs.d := $(srcs.o:%.o=%.d)
# Tell make that srcs.d are all up to date.  Without this, the include
# below has quadratic complexity, taking more than one second for a
# do-nothing build of PETSc (much worse for larger projects)
$(srcs.d) : ;

-include $(srcs.d)
