HPGMG_ARCH := $(if $(PETSC_ARCH),$(PETSC_ARCH),build)

all : config
	$(MAKE) -C $(HPGMG_ARCH)
	@echo "Build complete in $(HPGMG_ARCH).  Use 'make -C $(HPGMG_ARCH) test' to test."

config : $(HPGMG_ARCH)/Makefile

$(HPGMG_ARCH)/Makefile :
	./configure --arch=$(HPGMG_ARCH)

test : all
	$(MAKE) -C $(HPGMG_ARCH) test

clean :
	$(MAKE) -C $(HPGMG_ARCH) clean

.PHONY: all test clean config
