HPGMG_ARCH := $(if $(PETSC_ARCH),$(PETSC_ARCH),build)

all :
	./configure --arch=$(HPGMG_ARCH)
	$(MAKE) -C $(HPGMG_ARCH)
	@echo "Build complete in $(HPGMG_ARCH).  Use make -C $(HPGMG_ARCH) test to test."

test : all
	$(MAKE) -C $(HPGMG_ARCH) test

clean :
	$(MAKE) -C $(HPGMG_ARCH) clean

.PHONY: all test clean
