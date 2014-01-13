all:pms
#
OBJ= modules.o main.o kernels_simple.o bdryexchange.o mg.o
# f2kcli.o
#
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
#

#FFFLAGS= -DDYNAMIC -DZPERIODIC -DTWO_D -DPARALLEL -DXPERIODIC
#FFLAGS= $(FFFLAGS) -target=linux -fast -fastsse -O3
#FFLAGS= --warn --chk[aesu] --chkglobal -g --trace --trap $(FFFLAGS)
#F90=${FLINKER} 

#.SUFFIXES: .o .F

FFLAGS += -DHAVE_COMM_ARGS ${FC_FLAGS} -DHAVE_PETSC 
pms: ${OBJ} 
	-@${FLINKER} -o pms.${PETSC_ARCH}.ex ${OBJ} ${PETSC_SYS_LIB}
#	${RM} ${OBJ} *.mod	

#$(OBJ): modules.o

#%.o : %.F
#	-${FLINKER} -c -fPIC  -Wall -Wno-unused-variable -Wno-unused-dummy-argument -g -O0  -I/Users/markadams/Codes/petsc/include -I/Users/markadams/Codes/petsc/arch-macosx-gnu-g/include $<

realclean:
	rm -f *.o *.mod *.history

run1:
	-@./pms.${PETSC_ARCH}.ex -problem 0 -nxloc 64 -nyloc 64 -nzloc 64 -nvcycles 0 -nfcycles 1 -verbose 1

debug:
	-@/usr/bin/gdb --args ./pms.${PETSC_ARCH}.ex -nxloc 64 -nyloc 64 -nzloc 64 -nvcycles 0 -nfcycles 1 -verbose 4

val:
	-@${MPIEXEC} -n 1 valgrind --dsymutil=yes --leak-check=yes --gen-suppressions=yes --num-callers=20 ./pms.${PETSC_ARCH}.ex -nxloc 16 -nyloc 16 -nzloc 16 -sr_max_loc_sz 32 

run:
	-@${MPIEXEC} -n 64 ./pms.${PETSC_ARCH}.ex -verbose 1 -nxloc 16 -nyloc 16 -nzloc 16

runsr:
	-@${MPIEXEC} -n 64 ./pms.${PETSC_ARCH}.ex -verbose 1 -nxloc 16 -nyloc 16 -nzloc 16 -sr_max_loc_sz 32

runv:
	-@${MPIEXEC} -n 8 ./pms.${PETSC_ARCH}.ex -verbose 1 -nxloc 16 -nyloc 16 -nzloc 16 -nvcycles 10 -nfcycles 0 

runssr:
	-@${MPIEXEC} -n 64 ./pms.${PETSC_ARCH}.ex -verbose 1 -nxloc 16 -nyloc 16 -nzloc 16 -sr_min_sz 8
