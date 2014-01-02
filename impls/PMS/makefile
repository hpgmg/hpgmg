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

FFLAGS += -DHAVE_COMM_ARGS ${FC_FLAGS} -O0
pms: ${OBJ} 
	-${FLINKER} -o pms.ex ${OBJ}
#	-/opt/local/bin/gfortran -o pms.ex ${OBJ}
#	${RM} ${OBJ}	

#$(OBJ): modules.o

#%.o : %.F
#	-${FLINKER} -c -fPIC  -Wall -Wno-unused-variable -Wno-unused-dummy-argument -g -O0  -I/Users/markadams/Codes/petsc/include -I/Users/markadams/Codes/petsc/arch-macosx-gnu-g/include $<

realclean:
	rm -f *.o *.mod *.history

run:
	-@./pms.ex -problem 0 -nx 8 -ny 8 -nz 8

debug:
	-@/usr/bin/gdb --args ./pms.ex -problem 0 -nx 8 -ny 8 -nz 8

val:
	-@valgrind --dsymutil=yes --leak-check=no --gen-suppressions=yes --num-callers=20 ./pms.ex -problem 0 -nx 16 -ny 16 -nz 16 -nvcycles 0 -nfcycles 1

runedison:
	-@aprun -n 64 ./pms.ex -nxloc 64 -nyloc 64 -nzloc 64 -nvcycles 0 -nfcycles 1

runmpi:
	-@${MPIEXEC} -n 64 ./pms.ex -nxloc 8 -nyloc 8 -nzloc 8 -nvcycles 0 -nfcycles 1
