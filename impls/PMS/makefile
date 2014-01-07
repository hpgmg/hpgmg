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

FFLAGS += -DHAVE_COMM_ARGS ${FC_FLAGS} -O3
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
	-@./pms.ex -problem 0 -nxloc 64 -nyloc 64 -nzloc 64 -nvcycles 0 -nfcycles 1 -verbose 1

debug:
	-@/usr/bin/gdb --args ./pms.ex -nxloc 64 -nyloc 64 -nzloc 64 -nvcycles 0 -nfcycles 1 -verbose 4

val:
	-@valgrind --dsymutil=yes --leak-check=no --gen-suppressions=yes --num-callers=20 ./pms.ex -nxloc 64 -nyloc 64 -nzloc 64 -nvcycles 0 -nfcycles 1 -verbose 4

runcray:
	-@aprun -n 8 ./pms.ex -nxloc 8 -nyloc 64 -nzloc 8 -nvcycles 0 -nfcycles 1

runmpi:
	-@${MPIEXEC} -n 8 ./pms.ex -nxloc 32 -nyloc 32 -nzloc 32 -nvcycles 10 -nfcycles 1 -verbose 1 -nsmooths 2
