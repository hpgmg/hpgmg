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
	-@${FLINKER} -o pms.ex ${OBJ} ${PETSC_SYS_LIB}
	${RM} ${OBJ}	

edison: ${OBJ} 
	-@${FLINKER}  -o pms.edison.ex ${OBJ} ${PETSC_SYS_LIB}
	${RM} ${OBJ}	

hopper: ${OBJ} 
	-@${FLINKER}  -o pms.hopper.ex ${OBJ} ${PETSC_SYS_LIB}
	${RM} ${OBJ}	

osx: ${OBJ} 
	-@${FLINKER} -o pms.ex ${OBJ} ${PETSC_SYS_LIB}
	${RM} ${OBJ}	

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
	-@${MPIEXEC} -n 8 valgrind --dsymutil=yes --leak-check=yes --gen-suppressions=yes --num-callers=20 ./pms.ex -nxloc 64 -nyloc 64 -nzloc 64 -nvcycles 0 -nfcycles 1 -verbose 2

runmpi:
	-@${MPIEXEC} -n 8 ./pms.ex -verbose 1

runv:
	-@${MPIEXEC} -n 64 ./pms.ex -verbose 1 -nxloc 16 -nyloc 16 -nzloc 16 -nvcycles 10 -nfcycles 0 
