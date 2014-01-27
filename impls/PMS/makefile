all:pms
#
MOBJ=modules.o
OBJ=main.o mg.o kernels_simple.o bdryexchange.o
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
pms: ${MOBJ} ${OBJ} 
	-@${FLINKER} -o pms.${PETSC_ARCH}.ex ${MOBJ} ${OBJ} ${PETSC_SYS_LIB}
#	${RM} ${OBJ} *.mod	

$(OBJ): ${MOBJ}

#%.o : %.F
#	-${FLINKER} -c -fPIC  -Wall -Wno-unused-variable -Wno-unused-dummy-argument -g -O0  -I/Users/markadams/Codes/petsc/include -I/Users/markadams/Codes/petsc/arch-macosx-gnu-g/include $<

realclean:
	rm -f *.o *.mod *.history

# regression test -- this is what should work
reg:
	-@${MPIEXEC} -n 8 ./pms.${PETSC_ARCH}.ex -verbose 2 -nxloc 128 -nyloc 128 -nzloc 128 -verbose 2 -show_options true -nsmoothsup 2 -nsmoothsdown 2

regsr:
	-@${MPIEXEC} -n 8 ./pms.${PETSC_ARCH}.ex -nxloc 4 -nyloc 4 -nzloc 4 -verbose 2 -pe_min_sz 4 -sr_max_loc_sz 8 -nsmoothsup 2 -nsmoothsdown 2

valreg:
	-@${MPIEXEC} -n 64 ~/Codes/valgrind-3.9.0/bin/bin/valgrind --dsymutil=yes --leak-check=no --gen-suppressions=no --num-callers=20 ./pms.${PETSC_ARCH}.ex -verbose 2 -nxloc 8 -nyloc 8 -nzloc 8

runv:
	-@${MPIEXEC} -n 8 ./pms.${PETSC_ARCH}.ex -verbose 2 -nxloc 128 -nyloc 128 -nzloc 128 -nvcycles 10 -nfcycles 0 -nsmoothsup 3

SRLOCN=8
SRFINEST=32
SRNP=8
SRBUF_BASE=8
SRBUF_INC=0
NSMUP=2
NSMDOWN=0

sr:
	-@${MPIEXEC} -n ${SRNP} ./pms.${PETSC_ARCH}.ex -nxloc ${SRLOCN} -nyloc ${SRLOCN} -nzloc ${SRLOCN} -verbose 3 -pe_min_sz ${SRLOCN} -sr_max_loc_sz ${SRFINEST} -show_options true -sr_base_bufsz ${SRBUF_BASE} -sr_bufsz_inc  ${SRBUF_INC} -nsmoothsup ${NSMUP} -nsmoothsdown ${NSMDOWN}

run:
	-@${MPIEXEC} -n ${SRNP} ./pms.${PETSC_ARCH}.ex -problem 0 -nxloc ${SRFINEST} -nyloc ${SRFINEST} -nzloc ${SRFINEST} -nvcycles 0 -nfcycles 1 -verbose 3 -pe_min_sz ${SRLOCN} -show_options true -nsmoothsup ${NSMUP} -nsmoothsdown ${NSMDOWN}

val:
	-@${MPIEXEC} -n ${SRNP} ~/Codes/valgrind-3.9.0/bin/bin/valgrind --dsymutil=yes --leak-check=no --gen-suppressions=no --num-callers=20 ./pms.${PETSC_ARCH}.ex -nxloc ${SRLOCN} -nyloc ${SRLOCN} -nzloc ${SRLOCN} -verbose 2 -pe_min_sz ${SRLOCN} -sr_max_loc_sz ${SRFINEST} -show_options true -sr_base_bufsz ${SRBUF_BASE} -sr_bufsz_inc ${SRBUF_INC} -nsmoothsup ${NSMUP} -nsmoothsdown ${NSMDOWN}

debug:
	-@/usr/bin/gdb --args ./pms.${PETSC_ARCH}.ex -nxloc ${SRLOCN} -nyloc ${SRLOCN} -nzloc ${SRLOCN} -verbose 3 -pe_min_sz ${SRLOCN} -sr_max_loc_sz ${SRFINEST} -show_options true -sr_base_bufsz ${SRBUF_BASE} -sr_bufsz_inc ${SRBUF_INC} -nsmoothsup ${NSMUP} -nsmoothsdown ${NSMDOWN}
