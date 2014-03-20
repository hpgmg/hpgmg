#ifndef _fefas_op_h
#define _fefas_op_h

#include <petscdm.h>

typedef struct Op_private *Op;

PetscErrorCode OpCreateFromOptions(MPI_Comm,Op*);
PetscErrorCode OpDestroy(Op*);
PetscErrorCode OpSetDof(Op,PetscInt);
PetscErrorCode OpSetFEDegree(Op,PetscInt);
PetscErrorCode OpRegister(const char name[],PetscErrorCode (*f)(Op));
PetscErrorCode OpInitializePackage(void);
PetscErrorCode OpFinalizePackage(void);
PetscErrorCode OpDestroy(Op*);

#endif
