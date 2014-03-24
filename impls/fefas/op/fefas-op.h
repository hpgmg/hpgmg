#ifndef _fefas_op_h
#define _fefas_op_h

#include <petscdm.h>
#include "../fefas.h"

typedef struct Op_private *Op;

PetscErrorCode OpCreateFromOptions(MPI_Comm,Op*);
PetscErrorCode OpDestroy(Op*);
PetscErrorCode OpSetDof(Op,PetscInt);
PetscErrorCode OpGetDof(Op,PetscInt*);
PetscErrorCode OpSetFEDegree(Op,PetscInt);
PetscErrorCode OpGetFEDegree(Op,PetscInt*);
PetscErrorCode OpSetContext(Op,void*);
PetscErrorCode OpGetContext(Op,void*);
PetscErrorCode OpSetApply(Op,PetscErrorCode (*)(Op,DM,Vec,Vec));
PetscErrorCode OpSetPointwiseSolution(Op,PetscErrorCode (*)(Op,const PetscReal[],PetscScalar[]));
PetscErrorCode OpSetPointwiseForcing(Op,PetscErrorCode (*)(Op,const PetscReal[],PetscScalar[]));
typedef PetscErrorCode (*OpPointwiseElementFunction)(Op,PetscInt,PetscInt,const PetscScalar[],const PetscReal[],const PetscScalar[],PetscScalar[]);
PetscErrorCode OpSetPointwiseElement(Op,OpPointwiseElementFunction,PetscInt);
PetscErrorCode OpSetDestroy(Op,PetscErrorCode (*)(Op));
PetscErrorCode OpRegister(const char name[],PetscErrorCode (*f)(Op));
PetscErrorCode OpInitializePackage(void);
PetscErrorCode OpFinalizePackage(void);
PetscErrorCode OpDestroy(Op*);
PetscErrorCode OpApply(Op op,DM dm,Vec U,Vec F);
PetscErrorCode OpRestrictState(Op op,DM dm,Vec Uf,Vec Uc);
PetscErrorCode OpRestrictResidual(Op op,DM dm,Vec Uf,Vec Uc);
PetscErrorCode OpInterpolate(Op op,DM dm,Vec Uc,Vec Uf);
PetscErrorCode OpSolution(Op op,DM dm,Vec U);
PetscErrorCode OpForcing(Op op,DM dm,Vec F);
PetscErrorCode OpGetDiagonal(Op op,DM dm,Vec Diag);

#endif
