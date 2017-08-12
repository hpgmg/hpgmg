#ifndef _fefas_op_h
#define _fefas_op_h

#include <petscdm.h>
#include "../tensor.h"

typedef struct Op_private *Op;

MPI_Comm OpComm(Op);
PetscErrorCode OpCreateFromOptions(MPI_Comm,Op*);
PetscErrorCode OpDestroy(Op*);
PetscErrorCode OpSetDof(Op,PetscInt);
PetscErrorCode OpGetDof(Op,PetscInt*);
PetscErrorCode OpSetFEDegree(Op,PetscInt);
PetscErrorCode OpGetFEDegree(Op,PetscInt*);
PetscErrorCode OpSetAddQuadPts(Op,PetscInt);
PetscErrorCode OpGetAddQuadPts(Op,PetscInt*);
PetscErrorCode OpSetContext(Op,void*);
PetscErrorCode OpGetContext(Op,void*);
PetscErrorCode OpSetApply(Op,PetscErrorCode (*)(Op,DM,Vec,Vec));
PetscErrorCode OpSetPointwiseSolution(Op,PetscErrorCode (*)(Op,const PetscReal[],const PetscReal[],PetscScalar[]));
PetscErrorCode OpSetPointwiseForcing(Op,PetscErrorCode (*)(Op,const PetscReal[],const PetscReal[],PetscScalar[]));
typedef PetscErrorCode (*OpPointwiseElementFunction)(Op,PetscInt,PetscInt,const PetscScalar[],const PetscReal[],const PetscScalar[],const PetscScalar[],PetscScalar[],PetscScalar[]);
PetscErrorCode OpSetPointwiseElement(Op,OpPointwiseElementFunction,PetscInt,PetscInt);
PetscErrorCode OpSetAffineOnly(Op op,PetscBool affine);
PetscErrorCode OpGetAffineOnly(Op op,PetscBool *affine);
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
PetscErrorCode OpIntegrateNorms(Op op,DM dm,Vec U,PetscReal *normInfty,PetscReal *norm2);
PetscErrorCode OpGetDiagonal(Op op,DM dm,Vec Diag);
PetscErrorCode OpGetMat(Op op,DM dm,Mat *shell);
PetscErrorCode OpGetTensors(Op op,Tensor *TensorDOF,Tensor *Tensor3);

#endif
