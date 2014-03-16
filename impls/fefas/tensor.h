#ifndef _tensor_h
#define _tensor_h

#include <petscsys.h>

typedef enum {TENSOR_EVAL,TENSOR_TRANSPOSE} TensorMode;

PetscErrorCode TensorContract(PetscInt ne,PetscInt dof,PetscInt P,PetscInt Q,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]);

#endif
