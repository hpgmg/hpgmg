#ifndef _tensor_h
#define _tensor_h

#include <petscsys.h>
#include "fefas-align.h"

typedef enum {TENSOR_EVAL,TENSOR_TRANSPOSE} TensorMode;

typedef struct Tensor_private *Tensor;

PetscErrorCode TensorCreate(PetscInt ne,PetscInt dof,PetscInt P,PetscInt Q,Tensor *ten);
PetscErrorCode TensorDestroy(Tensor *ten);
PetscErrorCode TensorContract(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]);

#endif
