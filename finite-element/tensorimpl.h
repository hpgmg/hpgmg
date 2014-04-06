#ifndef _tensorimpl_h
#define _tensorimpl_h

#include "tensor.h"

struct Tensor_private {
  PetscInt ne;
  PetscInt dof;
  PetscInt P;
  PetscInt Q;
  PetscErrorCode (*Contract)(Tensor,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]);
};

#endif
