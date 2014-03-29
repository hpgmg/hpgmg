#ifndef _tensor_h
#define _tensor_h

#include <petscsys.h>

#ifdef __AVX__
#  define _align __attribute__((aligned(32))) /* AVX packed instructions need 32-byte alignment */
#else
#  define _align
#endif

typedef enum {TENSOR_EVAL,TENSOR_TRANSPOSE} TensorMode;

typedef struct Tensor_private *Tensor;

PetscErrorCode TensorCreate(PetscInt ne,PetscInt dof,PetscInt P,PetscInt Q,Tensor *ten);
PetscErrorCode TensorDestroy(Tensor *ten);
PetscErrorCode TensorContract(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]);

#endif
