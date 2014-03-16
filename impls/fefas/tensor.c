#include "tensor.h"

PetscErrorCode TensorContract(PetscInt ne,PetscInt dof,PetscInt P,PetscInt Q,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[])
{

  PetscFunctionBegin;
  if (tmode == TENSOR_TRANSPOSE) {PetscInt tmp = Q; Q = P; P = tmp;}
  {
    PetscReal R[Q][P],S[Q][P],T[Q][P];
    const PetscScalar (*x)[P*P*P][ne] = (const PetscScalar(*)[P*P*P][ne])xx;
    PetscScalar       (*y)[P*P*P][ne] =       (PetscScalar(*)[Q*Q*Q][ne])yy;
    PetscScalar u[dof][Q*P*P][ne],v[dof][Q*Q*P][ne];

    for (PetscInt i=0; i<Q; i++) {
      for (PetscInt j=0; j<P; j++) {
        R[i][j] = tmode == TENSOR_EVAL ? Rf[i*P+j] : Rf[j*Q+i];
        S[i][j] = tmode == TENSOR_EVAL ? Sf[i*P+j] : Sf[j*Q+i];
        T[i][j] = tmode == TENSOR_EVAL ? Tf[i*P+j] : Tf[j*Q+i];
      }
    }

    // u[l,a,j,k] = R[a,i] x[l,i,j,k]
    PetscMemzero(u,sizeof u);
    for (PetscInt i=0; i<P; i++) {
      for (PetscInt l=0; l<dof; l++) {
        for (PetscInt a=0; a<Q; a++) {
          for (PetscInt jk=0; jk<P*P; jk++) {
            for (PetscInt e=0; e<ne; e++) u[l][a*P*P+jk][e] += R[a][i] * x[l][i*P*P+jk][e];
          }
        }
      }
    }

    // v[l,a,b,k] = S[b,j] u[l,a,j,k]
    PetscMemzero(v,sizeof v);
    for (PetscInt l=0; l<dof; l++) {
      for (PetscInt a=0; a<Q; a++) {
        for (PetscInt k=0; k<P; k++) {
          for (PetscInt j=0; j<P; j++) {
            for (PetscInt b=0; b<Q; b++) {
              for (PetscInt e=0; e<ne; e++) v[l][(a*Q+b)*P+k][e] += S[b][j] * u[l][(a*P+j)*P+k][e];
            }
          }
        }
      }
    }

    // y[l,a,b,c] = T[c,k] v[l,a,b,k]
    for (PetscInt l=0; l<dof; l++) {
      for (PetscInt ab=0; ab<Q*Q; ab++) {
        for (PetscInt k=0; k<P; k++) {
          for (PetscInt c=0; c<Q; c++) {
            for (PetscInt e=0; e<ne; e++) y[l][ab*Q+c][e] += T[c][k] * v[l][ab*P+k][e];
          }
        }
      }
    }
    PetscLogFlops(dof*(Q*P*P*P+Q*Q*P*P+Q*Q*Q*P)*ne*2);
  }
  PetscFunctionReturn(0);
}
