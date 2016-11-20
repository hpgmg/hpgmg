#include "tensorimpl.h"

static inline PetscErrorCode TensorContract_Inline(PetscInt ne,PetscInt dof,PetscInt P,PetscInt Q,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[])
{

  PetscFunctionBegin;
  if (tmode == TENSOR_TRANSPOSE) {PetscInt tmp = Q; Q = P; P = tmp;}
  {
    PetscReal R[Q][P],S[Q][P],T[Q][P];
    const PetscScalar (*restrict x)[P*P*P][ne]_align = (const PetscScalar(*)[P*P*P][ne])xx;
    PetscScalar       (*restrict y)[P*P*P][ne]_align =       (PetscScalar(*)[Q*Q*Q][ne])yy;
    PetscScalar u[dof][Q*P*P][ne]_align,v[dof][Q*Q*P][ne]_align;

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

static inline PetscErrorCode TensorContract_Ref(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  PetscInt ne = ten->ne,dof = ten->dof,P = ten->P,Q = ten->Q;
  return TensorContract_Inline(ne,dof,P,Q,Rf,Sf,Tf,tmode,xx,yy);
}

static PetscErrorCode TensorContract_Ref_4_1_2_2(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  return TensorContract_Inline(4,1,2,2,Rf,Sf,Tf,tmode,xx,yy);
}
static PetscErrorCode TensorContract_Ref_4_3_2_2(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  return TensorContract_Inline(4,3,2,2,Rf,Sf,Tf,tmode,xx,yy);
}
static PetscErrorCode TensorContract_Ref_4_1_3_3(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  return TensorContract_Inline(4,1,3,3,Rf,Sf,Tf,tmode,xx,yy);
}
static PetscErrorCode TensorContract_Ref_4_3_3_3(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  return TensorContract_Inline(4,3,3,3,Rf,Sf,Tf,tmode,xx,yy);
}

PetscErrorCode TensorSelect_AVX(Tensor);
PetscErrorCode TensorSelect_AVX512(Tensor);
PetscErrorCode TensorSelect_QPX(Tensor);

PetscErrorCode TensorCreate(PetscInt ne,PetscInt dof,PetscInt P,PetscInt Q,Tensor *ten) {
  Tensor t;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNew(&t);CHKERRQ(ierr);
  t->ne  = ne;
  t->dof = dof;
  t->P   = P;
  t->Q   = Q;
  t->Contract = TensorContract_Ref;
  if (ne == 4) {
    switch (dof) {
    case 1: // Scalar problems with Q1 or Q2 elements
      if (P == 2 && Q == 2)      t->Contract = TensorContract_Ref_4_1_2_2;
      else if (P == 3 && Q == 3) t->Contract = TensorContract_Ref_4_1_3_3;
      break;
    case 3: // Coordinates or elasticity
      if (P == 2 && Q == 2)      t->Contract = TensorContract_Ref_4_3_2_2;
      else if (P == 3 && Q == 3) t->Contract = TensorContract_Ref_4_3_3_3;
      break;
    }
  }
  ierr = TensorSelect_AVX(t);CHKERRQ(ierr);
  ierr = TensorSelect_AVX512(t);CHKERRQ(ierr);
  ierr = TensorSelect_QPX(t);CHKERRQ(ierr);
  *ten = t;
  PetscFunctionReturn(0);
}

PetscErrorCode TensorDestroy(Tensor *ten) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFree(*ten);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TensorContract(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  return (*ten->Contract)(ten,Rf,Sf,Tf,tmode,xx,yy);
}
