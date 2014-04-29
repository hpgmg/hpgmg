#include "fefas-op.h"
#include "../fefas.h"
#include "../tensor.h"
#include "../pointwise.h"
#include <strings.h>

PetscErrorCode OpCreate__Poisson(Op op);

#ifdef __bgq__

#define NE 4

#define Pragma_(a) _Pragma(#a)
#define Pragma(a) Pragma_(a)

static inline PetscErrorCode OpPointwiseElement_PoissonN_QPX(Op op,PetscInt ne,PetscInt Q3,PetscScalar dx[3][3][Q3][NE],PetscReal wdxdet[Q3][NE],PetscScalar du[3][1][Q3][NE],PetscScalar dv[3][1][Q3][NE]) {
  Pragma(unroll(3))
  for (PetscInt i=0; i<Q3; i++) {
    vector4double dux[3],duX[3],dvx[3];
    vector4double dxi[3][3];
    for (PetscInt j=0; j<3; j++) {
      for (PetscInt k=0; k<3; k++) dxi[k][j] = vec_lda(0,dx[k][j][i]);
      duX[j] = vec_lda(0,du[j][0][i]);
    }
    for (PetscInt k=0; k<3; k++) {
      dux[k] = vec_madd(duX[2],dxi[k][2],vec_madd(duX[1],dxi[k][1],vec_mul(duX[0],dxi[k][0])));
    }
    for (PetscInt k=0; k<3; k++) dvx[k] = dux[k];
    for (PetscInt k=0; k<3; k++) {
      vector4double dvX = vec_madd(dvx[2],dxi[2][k],vec_madd(dvx[1],dxi[1][k],vec_mul(dvx[0],dxi[0][k])));
      vec_sta(vec_mul(dvX,vec_lda(0,wdxdet[i])),0,dv[k][0][i]);
    }
  }
  PetscLogFlops(Q3*ne*(5*3 + 0 + 6*3));
  return 0;
}
// Let the compiler specialize on known array sizes
static PetscErrorCode OpPointwiseElement_Poisson2_QPX(Op op,PetscInt ne,PetscInt Q3,PetscScalar dx[3][3][Q3][NE],PetscReal wdxdet[Q3][NE],PetscScalar du[3][1][Q3][NE],PetscScalar dv[3][1][Q3][NE]) {
  return OpPointwiseElement_PoissonN_QPX(op,NE,27,dx,wdxdet,du,dv);}


static PetscErrorCode OpApply_Poisson(Op op,DM dm,Vec U,Vec V,
                                      PetscErrorCode (*PointwiseElement)(Op op,PetscInt ne,PetscInt Q3,PetscScalar dx[3][3][Q3][NE],PetscReal wdxdet[Q3][NE],PetscScalar du[3][1][Q3][NE],PetscScalar dv[3][1][Q3][NE]))
{
  PetscErrorCode ierr;
  Vec X,Ul,Vl;
  PetscInt nelem,P,Q,P3,Q3;
  DM dmx;
  const PetscScalar *x,*u;
  PetscScalar *v;
  const PetscReal *B,*D,*w3;
  Tensor Tensor1,Tensor3;

  PetscFunctionBegin;
  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL,&w3);CHKERRQ(ierr);
  P3 = P*P*P;
  Q3 = Q*Q*Q;
  ierr = OpGetTensors(op,&Tensor1,&Tensor3);CHKERRQ(ierr);

  ierr = DMGetLocalVector(dm,&Ul);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&Vl);CHKERRQ(ierr);
  ierr = VecZeroEntries(Vl);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(dm,&dmx);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  ierr = DMFEGetNumElements(dm,&nelem);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,U,INSERT_VALUES,Ul);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,U,INSERT_VALUES,Ul);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Ul,&u);CHKERRQ(ierr);
  ierr = VecGetArray(Vl,&v);CHKERRQ(ierr);

  for (PetscInt e=0; e<nelem; e+=NE) {
    PetscScalar ve[1*P3*NE]_align,dv[3][1][Q3][NE]_align,ue[1*P3*NE]_align,du[3][1][Q3][NE]_align,xe[3*P3*NE]_align,dx[3][3][Q3][NE]_align,wdxdet[Q3][NE]_align;

    ierr = DMFEExtractElements(dmx,x,e,NE,xe);CHKERRQ(ierr);
    ierr = PetscMemzero(dx,sizeof dx);CHKERRQ(ierr);
    ierr = TensorContract(Tensor3,D,B,B,TENSOR_EVAL,xe,dx[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor3,B,D,B,TENSOR_EVAL,xe,dx[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor3,B,B,D,TENSOR_EVAL,xe,dx[2][0][0]);CHKERRQ(ierr);
    ierr = PointwiseJacobianInvert(NE,Q*Q*Q,w3,dx,wdxdet);CHKERRQ(ierr);
    ierr = DMFEExtractElements(dm,u,e,NE,ue);CHKERRQ(ierr);
    ierr = PetscMemzero(du,sizeof du);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,D,B,B,TENSOR_EVAL,ue,du[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,D,B,TENSOR_EVAL,ue,du[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,B,D,TENSOR_EVAL,ue,du[2][0][0]);CHKERRQ(ierr);
    ierr = PointwiseElement(op,NE,Q3,dx,wdxdet,du,dv);CHKERRQ(ierr);
    ierr = PetscMemzero(ve,sizeof ve);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,D,B,B,TENSOR_TRANSPOSE,dv[0][0][0],ve);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,D,B,TENSOR_TRANSPOSE,dv[1][0][0],ve);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,B,D,TENSOR_TRANSPOSE,dv[2][0][0],ve);CHKERRQ(ierr);
    ierr = DMFESetElements(dm,v,e,NE,ADD_VALUES,DOMAIN_INTERIOR,ve);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Ul,&u);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Ul);CHKERRQ(ierr);
  ierr = VecRestoreArray(Vl,&v);CHKERRQ(ierr);
  ierr = VecZeroEntries(V);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,Vl,ADD_VALUES,V);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,Vl,ADD_VALUES,V);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Vl);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
static PetscErrorCode OpApply_Poisson2(Op op,DM dm,Vec U,Vec V) { return OpApply_Poisson(op,dm,U,V,OpPointwiseElement_Poisson2_QPX); }

#endif

PetscErrorCode OpCreate_Poisson2_QPX(Op op)
{
#ifdef __bgq__
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = OpCreate__Poisson(op);CHKERRQ(ierr);
  ierr = OpSetFEDegree(op,2);CHKERRQ(ierr);
  ierr = OpSetApply(op,OpApply_Poisson2);CHKERRQ(ierr);
  ierr = OpSetPointwiseElement(op,(OpPointwiseElementFunction)OpPointwiseElement_Poisson2_QPX,NE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
#else
  PetscFunctionBegin;
  SETERRQ(OpComm(op),PETSC_ERR_SUP,"No support for poisson2-qpx on this machine");
  PetscFunctionReturn(0);
#endif
}
