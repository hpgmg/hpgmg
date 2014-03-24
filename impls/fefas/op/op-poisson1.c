#include "fefas-op.h"
#include "../tensor.h"
#include "../pointwise.h"

#define NE 4

typedef struct Ctx_private *Ctx;
struct Ctx_private {
  int dummy;
};

static PetscErrorCode OpPointwiseSolution_Poisson1(Op op,const PetscReal x[],PetscScalar u[]) {
  u[0] = PetscSinReal(1*PETSC_PI*x[0]) * PetscSinReal(2*PETSC_PI*x[1]) * PetscSinReal(3*PETSC_PI*x[2]);
  return 0;
}

static PetscErrorCode OpPointwiseForcing_Poisson1(Op op,const PetscReal x[],PetscScalar f[]) {
  f[0] = (  PetscSqr(1*PETSC_PI)*PetscSinReal(1*PETSC_PI*x[0]) * PetscSinReal(2*PETSC_PI*x[1]) * PetscSinReal(3*PETSC_PI*x[2])
          + PetscSqr(2*PETSC_PI)*PetscSinReal(1*PETSC_PI*x[0]) * PetscSinReal(2*PETSC_PI*x[1]) * PetscSinReal(3*PETSC_PI*x[2])
          + PetscSqr(3*PETSC_PI)*PetscSinReal(1*PETSC_PI*x[0]) * PetscSinReal(2*PETSC_PI*x[1]) * PetscSinReal(3*PETSC_PI*x[2]));
  return 0;
}

static PetscErrorCode PointwiseElement_Poisson1(PetscInt Q3,PetscScalar dx[3][3][Q3][NE],PetscReal wdxdet[Q3][NE],PetscScalar du[3][1][Q3][NE],PetscScalar dv[3][1][Q3][NE]) {
  for (PetscInt i=0; i<Q3; i++) {
    for (PetscInt e=0; e<NE; e++) {
      PetscScalar dux[3][1],dvx[3][1];
      for (PetscInt k=0; k<3; k++) {
        dux[k][0] = du[0][0][i][e] * dx[k][0][i][e] + du[1][0][i][e] * dx[k][1][i][e] + du[2][0][i][e] * dx[k][2][i][e];
      }
      for (PetscInt k=0; k<3; k++) dvx[k][0] = dux[k][0];
      for (PetscInt k=0; k<3; k++) {
        dv[k][0][i][e] = wdxdet[i][e] * (dvx[0][0] * dx[0][k][i][e] + dvx[1][0] * dx[1][k][i][e] + dvx[2][0] * dx[2][k][i][e]);
      }
    }
  }
  PetscLogFlops(Q3*NE*(5*3 + 0 + 6*3));
  return 0;
}

static PetscErrorCode OpApply_Poisson1(Op op,DM dm,Vec U,Vec V)
{
  PetscErrorCode ierr;
  Vec X,Ul,Vl;
  PetscInt nelem,P,Q,P3,Q3;
  DM dmx;
  const PetscScalar *x,*u;
  PetscScalar *v;
  const PetscReal *B,*D,*w3;

  PetscFunctionBegin;
  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL,&w3);CHKERRQ(ierr);
  P3 = P*P*P;
  Q3 = Q*Q*Q;

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
    PetscScalar ve[1*P3*NE],dv[3][1][Q3][NE],ue[1*P3*NE],du[3][1][Q3][NE],xe[3*P3*NE],dx[3][3][Q3][NE],wdxdet[Q3][NE];

    ierr = DMFEExtractElements(dmx,x,e,NE,xe);CHKERRQ(ierr);
    ierr = PetscMemzero(dx,sizeof dx);CHKERRQ(ierr);
    ierr = TensorContract(NE,3,P,Q,D,B,B,TENSOR_EVAL,xe,dx[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(NE,3,P,Q,B,D,B,TENSOR_EVAL,xe,dx[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(NE,3,P,Q,B,B,D,TENSOR_EVAL,xe,dx[2][0][0]);CHKERRQ(ierr);
    ierr = PointwiseJacobianInvert(NE,Q*Q*Q,w3,dx,wdxdet);CHKERRQ(ierr);
    ierr = DMFEExtractElements(dm,u,e,NE,ue);CHKERRQ(ierr);
    ierr = PetscMemzero(du,sizeof du);CHKERRQ(ierr);
    ierr = TensorContract(NE,1,P,Q,D,B,B,TENSOR_EVAL,ue,du[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(NE,1,P,Q,B,D,B,TENSOR_EVAL,ue,du[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(NE,1,P,Q,B,B,D,TENSOR_EVAL,ue,du[2][0][0]);CHKERRQ(ierr);
    ierr = PointwiseElement_Poisson1(Q3,dx,wdxdet,du,dv);CHKERRQ(ierr);
    ierr = PetscMemzero(ve,sizeof ve);CHKERRQ(ierr);
    ierr = TensorContract(NE,1,P,Q,D,B,B,TENSOR_TRANSPOSE,dv[0][0][0],ve);CHKERRQ(ierr);
    ierr = TensorContract(NE,1,P,Q,B,D,B,TENSOR_TRANSPOSE,dv[1][0][0],ve);CHKERRQ(ierr);
    ierr = TensorContract(NE,1,P,Q,B,B,D,TENSOR_TRANSPOSE,dv[2][0][0],ve);CHKERRQ(ierr);
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

static PetscErrorCode OpDestroy_Poisson1(Op op)
{
  PetscErrorCode ierr;
  Ctx ctx;

  PetscFunctionBegin;
  ierr = OpGetContext(op,&ctx);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpCreate_Poisson1(Op op)
{
  PetscErrorCode ierr;
  Ctx ctx;

  PetscFunctionBegin;
  ierr = PetscNew(&ctx);CHKERRQ(ierr);
  ierr = OpSetDof(op,1);CHKERRQ(ierr);
  ierr = OpSetFEDegree(op,1);CHKERRQ(ierr);
  ierr = OpSetContext(op,ctx);CHKERRQ(ierr);
  ierr = OpSetPointwiseSolution(op,OpPointwiseSolution_Poisson1);CHKERRQ(ierr);
  ierr = OpSetPointwiseForcing(op,OpPointwiseForcing_Poisson1);CHKERRQ(ierr);
  ierr = OpSetApply(op,OpApply_Poisson1);CHKERRQ(ierr);
  ierr = OpSetDestroy(op,OpDestroy_Poisson1);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
