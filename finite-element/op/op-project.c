#include "fefas-op.h"
#include "../fefas.h"
#include "../tensor.h"
#include <strings.h>

#define NE 8

typedef struct Ctx_private *Ctx;
struct Ctx_private {
  int dummy;
};

static PetscErrorCode OpPointwiseSolution_Project_Sine(Op op,const PetscReal x[],const PetscReal L[],PetscScalar u[]) {
  u[0] = PetscSinReal(1*PETSC_PI*x[0]/L[0]) * PetscSinReal(2*PETSC_PI*x[1]/L[1]) * PetscSinReal(3*PETSC_PI*x[2]/L[2]);
  return 0;
}

static inline PetscErrorCode OpPointwiseElement_ProjectN(Op op,PetscInt ne,PetscInt Q3,PetscScalar dx[3][3][Q3][NE],PetscReal wdxdet[Q3][NE],PetscScalar uu[1][1][Q3][NE],PetscScalar du[3][1][Q3][NE],PetscScalar vv[1][1][Q3][NE],PetscScalar dv[3][1][Q3][NE]) {
  for (PetscInt i=0; i<Q3; i++) {
    for (PetscInt e=0; e<ne; e++) {
      vv[0][0][i][e] = wdxdet[i][e] * uu[0][0][i][e];
    }
  }
  PetscLogFlops(Q3*ne);
  return 0;
}

static PetscErrorCode OpApply_Project2Affine(Op op,DM dm,Vec U,Vec V)
{
  PetscErrorCode ierr;
  Vec Ul,Vl;
  PetscInt nelem,P,Q,P3,Q3,Mglobal[3];
  PetscReal L[3];
  const PetscScalar *u;
  PetscScalar *v;
  const PetscReal *B,*D,*w3;
  Tensor Tensor1,Tensor3;

  PetscFunctionBegin;
  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL,&w3);CHKERRQ(ierr);
  P3 = P*P*P;
  Q3 = Q*Q*Q;
  ierr = OpGetTensors(op,&Tensor1,&Tensor3);CHKERRQ(ierr);
  ierr = DMFEGetInfo(dm,NULL,NULL,NULL,Mglobal,NULL);CHKERRQ(ierr);
  ierr = DMFEGetUniformCoordinates(dm,L);CHKERRQ(ierr);

  ierr = DMGetLocalVector(dm,&Ul);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&Vl);CHKERRQ(ierr);
  ierr = VecZeroEntries(Vl);CHKERRQ(ierr);
  ierr = DMFEGetNumElements(dm,&nelem);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,U,INSERT_VALUES,Ul);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,U,INSERT_VALUES,Ul);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Ul,&u);CHKERRQ(ierr);
  ierr = VecGetArray(Vl,&v);CHKERRQ(ierr);

  for (PetscInt e=0; e<nelem; e+=NE) {
    PetscScalar ve[1*P3*NE]_align,vv[3][1][Q3][NE]_align,ue[1*P3*NE]_align,uu[1][1][Q3][NE]_align,dx[3],wdxdet[Q3];

    for (PetscInt i=0; i<3; i++) dx[i] = 2.*Mglobal[i]/L[i];
    for (PetscInt i=0; i<Q3; i++) wdxdet[i] = w3[i] / (dx[0]*dx[1]*dx[2]);
    ierr = DMFEExtractElements(dm,u,e,NE,ue);CHKERRQ(ierr);
    ierr = PetscMemzero(uu,sizeof uu);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,B,B,TENSOR_EVAL,ue,uu[0][0][0]);CHKERRQ(ierr);
    for (PetscInt i=0; i<Q3; i++) {
      for (PetscInt j=0; j<3; j++) {
        for (PetscInt l=0; l<NE; l++) {
          vv[j][0][i][l] = wdxdet[i] * uu[j][0][i][l];
        }
      }
    }
    ierr = PetscMemzero(ve,sizeof ve);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,B,B,TENSOR_TRANSPOSE,vv[0][0][0],ve);CHKERRQ(ierr);
    ierr = DMFESetElements(dm,v,e,NE,ADD_VALUES,DOMAIN_INTERIOR,ve);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(Ul,&u);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Ul);CHKERRQ(ierr);
  ierr = VecRestoreArray(Vl,&v);CHKERRQ(ierr);
  ierr = VecZeroEntries(V);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,Vl,ADD_VALUES,V);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,Vl,ADD_VALUES,V);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Vl);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode OpDestroy_Project(Op op)
{
  PetscErrorCode ierr;
  Ctx ctx;

  PetscFunctionBegin;
  ierr = OpGetContext(op,&ctx);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpCreate_Project2Affine(Op op)
{
  PetscErrorCode ierr;
  Ctx ctx;

  PetscFunctionBegin;
  ierr = OpSetPointwiseSolution(op,OpPointwiseSolution_Project_Sine);CHKERRQ(ierr);
  ierr = OpSetPointwiseForcing(op,OpPointwiseSolution_Project_Sine);CHKERRQ(ierr);
  ierr = OpSetDof(op,1);CHKERRQ(ierr);
  ierr = PetscNew(&ctx);CHKERRQ(ierr);
  ierr = OpSetContext(op,ctx);CHKERRQ(ierr);
  ierr = OpSetDestroy(op,OpDestroy_Project);CHKERRQ(ierr);
  ierr = OpSetFEDegree(op,2);CHKERRQ(ierr);
  ierr = OpSetApply(op,OpApply_Project2Affine);CHKERRQ(ierr);
  ierr = OpSetAffineOnly(op,PETSC_TRUE);CHKERRQ(ierr);
  ierr = OpSetPointwiseElement(op,(OpPointwiseElementFunction)OpPointwiseElement_ProjectN,NE,011);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
