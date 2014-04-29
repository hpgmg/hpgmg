#include "fefas-op.h"
#include "../fefas.h"
#include "../tensor.h"
#include "../pointwise.h"
#include <strings.h>

#define NE 4

typedef struct Ctx_private *Ctx;
struct Ctx_private {
  int dummy;
};

static PetscReal Sqr(PetscReal x) {return x*x;}

static PetscErrorCode OpPointwiseSolution_Poisson_Sine(Op op,const PetscReal x[],const PetscReal L[],PetscScalar u[]) {
  u[0] = PetscSinReal(1*PETSC_PI*x[0]/L[0]) * PetscSinReal(2*PETSC_PI*x[1]/L[1]) * PetscSinReal(3*PETSC_PI*x[2]/L[2]);
  return 0;
}
static PetscErrorCode OpPointwiseForcing_Poisson_Sine(Op op,const PetscReal x[],const PetscReal L[],PetscScalar f[]) {
  f[0] = (  Sqr(1*PETSC_PI/L[0])*PetscSinReal(1*PETSC_PI*x[0]/L[0]) * PetscSinReal(2*PETSC_PI*x[1]/L[1]) * PetscSinReal(3*PETSC_PI*x[2]/L[2])
          + Sqr(2*PETSC_PI/L[1])*PetscSinReal(1*PETSC_PI*x[0]/L[0]) * PetscSinReal(2*PETSC_PI*x[1]/L[1]) * PetscSinReal(3*PETSC_PI*x[2]/L[2])
          + Sqr(3*PETSC_PI/L[2])*PetscSinReal(1*PETSC_PI*x[0]/L[0]) * PetscSinReal(2*PETSC_PI*x[1]/L[1]) * PetscSinReal(3*PETSC_PI*x[2]/L[2]));
  return 0;
}

static PetscReal Hump(const PetscReal x[],const PetscReal L[]){
  return PetscSinReal(PETSC_PI*x[0]/L[0]) * PetscSinReal(PETSC_PI*x[1]/L[1]) * PetscSinReal(PETSC_PI*x[2]/L[2]);}
static PetscReal Hump_x0(const PetscReal x[],const PetscReal L[]) {
  return PETSC_PI/L[0]*PetscCosReal(PETSC_PI*x[0]/L[0]) * PetscSinReal(PETSC_PI*x[1]/L[1]) * PetscSinReal(PETSC_PI*x[2]/L[2]);}
static PetscReal Hump_x1(const PetscReal x[],const PetscReal L[]) {
  return PETSC_PI/L[1]*PetscSinReal(PETSC_PI*x[0]/L[0]) * PetscCosReal(PETSC_PI*x[1]/L[1]) * PetscSinReal(PETSC_PI*x[2]/L[2]);}
static PetscReal Hump_x2(const PetscReal x[],const PetscReal L[]) {
  return PETSC_PI/L[2]*PetscSinReal(PETSC_PI*x[0]/L[0]) * PetscSinReal(PETSC_PI*x[1]/L[1]) * PetscCosReal(PETSC_PI*x[2]/L[2]);}
static PetscReal Hump_xx0(const PetscReal x[],const PetscReal L[]) {
  return -Sqr(PETSC_PI/L[0])*Hump(x,L);}
static PetscReal Hump_xx1(const PetscReal x[],const PetscReal L[]) {
  return -Sqr(PETSC_PI/L[1])*Hump(x,L);}
static PetscReal Hump_xx2(const PetscReal x[],const PetscReal L[]) {
  return -Sqr(PETSC_PI/L[2])*Hump(x,L);}

static PetscReal Bend(const PetscReal x[],const PetscReal L[]) {
  return PetscTanhReal(x[0]/L[0]) + PetscLogReal(1 + x[1]/L[1]) + PetscExpReal(-x[2]/L[2]);}
static PetscReal Bend_x0(const PetscReal x[],const PetscReal L[]) {
  return (-Sqr(PetscTanhReal(x[0]/L[0])) + 1)/L[0];}
static PetscReal Bend_xx0(const PetscReal x[],const PetscReal L[]) {
  return -(-2*Sqr(PetscTanhReal(x[0]/L[0])) + 2)*PetscTanhReal(x[0]/L[0])/Sqr(L[0]);}
static PetscReal Bend_x1(const PetscReal x[],const PetscReal L[]) {
  return 1/(L[1]*(1 + x[1]/L[1]));}
static PetscReal Bend_xx1(const PetscReal x[],const PetscReal L[]) {
  return -1/(Sqr(L[1])*Sqr(1 + x[1]/L[1]));}
static PetscReal Bend_x2(const PetscReal x[],const PetscReal L[]) {
  return -PetscExpReal(-x[2]/L[2])/L[2];}
static PetscReal Bend_xx2(const PetscReal x[],const PetscReal L[]) {
  return PetscExpReal(-x[2]/L[2])/Sqr(L[2]);}
static PetscErrorCode OpPointwiseSolution_Poisson_Hump(Op op,const PetscReal x[],const PetscReal L[],PetscScalar u[]) {
  u[0] = Hump(x,L) * Bend(x,L);
  return 0;
}
// (f g)'' = f'' g + 2 f' g' + f g''
static PetscErrorCode OpPointwiseForcing_Poisson_Hump(Op op,const PetscReal x[],const PetscReal L[],PetscScalar f[]) {
  f[0] = -(  Hump_xx0(x,L)*Bend(x,L) + 2*Hump_x0(x,L)*Bend_x0(x,L) + Hump(x,L)*Bend_xx0(x,L)
           + Hump_xx1(x,L)*Bend(x,L) + 2*Hump_x1(x,L)*Bend_x1(x,L) + Hump(x,L)*Bend_xx1(x,L)
           + Hump_xx2(x,L)*Bend(x,L) + 2*Hump_x2(x,L)*Bend_x2(x,L) + Hump(x,L)*Bend_xx2(x,L));
  return 0;
}

static PetscReal Wave(const PetscReal x) {
  return x*x*x*x - x*x + 2*x*x*x - 2*x*x*x*x*x;}
static PETSC_UNUSED PetscReal Wave_x(const PetscReal x) {
  return 4*x*x*x - 2*x + 6*x*x - 10*x*x*x*x;}
static PetscReal Wave_xx(const PetscReal x) {
  return 12*x*x - 2 + 12*x - 40*x*x*x;}
static PetscErrorCode OpPointwiseSolution_Poisson_Wave(Op op,const PetscReal x[],const PetscReal L[],PetscScalar u[]) {
  u[0] = Wave(x[0]/L[0]) * Wave(x[1]/L[1]) * Wave(x[2]/L[2]);
  return 0;
}
static PetscErrorCode OpPointwiseForcing_Poisson_Wave(Op op,const PetscReal x[],const PetscReal L[],PetscScalar f[]) {
  f[0] = -(  Wave_xx(x[0]/L[0]) * Wave   (x[1]/L[1]) * Wave   (x[2]/L[2])/Sqr(L[0])
           + Wave   (x[0]/L[0]) * Wave_xx(x[1]/L[1]) * Wave   (x[2]/L[2])/Sqr(L[1])
           + Wave   (x[0]/L[0]) * Wave   (x[1]/L[1]) * Wave_xx(x[2]/L[2])/Sqr(L[2]));
  return 0;
}

static inline PetscErrorCode OpPointwiseElement_PoissonN(Op op,PetscInt ne,PetscInt Q3,PetscScalar dx[3][3][Q3][NE],PetscReal wdxdet[Q3][NE],PetscScalar du[3][1][Q3][NE],PetscScalar dv[3][1][Q3][NE]) {
  for (PetscInt i=0; i<Q3; i++) {
    for (PetscInt e=0; e<ne; e++) {
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
  PetscLogFlops(Q3*ne*(5*3 + 0 + 6*3));
  return 0;
}
// Let the compiler specialize on known array sizes
static PetscErrorCode OpPointwiseElement_Poisson1(Op op,PetscInt ne,PetscInt Q3,PetscScalar dx[3][3][Q3][NE],PetscReal wdxdet[Q3][NE],PetscScalar du[3][1][Q3][NE],PetscScalar dv[3][1][Q3][NE]) {
  return OpPointwiseElement_PoissonN(op,NE,8,dx,wdxdet,du,dv);}
static PetscErrorCode OpPointwiseElement_Poisson2(Op op,PetscInt ne,PetscInt Q3,PetscScalar dx[3][3][Q3][NE],PetscReal wdxdet[Q3][NE],PetscScalar du[3][1][Q3][NE],PetscScalar dv[3][1][Q3][NE]) {
  return OpPointwiseElement_PoissonN(op,NE,27,dx,wdxdet,du,dv);}


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
// Not trying to inline yet because Q1 is faster without tensor product, so it's really worth specializing the whole
// element loop.  The non-tensor code is slightly simpler, but not implemented yet, so delay.  I don't think we'll end
// up using Q1 anyway at the end of the day.
static PetscErrorCode OpApply_Poisson1(Op op,DM dm,Vec U,Vec V) { return OpApply_Poisson(op,dm,U,V,OpPointwiseElement_Poisson1); }
static PetscErrorCode OpApply_Poisson2(Op op,DM dm,Vec U,Vec V) { return OpApply_Poisson(op,dm,U,V,OpPointwiseElement_Poisson2); }

static PetscErrorCode OpApply_Poisson2Affine(Op op,DM dm,Vec U,Vec V)
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
    PetscScalar ve[1*P3*NE]_align,dv[3][1][Q3][NE]_align,ue[1*P3*NE]_align,du[3][1][Q3][NE]_align,dx[3],wdxdet[Q3];

    for (PetscInt i=0; i<3; i++) dx[i] = 2.*Mglobal[i]/L[i];
    for (PetscInt i=0; i<Q3; i++) wdxdet[i] = w3[i] / (dx[0]*dx[1]*dx[2]);
    ierr = DMFEExtractElements(dm,u,e,NE,ue);CHKERRQ(ierr);
    ierr = PetscMemzero(du,sizeof du);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,D,B,B,TENSOR_EVAL,ue,du[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,D,B,TENSOR_EVAL,ue,du[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,B,D,TENSOR_EVAL,ue,du[2][0][0]);CHKERRQ(ierr);
    for (PetscInt i=0; i<Q3; i++) {
      for (PetscInt j=0; j<3; j++) {
        for (PetscInt l=0; l<NE; l++) {
          dv[j][0][i][l] = wdxdet[i] * dx[j] * dx[j] * du[j][0][i][l];
        }
      }
    }
    ierr = PetscMemzero(ve,sizeof ve);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,D,B,B,TENSOR_TRANSPOSE,dv[0][0][0],ve);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,D,B,TENSOR_TRANSPOSE,dv[1][0][0],ve);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,B,D,TENSOR_TRANSPOSE,dv[2][0][0],ve);CHKERRQ(ierr);
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

static PetscErrorCode OpDestroy_Poisson(Op op)
{
  PetscErrorCode ierr;
  Ctx ctx;

  PetscFunctionBegin;
  ierr = OpGetContext(op,&ctx);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpCreate__Poisson(Op op)
{
  PetscErrorCode ierr;
  Ctx ctx;
  char solname[256] = "wave";

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(OpComm(op),NULL,"Poisson Options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-poisson_solution","Name of manufactured/analytic solution","",solname,solname,sizeof solname,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  if (!strcasecmp(solname,"sine")) {
    ierr = OpSetPointwiseSolution(op,OpPointwiseSolution_Poisson_Sine);CHKERRQ(ierr);
    ierr = OpSetPointwiseForcing(op,OpPointwiseForcing_Poisson_Sine);CHKERRQ(ierr);
  } else if (!strcasecmp(solname,"hump")) {
    ierr = OpSetPointwiseSolution(op,OpPointwiseSolution_Poisson_Hump);CHKERRQ(ierr);
    ierr = OpSetPointwiseForcing(op,OpPointwiseForcing_Poisson_Hump);CHKERRQ(ierr);
  } else if (!strcasecmp(solname,"wave")) {
    ierr = OpSetPointwiseSolution(op,OpPointwiseSolution_Poisson_Wave);CHKERRQ(ierr);
    ierr = OpSetPointwiseForcing(op,OpPointwiseForcing_Poisson_Wave);CHKERRQ(ierr);
  }
  ierr = PetscNew(&ctx);CHKERRQ(ierr);
  ierr = OpSetDof(op,1);CHKERRQ(ierr);
  ierr = OpSetContext(op,ctx);CHKERRQ(ierr);
  ierr = OpSetDestroy(op,OpDestroy_Poisson);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpCreate_Poisson1(Op op)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = OpCreate__Poisson(op);CHKERRQ(ierr);
  ierr = OpSetFEDegree(op,1);CHKERRQ(ierr);
  ierr = OpSetApply(op,OpApply_Poisson1);CHKERRQ(ierr);
  ierr = OpSetPointwiseElement(op,(OpPointwiseElementFunction)OpPointwiseElement_Poisson1,NE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode OpCreate_Poisson2(Op op)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = OpCreate__Poisson(op);CHKERRQ(ierr);
  ierr = OpSetFEDegree(op,2);CHKERRQ(ierr);
  ierr = OpSetApply(op,OpApply_Poisson2);CHKERRQ(ierr);
  ierr = OpSetPointwiseElement(op,(OpPointwiseElementFunction)OpPointwiseElement_Poisson2,NE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
PetscErrorCode OpCreate_Poisson2Affine(Op op)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = OpCreate__Poisson(op);CHKERRQ(ierr);
  ierr = OpSetFEDegree(op,2);CHKERRQ(ierr);
  // This implementation requires both affine and non-rotated.
  ierr = OpSetApply(op,OpApply_Poisson2Affine);CHKERRQ(ierr);
  ierr = OpSetAffineOnly(op,PETSC_TRUE);CHKERRQ(ierr);
  ierr = OpSetPointwiseElement(op,(OpPointwiseElementFunction)OpPointwiseElement_Poisson2,NE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
