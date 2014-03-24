#include "fefas-op.h"
#include "../tensor.h"
#include "../pointwise.h"

PetscErrorCode OpRegisterAll_Generated(void);

static PetscFunctionList OpList;
static PetscBool OpPackageInitialized;

struct Op_private {
  MPI_Comm comm;                /* Finest level comm (only for diagnostics at setup time) */
  PetscInt fedegree;
  PetscInt dof;
  PetscErrorCode (*Apply)(Op,DM,Vec,Vec);
  PetscErrorCode (*RestrictState)(Op,DM,Vec,Vec);
  PetscErrorCode (*RestrictResidual)(Op,DM,Vec,Vec);
  PetscErrorCode (*Interpolate)(Op,DM,Vec,Vec);
  PetscErrorCode (*PointwiseSolution)(Op,const PetscReal[],PetscScalar[]);
  PetscErrorCode (*PointwiseForcing)(Op,const PetscReal[],PetscScalar[]);
  PetscErrorCode (*Destroy)(Op);
  void *ctx;
};

PetscErrorCode OpSetFEDegree(Op op,PetscInt degree) {
  op->fedegree = degree;
  return 0;
}
PetscErrorCode OpGetFEDegree(Op op,PetscInt *degree) {
  *degree = op->fedegree;
  return 0;
}
PetscErrorCode OpSetDof(Op op,PetscInt dof) {
  op->dof = dof;
  return 0;
}
PetscErrorCode OpGetDof(Op op,PetscInt *dof) {
  *dof = op->dof;
  return 0;
}
PetscErrorCode OpSetContext(Op op,void *ctx) {
  op->ctx = ctx;
  return 0;
}
PetscErrorCode OpGetContext(Op op,void *ctx) {
  *(void**)ctx = op->ctx;
  return 0;
}
PetscErrorCode OpSetDestroy(Op op,PetscErrorCode (*f)(Op)) {
  op->Destroy = f;
  return 0;
}
PetscErrorCode OpSetApply(Op op,PetscErrorCode (*f)(Op,DM,Vec,Vec)) {
  op->Apply = f;
  return 0;
}
PetscErrorCode OpSetPointwiseSolution(Op op,PetscErrorCode (*f)(Op,const PetscReal[],PetscScalar[])) {
  op->PointwiseSolution = f;
  return 0;
}
PetscErrorCode OpSetPointwiseForcing(Op op,PetscErrorCode (*f)(Op,const PetscReal[],PetscScalar[])) {
  op->PointwiseForcing = f;
  return 0;
}

PetscErrorCode OpSolution(Op op,DM dm,Vec U) {
  PetscErrorCode ierr;
  Vec X;
  const PetscScalar *x;
  PetscScalar *u;
  PetscInt i,m,bs;

  PetscFunctionBegin;
  ierr = DMGetCoordinates(dm,&X);CHKERRQ(ierr);
  ierr = VecGetLocalSize(U,&m);CHKERRQ(ierr);
  ierr = VecGetBlockSize(U,&bs);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  for (i=0; i<m/bs; i++) {
    ierr = (op->PointwiseSolution)(op,&x[i*3],&u[i*bs]);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpForcing(Op op,DM dm,Vec F) {
  PetscErrorCode ierr;
  Vec X,Floc;
  DM dmx;
  const PetscScalar *x;
  PetscScalar *f;
  const PetscReal *B,*D,*w3;
  PetscInt nelem,ne = 1,P,Q,P3,Q3;

  PetscFunctionBegin;
  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL,&w3);CHKERRQ(ierr);
  P3 = P*P*P;
  Q3 = Q*Q*Q;

  ierr = DMGetLocalVector(dm,&Floc);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(dm,&dmx);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  ierr = DMFEGetNumElements(dm,&nelem);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecZeroEntries(Floc);CHKERRQ(ierr);
  ierr = VecGetArray(Floc,&f);CHKERRQ(ierr);

  for (PetscInt e=0; e<nelem; e+=ne) {
    PetscScalar fe[op->dof*P3*ne],fq[op->dof][Q3][ne],xe[3*P3*ne],xq[3][Q3][ne],dx[3][3][Q3][ne],wdxdet[Q3][ne];

    ierr = DMFEExtractElements(dmx,x,e,ne,xe);CHKERRQ(ierr);
    ierr = PetscMemzero(xq,sizeof xq);CHKERRQ(ierr);
    ierr = TensorContract(ne,3,P,Q,B,B,B,TENSOR_EVAL,xe,xq[0][0]);CHKERRQ(ierr);
    ierr = PetscMemzero(dx,sizeof dx);CHKERRQ(ierr);
    ierr = TensorContract(ne,3,P,Q,D,B,B,TENSOR_EVAL,xe,dx[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(ne,3,P,Q,B,D,B,TENSOR_EVAL,xe,dx[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(ne,3,P,Q,B,B,D,TENSOR_EVAL,xe,dx[2][0][0]);CHKERRQ(ierr);
    ierr = PointwiseJacobianInvert(ne,Q*Q*Q,w3,dx,wdxdet);CHKERRQ(ierr);

    for (PetscInt i=0; i<Q3; i++) {
      for (PetscInt l=0; l<ne; l++) {
        PetscReal xx[] = {xq[0][i][l],xq[1][i][l],xq[2][i][l]};
        PetscScalar fql[op->dof];
        ierr = (op->PointwiseForcing)(op,xx,fql);CHKERRQ(ierr);
        for (PetscInt d=0; d<op->dof; d++) fq[d][i][l] = wdxdet[i][l] * fql[d];
      }
    }
    ierr = PetscMemzero(fe,sizeof fe);CHKERRQ(ierr);
    ierr = TensorContract(ne,op->dof,P,Q,B,B,B,TENSOR_TRANSPOSE,fq[0][0],fe);CHKERRQ(ierr);
    ierr = DMFESetElements(dm,f,e,ne,ADD_VALUES,DOMAIN_INTERIOR,fe);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Floc,&f);CHKERRQ(ierr);
  ierr = VecZeroEntries(F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,Floc,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,Floc,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Floc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpApply(Op op,DM dm,Vec U,Vec F) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!op->Apply) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"No Apply implemented, use OpSetApply()");
  ierr = (*op->Apply)(op,dm,U,F);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpRestrictState(Op op,DM dm,Vec Uf,Vec Uc) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (op->RestrictState) {
    ierr = (*op->RestrictState)(op,dm,Uf,Uc);CHKERRQ(ierr);
  } else {
    ierr = DMFEInject(dm,Uf,Uc);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode OpRestrictResidual(Op op,DM dm,Vec Uf,Vec Uc) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (op->RestrictResidual) {
    ierr = (*op->RestrictResidual)(op,dm,Uf,Uc);CHKERRQ(ierr);
  } else {
    ierr = DMFERestrict(dm,Uf,Uc);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode OpInterpolate(Op op,DM dm,Vec Uc,Vec Uf) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (op->Interpolate) {
    ierr = (*op->Interpolate)(op,dm,Uc,Uf);CHKERRQ(ierr);
  } else {
    ierr = DMFEInterpolate(dm,Uc,Uf);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode OpCreateFromOptions(MPI_Comm comm,Op *op)
{
  PetscErrorCode ierr,(*f)(Op);
  Op o;
  char opname[256] = "poisson1";

  PetscFunctionBegin;
  ierr = OpInitializePackage();CHKERRQ(ierr);
  ierr = PetscNew(&o);CHKERRQ(ierr);
  ierr = PetscCommDuplicate(comm,&o->comm,NULL);CHKERRQ(ierr);
  ierr = OpSetFEDegree(o,1);CHKERRQ(ierr);
  ierr = OpSetDof(o,1);CHKERRQ(ierr);
  ierr = PetscOptionsBegin(comm,NULL,"Operator options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsFList("-op_type","Operator type","",OpList,opname,opname,sizeof opname,NULL);CHKERRQ(ierr);
  ierr = PetscFunctionListFind(OpList,opname,&f);CHKERRQ(ierr);
  ierr = (*f)(o);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  *op = o;
  PetscFunctionReturn(0);
}

PetscErrorCode OpDestroy(Op *op)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*op) PetscFunctionReturn(0);
  if ((*op)->Destroy) {
    ierr = (*op)->Destroy(*op);CHKERRQ(ierr);
  }
  ierr = PetscCommDestroy(&(*op)->comm);CHKERRQ(ierr);
  ierr = PetscFree(*op);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpRegister(const char *name,PetscErrorCode (*f)(Op))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListAdd(&OpList,name,f);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpFinalizePackage()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListDestroy(&OpList);CHKERRQ(ierr);
  OpPackageInitialized = PETSC_FALSE;
  PetscFunctionReturn(0);
}

PetscErrorCode OpInitializePackage()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (OpPackageInitialized) PetscFunctionReturn(0);
  ierr = OpRegisterAll_Generated();CHKERRQ(ierr);
  ierr = PetscRegisterFinalize(OpFinalizePackage);CHKERRQ(ierr);
  OpPackageInitialized = PETSC_TRUE;
  PetscFunctionReturn(0);
}
