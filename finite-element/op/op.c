#include "fefas-op.h"
#include "../fefas.h"
#include "../tensor.h"
#include "../pointwise.h"

PetscErrorCode OpRegisterAll_Generated(void);

static PetscFunctionList OpList;
static PetscBool OpPackageInitialized;
static PetscLogEvent OP_Apply,OP_RestrictState,OP_RestrictResidual,OP_Interpolate,OP_Solution,OP_Forcing,OP_IntegrateNorms,OP_GetDiagonal;

struct Op_private {
  MPI_Comm comm;                /* Finest level comm (only for diagnostics at setup time) */
  PetscInt fedegree;
  PetscInt dof;
  PetscInt addquadpts;
  PetscInt ne;                  /* Preferred number of elements over which to vectorize */
  PetscInt opcode;              /* 01=uu | 02=du, 010=vv | 020=dv */
  PetscBool affineonly;
  Tensor   TensorDOF,Tensor3;
  PetscErrorCode (*Apply)(Op,DM,Vec,Vec);
  PetscErrorCode (*RestrictState)(Op,DM,Vec,Vec);
  PetscErrorCode (*RestrictResidual)(Op,DM,Vec,Vec);
  PetscErrorCode (*Interpolate)(Op,DM,Vec,Vec);
  PetscErrorCode (*PointwiseSolution)(Op,const PetscReal[],const PetscReal[],PetscScalar[]);
  PetscErrorCode (*PointwiseForcing)(Op,const PetscReal[],const PetscReal[],PetscScalar[]);
  PetscErrorCode (*PointwiseElement)(Op,PetscInt,PetscInt,const PetscScalar[],const PetscReal[],const PetscScalar[],const PetscScalar[],PetscScalar[],PetscScalar[]);
  PetscErrorCode (*Destroy)(Op);
  void *ctx;
};

MPI_Comm OpComm(Op op) { return op->comm; }
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
PetscErrorCode OpGetAddQuadPts(Op op,PetscInt *addquadpts) {
  *addquadpts = op->addquadpts;
  return 0;
}
PetscErrorCode OpSetAddQuadPts(Op op,PetscInt addquadpts) {
  op->addquadpts = addquadpts;
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
PetscErrorCode OpSetPointwiseSolution(Op op,PetscErrorCode (*f)(Op,const PetscReal[],const PetscReal[],PetscScalar[])) {
  op->PointwiseSolution = f;
  return 0;
}
PetscErrorCode OpSetPointwiseForcing(Op op,PetscErrorCode (*f)(Op,const PetscReal[],const PetscReal[],PetscScalar[])) {
  op->PointwiseForcing = f;
  return 0;
}
PetscErrorCode OpSetPointwiseElement(Op op,OpPointwiseElementFunction f,PetscInt ne,PetscInt opcode) {
  op->PointwiseElement = f;
  op->ne = ne;
  op->opcode = opcode;
  return 0;
}
PetscErrorCode OpSetAffineOnly(Op op,PetscBool affine) {
  op->affineonly = affine;
  return 0;
}
PetscErrorCode OpGetAffineOnly(Op op,PetscBool *affine) {
  *affine = op->affineonly;
  return 0;
}

PetscErrorCode OpSolution(Op op,DM dm,Vec U) {
  PetscErrorCode ierr;
  Vec X;
  const PetscScalar *x;
  PetscScalar *u;
  PetscReal L[3];
  PetscInt i,m,bs;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(OP_Solution,dm,U,0,0);CHKERRQ(ierr);
  ierr = DMGetCoordinates(dm,&X);CHKERRQ(ierr);
  ierr = VecStrideMax(X,0,NULL,&L[0]);CHKERRQ(ierr);
  ierr = VecStrideMax(X,1,NULL,&L[1]);CHKERRQ(ierr);
  ierr = VecStrideMax(X,2,NULL,&L[2]);CHKERRQ(ierr);
  ierr = VecGetLocalSize(U,&m);CHKERRQ(ierr);
  ierr = VecGetBlockSize(U,&bs);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  for (i=0; i<m/bs; i++) {
    ierr = (op->PointwiseSolution)(op,&x[i*3],L,&u[i*bs]);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(OP_Solution,dm,U,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpForcing(Op op,DM dm,Vec F) {
  PetscErrorCode ierr;
  Vec X,Floc;
  DM dmx;
  const PetscScalar *x;
  PetscScalar *f;
  const PetscReal *B,*D,*w3;
  PetscReal L[3];
  PetscInt nelem,ne = op->ne,P,Q,P3,Q3;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(OP_Forcing,dm,F,0,0);CHKERRQ(ierr);
  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL,&w3);CHKERRQ(ierr);
  P3 = P*P*P;
  Q3 = Q*Q*Q;

  ierr = DMFEGetUniformCoordinates(dm,L);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&Floc);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(dm,&dmx);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  ierr = DMFEGetNumElements(dm,&nelem);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecZeroEntries(Floc);CHKERRQ(ierr);
  ierr = VecGetArray(Floc,&f);CHKERRQ(ierr);

  for (PetscInt e=0; e<nelem; e+=ne) {
    PetscScalar fe[op->dof*P3*ne]_align,fq[op->dof][Q3][ne]_align,xe[3*P3*ne]_align,xq[3][Q3][ne]_align,dx[3][3][Q3][ne]_align,wdxdet[Q3][ne]_align;

    ierr = DMFEExtractElements(dmx,x,e,ne,xe);CHKERRQ(ierr);
    ierr = PetscMemzero(xq,sizeof xq);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,B,B,B,TENSOR_EVAL,xe,xq[0][0]);CHKERRQ(ierr);
    ierr = PetscMemzero(dx,sizeof dx);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,D,B,B,TENSOR_EVAL,xe,dx[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,B,D,B,TENSOR_EVAL,xe,dx[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,B,B,D,TENSOR_EVAL,xe,dx[2][0][0]);CHKERRQ(ierr);
    ierr = PointwiseJacobianInvert(ne,Q*Q*Q,w3,dx,wdxdet);CHKERRQ(ierr);

    for (PetscInt i=0; i<Q3; i++) {
      for (PetscInt l=0; l<ne; l++) {
        PetscReal xx[] = {xq[0][i][l],xq[1][i][l],xq[2][i][l]};
        PetscScalar fql[op->dof];
        ierr = (op->PointwiseForcing)(op,xx,L,fql);CHKERRQ(ierr);
        for (PetscInt d=0; d<op->dof; d++) fq[d][i][l] = wdxdet[i][l] * fql[d];
      }
    }
    ierr = PetscMemzero(fe,sizeof fe);CHKERRQ(ierr);
    ierr = TensorContract(op->TensorDOF,B,B,B,TENSOR_TRANSPOSE,fq[0][0],fe);CHKERRQ(ierr);
    ierr = DMFESetElements(dm,f,e,ne,ADD_VALUES,DOMAIN_INTERIOR,fe);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Floc,&f);CHKERRQ(ierr);
  ierr = VecZeroEntries(F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,Floc,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,Floc,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Floc);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(OP_Forcing,dm,F,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Approximates continuous L^\infty and L^2 norms of error, normalized by L^\infty and L^2 norms of analytic solution
PetscErrorCode OpIntegrateNorms(Op op,DM dm,Vec U,PetscReal *normInfty,PetscReal *normL2) {
  PetscErrorCode ierr;
  Vec X,Uloc;
  DM dmx;
  const PetscScalar *x,*u;
  const PetscReal *B,*D,*w3;
  PetscReal L[3];
  struct {PetscReal error,u;} sumInfty={},sum2={};
  PetscInt nelem,ne = op->ne,P,Q,P3,Q3;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(OP_IntegrateNorms,dm,U,0,0);CHKERRQ(ierr);
  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL,&w3);CHKERRQ(ierr);
  P3 = P*P*P;
  Q3 = Q*Q*Q;

  ierr = DMFEGetUniformCoordinates(dm,L);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(dm,&dmx);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm,&Uloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,U,INSERT_VALUES,Uloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,U,INSERT_VALUES,Uloc);CHKERRQ(ierr);
  ierr = DMFEGetNumElements(dm,&nelem);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Uloc,&u);CHKERRQ(ierr);

  for (PetscInt e=0; e<nelem; e+=ne) {
    PetscScalar ue[op->dof*P3*ne]_align,uq[op->dof][Q3][ne]_align,xe[3*P3*ne]_align,xq[3][Q3][ne]_align,dx[3][3][Q3][ne]_align,wdxdet[Q3][ne]_align;

    ierr = DMFEExtractElements(dmx,x,e,ne,xe);CHKERRQ(ierr);
    ierr = PetscMemzero(xq,sizeof xq);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,B,B,B,TENSOR_EVAL,xe,xq[0][0]);CHKERRQ(ierr);
    ierr = PetscMemzero(dx,sizeof dx);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,D,B,B,TENSOR_EVAL,xe,dx[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,B,D,B,TENSOR_EVAL,xe,dx[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,B,B,D,TENSOR_EVAL,xe,dx[2][0][0]);CHKERRQ(ierr);
    ierr = PointwiseJacobianInvert(ne,Q3,w3,dx,wdxdet);CHKERRQ(ierr);

    ierr = DMFEExtractElements(dm,u,e,ne,ue);CHKERRQ(ierr);
    ierr = PetscMemzero(uq,sizeof uq);CHKERRQ(ierr);
    ierr = TensorContract(op->TensorDOF,B,B,B,TENSOR_EVAL,ue,uq[0][0]);CHKERRQ(ierr);

    for (PetscInt i=0; i<Q3; i++) {
      for (PetscInt l=0; l<ne; l++) {
        PetscReal xx[] = {xq[0][i][l],xq[1][i][l],xq[2][i][l]};
        PetscScalar uql[op->dof],fql[op->dof];
        ierr = (op->PointwiseSolution)(op,xx,L,uql);CHKERRQ(ierr);
        ierr = (op->PointwiseForcing)(op,xx,L,fql);CHKERRQ(ierr);
        for (PetscInt d=0; d<op->dof; d++) {
          PetscReal error = uq[d][i][l] - uql[d];
          sumInfty.error = PetscMax(sumInfty.error,PetscAbs(error));
          sumInfty.u     = PetscMax(sumInfty.u    ,PetscAbs(uql[d]));
          sum2.error    += PetscSqr(error) * wdxdet[i][l];
          sum2.u        += PetscSqr(uql[d]) * wdxdet[i][l];
        }
      }
    }
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(Uloc,&u);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Uloc);CHKERRQ(ierr);
  ierr = MPI_Allreduce(MPI_IN_PLACE,(void*)&sumInfty,2,MPIU_REAL,MPIU_MAX,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
  ierr = MPI_Allreduce(MPI_IN_PLACE,(void*)&sum2,2,MPIU_REAL,MPIU_SUM,PetscObjectComm((PetscObject)dm));CHKERRQ(ierr);
  *normInfty = sumInfty.error/sumInfty.u;
  *normL2    = PetscSqrtReal(sum2.error)/PetscSqrtReal(sum2.u);
  ierr = PetscLogEventEnd(OP_IntegrateNorms,dm,U,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpApply(Op op,DM dm,Vec U,Vec F) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!op->Apply) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"No Apply implemented, use OpSetApply()");
  ierr = PetscLogEventBegin(OP_Apply,dm,U,F,0);CHKERRQ(ierr);
  ierr = (*op->Apply)(op,dm,U,F);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(OP_Apply,dm,U,F,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpGetDiagonal(Op op,DM dm,Vec Diag) {
  PetscErrorCode ierr;
  Vec X,Vl;
  PetscInt nelem,P,Q,P3,Q3,NE;
  DM dmx;
  const PetscScalar *x;
  PetscScalar *diag;
  const PetscReal *B,*D,*w3;

  PetscFunctionBegin;
  if (!op->PointwiseElement) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_USER,"No PointwiseElement implemented, use OpSetPointwiseElement()");
  if (op->dof != 1) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"dof != 1");
  ierr = PetscLogEventBegin(OP_GetDiagonal,dm,Diag,0,0);CHKERRQ(ierr);
  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL,&w3);CHKERRQ(ierr);
  P3 = P*P*P;
  Q3 = Q*Q*Q;
  NE = op->ne;

  ierr = DMGetLocalVector(dm,&Vl);CHKERRQ(ierr);
  ierr = VecZeroEntries(Vl);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(dm,&dmx);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  ierr = DMFEGetNumElements(dm,&nelem);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecGetArray(Vl,&diag);CHKERRQ(ierr);

  for (PetscInt e=0; e<nelem; e+=NE) {
    PetscScalar diage[1*P3*NE]_align,ve[1*P3*NE]_align,vv[1][1][Q3][NE]_align,dv[3][1][Q3][NE]_align,ue[1*P3*NE]_align,uu[1][1][Q3][NE]_align,du[3][1][Q3][NE]_align,xe[3*P3*NE]_align,dx[3][3][Q3][NE]_align,wdxdet[Q3][NE]_align;

    ierr = DMFEExtractElements(dmx,x,e,NE,xe);CHKERRQ(ierr);
    ierr = PetscMemzero(dx,sizeof dx);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,D,B,B,TENSOR_EVAL,xe,dx[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,B,D,B,TENSOR_EVAL,xe,dx[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(op->Tensor3,B,B,D,TENSOR_EVAL,xe,dx[2][0][0]);CHKERRQ(ierr);
    ierr = PointwiseJacobianInvert(NE,Q*Q*Q,w3,dx,wdxdet);CHKERRQ(ierr);

    for (PetscInt i=0; i<P3; i++) {
      ierr = PetscMemzero(ue,sizeof ue);CHKERRQ(ierr);
      for (PetscInt k=0; k<op->ne; k++) ue[i*NE+k] = 1;
      if (op->opcode & 01) {
        ierr = PetscMemzero(uu,sizeof uu);CHKERRQ(ierr);
        ierr = TensorContract(op->TensorDOF,B,B,B,TENSOR_EVAL,ue,uu[0][0][0]);CHKERRQ(ierr);
      }
      if (op->opcode & 02) {
        ierr = PetscMemzero(du,sizeof du);CHKERRQ(ierr);
        ierr = TensorContract(op->TensorDOF,D,B,B,TENSOR_EVAL,ue,du[0][0][0]);CHKERRQ(ierr);
        ierr = TensorContract(op->TensorDOF,B,D,B,TENSOR_EVAL,ue,du[1][0][0]);CHKERRQ(ierr);
        ierr = TensorContract(op->TensorDOF,B,B,D,TENSOR_EVAL,ue,du[2][0][0]);CHKERRQ(ierr);
      }
      ierr = (*op->PointwiseElement)(op,NE,Q3,dx[0][0][0],wdxdet[0],uu[0][0][0],du[0][0][0],vv[0][0][0],dv[0][0][0]);CHKERRQ(ierr);
      ierr = PetscMemzero(ve,sizeof ve);CHKERRQ(ierr);
      if (op->opcode & 010) {
        ierr = TensorContract(op->TensorDOF,B,B,B,TENSOR_TRANSPOSE,vv[0][0][0],ve);CHKERRQ(ierr);
      }
      if (op->opcode & 020) {
        ierr = TensorContract(op->TensorDOF,D,B,B,TENSOR_TRANSPOSE,dv[0][0][0],ve);CHKERRQ(ierr);
        ierr = TensorContract(op->TensorDOF,B,D,B,TENSOR_TRANSPOSE,dv[1][0][0],ve);CHKERRQ(ierr);
        ierr = TensorContract(op->TensorDOF,B,B,D,TENSOR_TRANSPOSE,dv[2][0][0],ve);CHKERRQ(ierr);
      }
      for (PetscInt k=0; k<op->ne; k++) diage[i*NE+k] = ve[i*NE+k];
    }
    ierr = DMFESetElements(dm,diag,e,NE,ADD_VALUES,DOMAIN_INTERIOR,diage);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(Vl,&diag);CHKERRQ(ierr);
  ierr = VecZeroEntries(Diag);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,Vl,ADD_VALUES,Diag);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,Vl,ADD_VALUES,Diag);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Vl);CHKERRQ(ierr);
  ierr = PetscLogEventEnd(OP_GetDiagonal,dm,Diag,0,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

typedef struct {
  Op op;
  DM dm;
} Mat_Op;

static PetscErrorCode MatDestroy_Op(Mat A) {
  PetscErrorCode ierr;
  Mat_Op *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,&ctx);CHKERRQ(ierr);
  ierr = PetscFree(ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatMult_Op(Mat A,Vec X,Vec Y) {
  PetscErrorCode ierr;
  Mat_Op *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,&ctx);CHKERRQ(ierr);
  ierr = OpApply(ctx->op,ctx->dm,X,Y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode MatGetDiagonal_Op(Mat A,Vec D) {
  PetscErrorCode ierr;
  Mat_Op *ctx;

  PetscFunctionBegin;
  ierr = MatShellGetContext(A,&ctx);CHKERRQ(ierr);
  ierr = OpGetDiagonal(ctx->op,ctx->dm,D);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpGetMat(Op op,DM dm,Mat *shell) {
  PetscErrorCode ierr;
  Mat A;
  Vec U;
  PetscInt m;
  Mat_Op *ctx;

  PetscFunctionBegin;
  ierr = DMGetGlobalVector(dm,&U);CHKERRQ(ierr);
  ierr = VecGetLocalSize(U,&m);CHKERRQ(ierr);
  ierr = DMRestoreGlobalVector(dm,&U);CHKERRQ(ierr);
  ierr = PetscMalloc1(1,&ctx);CHKERRQ(ierr);
  // Danger: no reference counting
  ctx->op = op;
  ctx->dm = dm;
  ierr = MatCreateShell(PetscObjectComm((PetscObject)dm),m,m,PETSC_DETERMINE,PETSC_DETERMINE,ctx,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_DESTROY,(void(*)(void))MatDestroy_Op);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatMult_Op);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonal_Op);CHKERRQ(ierr);
  *shell = A;
  PetscFunctionReturn(0);
}

PetscErrorCode OpGetTensors(Op op,Tensor *TensorDOF,Tensor *Tensor3) {
  PetscFunctionBegin;
  if (TensorDOF) *TensorDOF = op->TensorDOF;
  if (Tensor3)   *Tensor3   = op->Tensor3;
  PetscFunctionReturn(0);
}

PetscErrorCode OpRestrictState(Op op,DM dm,Vec Uf,Vec Uc) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(OP_RestrictState,dm,Uf,Uc,0);CHKERRQ(ierr);
  if (op->RestrictState) {
    ierr = (*op->RestrictState)(op,dm,Uf,Uc);CHKERRQ(ierr);
  } else {
    ierr = DMFEInject(dm,Uf,Uc);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(OP_RestrictState,dm,Uf,Uc,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpRestrictResidual(Op op,DM dm,Vec Uf,Vec Uc) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(OP_RestrictResidual,dm,Uf,Uc,0);CHKERRQ(ierr);
  if (op->RestrictResidual) {
    ierr = (*op->RestrictResidual)(op,dm,Uf,Uc);CHKERRQ(ierr);
  } else {
    ierr = DMFERestrict(dm,Uf,Uc);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(OP_RestrictResidual,dm,Uf,Uc,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpInterpolate(Op op,DM dm,Vec Uc,Vec Uf) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogEventBegin(OP_Interpolate,dm,Uc,Uf,0);CHKERRQ(ierr);
  if (op->Interpolate) {
    ierr = (*op->Interpolate)(op,dm,Uc,Uf);CHKERRQ(ierr);
  } else {
    ierr = DMFEInterpolate(dm,Uc,Uf);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(OP_Interpolate,dm,Uc,Uf,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpCreateFromOptions(MPI_Comm comm,Op *op)
{
  PetscErrorCode ierr,(*f)(Op);
  Op o;
  char opname[256] = "poisson2";
  PetscInt addquadpts=0;

  PetscFunctionBegin;
  ierr = OpInitializePackage();CHKERRQ(ierr);
  ierr = PetscNew(&o);CHKERRQ(ierr);
  ierr = PetscCommDuplicate(comm,&o->comm,NULL);CHKERRQ(ierr);
  ierr = OpSetFEDegree(o,1);CHKERRQ(ierr);
  ierr = OpSetDof(o,1);CHKERRQ(ierr);
  ierr = PetscOptionsBegin(comm,NULL,"Operator options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsFList("-op_type","Operator type","",OpList,opname,opname,sizeof opname,NULL);CHKERRQ(ierr);
  ierr = PetscFunctionListFind(OpList,opname,&f);CHKERRQ(ierr);
  if (!f) SETERRQ1(comm,PETSC_ERR_USER,"Operator type '%s' not found",opname);
  ierr = (*f)(o);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-add_quad_pts","Number of additional quadrature points","",addquadpts,&addquadpts,NULL);CHKERRQ(ierr);
  ierr = OpSetAddQuadPts(o,addquadpts); CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  ierr = TensorCreate(o->ne,o->dof,o->fedegree+1,o->fedegree+1+o->addquadpts,&o->TensorDOF);CHKERRQ(ierr);
  ierr = TensorCreate(o->ne,3,o->fedegree+1,o->fedegree+1+o->addquadpts,&o->Tensor3);CHKERRQ(ierr);
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
  ierr = TensorDestroy(&(*op)->TensorDOF);CHKERRQ(ierr);
  ierr = TensorDestroy(&(*op)->Tensor3);CHKERRQ(ierr);
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
  ierr = PetscLogEventRegister("OpApply"         ,DM_CLASSID,&OP_Apply);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OpRestrictState" ,DM_CLASSID,&OP_RestrictState);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OpRestrictResid" ,DM_CLASSID,&OP_RestrictResidual);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OpInterpolate"   ,DM_CLASSID,&OP_Interpolate);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OpForcing"       ,DM_CLASSID,&OP_Forcing);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OpSolution"      ,DM_CLASSID,&OP_Solution);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OpIntegNorms"    ,DM_CLASSID,&OP_IntegrateNorms);CHKERRQ(ierr);
  ierr = PetscLogEventRegister("OpGetDiagonal"   ,DM_CLASSID,&OP_GetDiagonal);CHKERRQ(ierr);
  ierr = PetscRegisterFinalize(OpFinalizePackage);CHKERRQ(ierr);
  OpPackageInitialized = PETSC_TRUE;
  PetscFunctionReturn(0);
}
