#include "fefas.h"
#include "tensor.h"
#include "pointwise.h"
#include "op/fefas-op.h"
#include <petscksp.h>
#include <inttypes.h>

typedef struct Options_private *Options;
struct Options_private {
  PetscInt M[3];
  PetscInt p[3];
  PetscInt cmax;
  PetscReal L[3];
  PetscInt addquadpts;
};

static PetscErrorCode OptionsParse(const char *header,Options *opt)
{
  PetscErrorCode ierr;
  Options o;
  PetscInt three,M_max;

  PetscFunctionBegin;
  ierr = PetscNew(&o);CHKERRQ(ierr);
  o->M[0] = 10;
  o->M[1] = 10;
  o->M[2] = 10;
  o->p[0] = 1;
  o->p[1] = 1;
  o->p[2] = 1;
  o->cmax = 3*4*4;
  o->addquadpts = 0;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,header,NULL);CHKERRQ(ierr);
  three = 3;
  ierr = PetscOptionsIntArray("-M","Fine grid dimensions","",o->M,&three,NULL);CHKERRQ(ierr);
  three = 3;
  ierr = PetscOptionsIntArray("-p","Process grid dimensions","",o->p,&three,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-cmax","Max coarse grid size","",o->cmax,&o->cmax,NULL);CHKERRQ(ierr);
  M_max = PetscMax(o->M[0],PetscMax(o->M[1],o->M[2]));
  o->L[0] = o->M[0]*1./M_max;
  o->L[1] = o->M[1]*1./M_max;
  o->L[2] = o->M[2]*1./M_max;
  three = 3;
  ierr = PetscOptionsRealArray("-L","Grid dimensions","",o->L,&three,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-add_quad_pts","Additional Quadrature Points","",o->addquadpts,&o->addquadpts,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  *opt = o;
  PetscFunctionReturn(0);
}

PetscErrorCode TestAddQuad()
{
  PetscErrorCode ierr;
  Grid grid;
  DM dm;
  Vec WF;
  PetscInt fedegree=1,P,Q;
  PetscScalar *wf,quadint,trueval=14.7344823077179;
  const PetscReal *x,*w3;
  Options opt;
  
  PetscFunctionBegin;
  
  ierr = OptionsParse("Finite Element Test Additional Quadrature Points",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,1,opt->addquadpts,&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = DMFEGetTensorEval(dm,&P,&Q,NULL,NULL,&x,NULL,&w3);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&WF);CHKERRQ(ierr);
  ierr = VecSetType(WF,VECMPI);CHKERRQ(ierr);
  ierr = VecSetSizes(WF,PETSC_DECIDE,Q*Q*Q);CHKERRQ(ierr);
  ierr = VecGetArray(WF,&wf);CHKERRQ(ierr);

  for (PetscInt i=0; i<Q; i++) {
    for (PetscInt j=0; j<Q; j++) {
      for (PetscInt k=0; k<Q; k++) {
        wf[(i*Q+j)*Q+k]=w3[(i*Q+j)*Q+k]/(0.001+PetscSqr(x[i])+PetscSqr(x[j])+PetscSqr(x[k]));
      }
    }
  }
  ierr = VecRestoreArray(WF,&wf);CHKERRQ(ierr);
  ierr = VecSum(WF,&quadint);CHKERRQ(ierr);

  if (PetscAbs(quadint - trueval) > 1e-6) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"quadrature = %f expected %f\n",(double)quadint,(double)trueval);CHKERRQ(ierr);
  }

  ierr = VecDestroy(&WF);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestGrid()
{
  PetscErrorCode ierr;
  Grid grid;
  Options opt;

  PetscFunctionBegin;
  ierr = OptionsParse("Finite Element FAS Test Grid",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestFESpace()
{
  PetscErrorCode ierr;
  Grid grid;
  DM dm;
  Vec G,L;
  PetscInt i,rstart,rend;
  PetscScalar *g;
  Options opt;

  PetscFunctionBegin;
  ierr = OptionsParse("Finite Element FAS Test FE global-to-local",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,1,1,0,&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&G);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(dm,&L);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange(G,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetArray(G,&g);CHKERRQ(ierr);
  for (i=rstart; i<rend; i++) g[i-rstart] = i;
  ierr = VecRestoreArray(G,&g);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm,G,INSERT_VALUES,L);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm,G,INSERT_VALUES,L);CHKERRQ(ierr);
  ierr = VecView(L,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = VecDestroy(&G);CHKERRQ(ierr);
  ierr = VecDestroy(&L);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestFEGrad()
{
  PetscErrorCode ierr;
  Grid grid;
  DM dm;
  Vec X,L;
  PetscInt fedegree = 1,P,Q,ne = 2,m,nelems;
  const PetscReal *B,*D,TestGrad[] = {2,3,5};
  PetscScalar *u,*ue,*du;
  const PetscScalar *l,*x;
  Options opt;
  Tensor Tensor1;

  PetscFunctionBegin;
  ierr = OptionsParse("Finite Element FAS Test Element Gradients",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,1,0,&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);

  ierr = DMCreateLocalVector(dm,&L);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  ierr = VecGetLocalSize(L,&m);CHKERRQ(ierr);
  ierr = VecGetArray(L,&u);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  for (PetscInt i=0; i<m; i++) {
    u[i] = TestGrad[0]*x[i*3+0] + TestGrad[1]*x[i*3+1] + TestGrad[2]*x[i*3+2];
  }
  ierr = VecRestoreArray(L,&u);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL,NULL);CHKERRQ(ierr);
  ierr = TensorCreate(ne,1,P,Q,&Tensor1);CHKERRQ(ierr);
  ierr = DMFEGetNumElements(dm,&nelems);CHKERRQ(ierr);
  ierr = PetscMalloc2(P*P*P*ne,&ue,3*Q*Q*Q*ne,&du);CHKERRQ(ierr);

  ierr = VecGetArrayRead(L,&l);CHKERRQ(ierr);
  for (PetscInt e=0; e<nelems; e+=ne) {
    ierr = DMFEExtractElements(dm,l,e,ne,ue);CHKERRQ(ierr);
    ierr = PetscMemzero(du,3*Q*Q*Q*ne*sizeof(*du));CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,D,B,B,TENSOR_EVAL,ue,&du[0*Q*Q*Q*ne]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,D,B,TENSOR_EVAL,ue,&du[1*Q*Q*Q*ne]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,B,D,TENSOR_EVAL,ue,&du[2*Q*Q*Q*ne]);CHKERRQ(ierr);
    for (PetscInt ee=0; ee<ne; ee++) {
      for (PetscInt i=0; i<Q; i++) {
        for (PetscInt j=0; j<Q; j++) {
          for (PetscInt k=0; k<Q; k++) {
            PetscInt q = ((i*Q+j)*Q+k)*ne+ee;
            for (PetscInt d=0; d<3; d++) {
              PetscScalar dux = du[d*Q*Q*Q*ne+q]*2*opt->M[d]/opt->L[d];
              if (PetscAbs(dux - TestGrad[d]) > 1e-12) {
                ierr = PetscPrintf(PETSC_COMM_WORLD,"GradU[elem %D][qp %D,%D,%D][%D] = %f expected %f\n",e+ee,i,j,k,d,(double)dux,(double)TestGrad[d]);CHKERRQ(ierr);
              }
            }
          }
        }
      }
    }
  }
  ierr = VecRestoreArrayRead(L,&l);CHKERRQ(ierr);
  ierr = PetscFree2(ue,du);CHKERRQ(ierr);
  ierr = TensorDestroy(&Tensor1);CHKERRQ(ierr);
  ierr = VecDestroy(&L);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestFEInject()
{
  PetscErrorCode ierr;
  Grid grid;
  DM dm,dmcoarse;
  Vec X,L,G;
  PetscInt fedegree = 1,m;
  const PetscReal TestGrad[] = {1000000,1000,1};
  PetscScalar *u;
  const PetscScalar *x;
  Options opt;

  PetscFunctionBegin;
  ierr = OptionsParse("Finite Element FAS Test Injection (state restriction)",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,1,0,&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);

  ierr = DMCreateLocalVector(dm,&L);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  ierr = VecGetLocalSize(L,&m);CHKERRQ(ierr);
  ierr = VecGetArray(L,&u);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  for (PetscInt i=0; i<m; i++) {
    u[i] = TestGrad[0]*x[i*3+0] + TestGrad[1]*x[i*3+1] + TestGrad[2]*x[i*3+2];
  }
  ierr = VecRestoreArray(L,&u);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&G);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,L,INSERT_VALUES,G);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,L,INSERT_VALUES,G);CHKERRQ(ierr);

  ierr = DMFECoarsen(dm,&dmcoarse);CHKERRQ(ierr);
  if (dmcoarse) {
    Vec Gc;
    const PetscScalar *u;
    ierr = DMCreateGlobalVector(dmcoarse,&Gc);CHKERRQ(ierr);
    ierr = DMFEInject(dm,G,Gc);CHKERRQ(ierr);
    ierr = DMGetCoordinates(dmcoarse,&X);CHKERRQ(ierr);
    ierr = VecGetLocalSize(Gc,&m);CHKERRQ(ierr);
    ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
    ierr = VecGetArrayRead(Gc,&u);CHKERRQ(ierr);
    for (PetscInt i=0; i<m; i++) {
      ierr = PetscPrintf(PetscObjectComm((PetscObject)dmcoarse),"coarse u[%2D] = %10.1f at %4.1f %4.1f %4.1f\n",i,u[i],x[i*3+0],x[i*3+1],x[i*3+2]);CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(Gc,&u);CHKERRQ(ierr);
    ierr = VecDestroy(&Gc);CHKERRQ(ierr);
  } else {
    ierr = DMFEInject(dm,G,NULL);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&G);CHKERRQ(ierr);
  ierr = VecDestroy(&L);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestFEInterp()
{
  PetscErrorCode ierr;
  Grid grid;
  DM dm,dmcoarse;
  Vec X,G,G2,Gc;
  PetscInt fedegree = 1,m;
  const PetscReal TestGrad[] = {1000000,1000,1};
  PetscScalar *u;
  PetscReal norm;
  const PetscScalar *x;
  Options opt;

  PetscFunctionBegin;
  ierr = OptionsParse("Finite Element FAS Test Interpolation",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,1,0,&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&G);CHKERRQ(ierr);
  ierr = DMGetCoordinates(dm,&X);CHKERRQ(ierr);
  ierr = VecGetLocalSize(G,&m);CHKERRQ(ierr);
  ierr = VecGetArray(G,&u);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  for (PetscInt i=0; i<m; i++) {
    u[i] = TestGrad[0]*x[i*3+0] + TestGrad[1]*x[i*3+1] + TestGrad[2]*x[i*3+2];
  }
  ierr = VecRestoreArray(G,&u);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&G2);CHKERRQ(ierr);
  ierr = DMFECoarsen(dm,&dmcoarse);CHKERRQ(ierr);
  if (dmcoarse) {
    ierr = DMCreateGlobalVector(dmcoarse,&Gc);CHKERRQ(ierr);
  } else Gc = NULL;
  ierr = DMFEInject(dm,G,Gc);CHKERRQ(ierr);
  ierr = DMFEInterpolate(dm,Gc,G2);CHKERRQ(ierr);
  if (0) {
    const PetscScalar *u,*v;
    ierr = VecGetArrayRead(G,&u);CHKERRQ(ierr);
    ierr = VecGetArrayRead(G2,&v);CHKERRQ(ierr);
    ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
    for (PetscInt i=0; i<m; i++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"u[%2d] = %10.1f vs %10.1f at %4.1f %4.1f %4.1f\n",i,u[i],v[i],x[i*3+0],x[i*3+1],x[i*3+2]);CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(G,&u);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(G2,&v);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  }
  ierr = VecAXPY(G2,-1.,G);CHKERRQ(ierr);
  ierr = VecNorm(G2,NORM_MAX,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|u - I Ihat u|_max = %5.1g\n",norm);CHKERRQ(ierr);

  ierr = VecDestroy(&Gc);CHKERRQ(ierr);
  ierr = VecDestroy(&G);CHKERRQ(ierr);
  ierr = VecDestroy(&G2);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode IntegrateTestFunction(DM dm,Vec G)
{
  PetscErrorCode ierr;
  Vec L,X;
  PetscScalar *u;
  const PetscScalar *x;
  const PetscReal *B,*D,*w3;
  PetscInt nelem,ne = 2,P,Q,P3,Q3;
  DM dmx;
  Tensor Tensor1,Tensor3;

  PetscFunctionBegin;
  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL,&w3);CHKERRQ(ierr);
  P3 = P*P*P;
  Q3 = Q*Q*Q;
  ierr = TensorCreate(ne,1,P,Q,&Tensor1);CHKERRQ(ierr);
  ierr = TensorCreate(ne,3,P,Q,&Tensor3);CHKERRQ(ierr);

  ierr = DMGetLocalVector(dm,&L);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  ierr = DMGetCoordinateDM(dm,&dmx);CHKERRQ(ierr);
  ierr = DMFEGetNumElements(dm,&nelem);CHKERRQ(ierr);
  ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecZeroEntries(L);CHKERRQ(ierr);
  ierr = VecGetArray(L,&u);CHKERRQ(ierr);
  for (PetscInt e=0; e<nelem; e+=ne) {
    PetscScalar ue[P3*ne],uq[Q3][ne],xe[3*P3*ne],xq[3][Q3][ne],dx[3][3][Q3][ne],wdxdet[Q3][ne];
    ierr = DMFEExtractElements(dmx,x,e,ne,xe);CHKERRQ(ierr);
    ierr = PetscMemzero(xq,sizeof xq);CHKERRQ(ierr);
    ierr = TensorContract(Tensor3,B,B,B,TENSOR_EVAL,xe,xq[0][0]);CHKERRQ(ierr);
    ierr = PetscMemzero(dx,sizeof dx);CHKERRQ(ierr);
    ierr = TensorContract(Tensor3,D,B,B,TENSOR_EVAL,xe,dx[0][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor3,B,D,B,TENSOR_EVAL,xe,dx[1][0][0]);CHKERRQ(ierr);
    ierr = TensorContract(Tensor3,B,B,D,TENSOR_EVAL,xe,dx[2][0][0]);CHKERRQ(ierr);

    ierr = PointwiseJacobianInvert(ne,Q*Q*Q,w3,dx,wdxdet);CHKERRQ(ierr);
    for (PetscInt i=0; i<Q3; i++) {
      for (PetscInt e=0; e<ne; e++) {
        uq[i][e] = wdxdet[i][e] * (2.*xq[0][i][e] + 3.*xq[1][i][e] + 5.*xq[2][i][e]);
      }
    }

    ierr = PetscMemzero(ue,sizeof ue);CHKERRQ(ierr);
    ierr = TensorContract(Tensor1,B,B,B,TENSOR_TRANSPOSE,uq[0],ue);CHKERRQ(ierr);
    ierr = DMFESetElements(dm,u,e,ne,ADD_VALUES,DOMAIN_CLOSURE,ue);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
  ierr = VecRestoreArray(L,&u);CHKERRQ(ierr);
  ierr = TensorDestroy(&Tensor1);CHKERRQ(ierr);
  ierr = TensorDestroy(&Tensor3);CHKERRQ(ierr);
  ierr = VecZeroEntries(G);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(dm,L,ADD_VALUES,G);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,L,ADD_VALUES,G);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&L);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestFERestrict()
{
  PetscErrorCode ierr;
  Grid grid;
  DM dm,dmcoarse;
  Vec G,Gc;
  PetscInt fedegree = 1;
  Options opt;

  PetscFunctionBegin;
  ierr = OptionsParse("Finite Element FAS Test Restriction",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,1,0,&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&G);CHKERRQ(ierr);
  ierr = IntegrateTestFunction(dm,G);CHKERRQ(ierr);

  ierr = DMFECoarsen(dm,&dmcoarse);CHKERRQ(ierr);
  if (dmcoarse) {
    ierr = DMCreateGlobalVector(dmcoarse,&Gc);CHKERRQ(ierr);
  } else Gc = NULL;
  ierr = DMFERestrict(dm,G,Gc);CHKERRQ(ierr);
  if (dmcoarse) {
    PetscReal norm;
    Vec Gc2;
    ierr = DMCreateGlobalVector(dmcoarse,&Gc2);CHKERRQ(ierr);
    ierr = IntegrateTestFunction(dmcoarse,Gc2);CHKERRQ(ierr);
    if (0) {
      Vec X;
      const PetscScalar *u,*u2,*x;
      PetscInt m;
      ierr = DMGetCoordinates(dmcoarse,&X);CHKERRQ(ierr);
      ierr = VecGetLocalSize(Gc,&m);CHKERRQ(ierr);
      ierr = VecGetArrayRead(Gc,&u);CHKERRQ(ierr);
      ierr = VecGetArrayRead(Gc2,&u2);CHKERRQ(ierr);
      ierr = VecGetArrayRead(X,&x);CHKERRQ(ierr);
      for (PetscInt i=0; i<m; i++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"u[%2d] = %10.3f vs %10.3f at %4.1f %4.1f %4.1f\n",i,u[i],u2[i],x[i*3+0],x[i*3+1],x[i*3+2]);CHKERRQ(ierr);
      }
      ierr = VecRestoreArrayRead(Gc,&u);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(Gc2,&u2);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(X,&x);CHKERRQ(ierr);
    }

    ierr = VecAXPY(Gc2,-1.,Gc);CHKERRQ(ierr);
    ierr = VecNorm(Gc2,NORM_MAX,&norm);CHKERRQ(ierr);
    if (PetscAbs(norm) < 1e-10) norm = 0;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"|u_c - I_h^H u_f|_max = %5.1g\n",norm);CHKERRQ(ierr);
    ierr = VecDestroy(&Gc2);CHKERRQ(ierr);
  }

  ierr = VecDestroy(&Gc);CHKERRQ(ierr);
  ierr = VecDestroy(&G);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestOpApply()
{
  PetscErrorCode ierr;
  Grid grid;
  Options opt;
  Op op;
  PetscInt fedegree,dof,addquadpts;
  DM dm;
  Vec U,Y,F;
  PetscReal normF,norm;

  PetscFunctionBegin;
  ierr = OpCreateFromOptions(PETSC_COMM_WORLD,&op);CHKERRQ(ierr);
  ierr = OpGetFEDegree(op,&fedegree);CHKERRQ(ierr);
  ierr = OpGetDof(op,&dof);CHKERRQ(ierr);
  ierr = OpGetAddQuadPts(op,&addquadpts);CHKERRQ(ierr);
  ierr = OptionsParse("Finite Element FAS Test operator application",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,dof,addquadpts,&dm);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&U);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&F);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&Y);CHKERRQ(ierr);
  ierr = OpSolution(op,dm,U);CHKERRQ(ierr);   // Nodally-exact solution represented in finite-element space
  ierr = OpForcing(op,dm,F);CHKERRQ(ierr);    // Manufactured right hand side that makes U solve continuum equation
  ierr = OpApply(op,dm,U,Y);CHKERRQ(ierr);    // Y <- A U

  ierr = VecAXPY(Y,-1.,F);CHKERRQ(ierr);
  ierr = VecNorm(Y,NORM_MAX,&norm);CHKERRQ(ierr);
  ierr = VecNorm(F,NORM_MAX,&normF);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|A u - F|_max/|F|_max = %g\n",(double)norm/normF);CHKERRQ(ierr);

  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&Y);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = OpDestroy(&op);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestOpDiagonal()
{
  PetscErrorCode ierr;
  Grid grid;
  Options opt;
  Op op;
  PetscInt fedegree,dof,addquadpts;
  DM dm;
  Vec Diag;
  PetscReal norm1,norm2,normMax;

  PetscFunctionBegin;
  ierr = OpCreateFromOptions(PETSC_COMM_WORLD,&op);CHKERRQ(ierr);
  ierr = OpGetFEDegree(op,&fedegree);CHKERRQ(ierr);
  ierr = OpGetDof(op,&dof);CHKERRQ(ierr);
  ierr = OpGetAddQuadPts(op,&addquadpts);CHKERRQ(ierr);
  ierr = OptionsParse("Finite Element FAS Test diagonal extraction",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,dof,addquadpts,&dm);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&Diag);CHKERRQ(ierr);
  ierr = OpGetDiagonal(op,dm,Diag);CHKERRQ(ierr);

  ierr = VecNorm(Diag,NORM_1,&norm1);CHKERRQ(ierr);
  ierr = VecNorm(Diag,NORM_2,&norm2);CHKERRQ(ierr);
  ierr = VecNorm(Diag,NORM_MAX,&normMax);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|D|_1 = %g  |D|_2 = %g  |D|_max = %g\n",(double)norm1,(double)norm2,(double)normMax);CHKERRQ(ierr);
  // ierr = VecView(Diag,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = VecDestroy(&Diag);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = OpDestroy(&op);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestKSPSolve()
{
  PetscErrorCode ierr;
  Grid grid;
  Options opt;
  Op op;
  PetscInt fedegree,dof,addquadpts;
  DM dm;
  Vec U,V,F;
  Mat A;
  KSP ksp;
  PetscReal normU,normError;

  PetscFunctionBegin;
  ierr = OpCreateFromOptions(PETSC_COMM_WORLD,&op);CHKERRQ(ierr);
  ierr = OpGetFEDegree(op,&fedegree);CHKERRQ(ierr);
  ierr = OpGetDof(op,&dof);CHKERRQ(ierr);
  ierr = OpGetAddQuadPts(op,&addquadpts);CHKERRQ(ierr);
  ierr = OptionsParse("Finite Element FAS Test KSPSolve",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,dof,addquadpts,&dm);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&U);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&V);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&F);CHKERRQ(ierr);
  ierr = OpSolution(op,dm,U);CHKERRQ(ierr);   // Nodally-exact solution represented in finite-element space
  ierr = OpForcing(op,dm,F);CHKERRQ(ierr);    // Manufactured right hand side that makes U solve continuum equation

  ierr = OpGetMat(op,dm,&A);CHKERRQ(ierr);
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,F,V);CHKERRQ(ierr);
  ierr = VecAXPY(V,-1.,U);CHKERRQ(ierr);
  ierr = VecNorm(U,NORM_2,&normU);CHKERRQ(ierr);
  ierr = VecNorm(V,NORM_2,&normError);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"|v-u|_2/|u|_2 = %g\n",normError/normU);CHKERRQ(ierr);

  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&V);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = OpDestroy(&op);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode TestSampler()
{
  PetscErrorCode ierr;
  PetscInt tmp,two = 2,maxsamples = 5,nsamples,(*gridsize)[3],squarest[3];
  PetscMPIInt nranks = -1;
  PetscReal local[2] = {48,10000};

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Sampler options",NULL);CHKERRQ(ierr);
  tmp = 1;
  ierr = PetscOptionsInt("-nranks","number of processes (to emulate)","",tmp,&tmp,NULL);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(tmp,&nranks);CHKERRQ(ierr);
  ierr = PetscOptionsRealArray("-local","range of local problem sizes","",local,&two,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-maxsamples","maximum number of samples across range","",maxsamples,&maxsamples,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  ierr = ProcessGridFindSquarest(nranks,squarest);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Processors: [%4D %4D %4D] = %d\n",squarest[0],squarest[1],squarest[2],nranks);CHKERRQ(ierr);

  ierr = SampleGridRangeCreate(nranks,(PetscReal)local[0],(PetscReal)local[1],maxsamples,&nsamples,(PetscInt**)&gridsize);CHKERRQ(ierr);
  for (PetscInt i=0; i<nsamples; i++) {
    PetscInt nlevels = GridLevelFromM(gridsize[i]);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Filtered Grid: L%2D [%4D %4D %4D] = %12"PRId64"\n",nlevels,gridsize[i][0],gridsize[i][1],gridsize[i][2],SampleGridNumElements(gridsize[i]));CHKERRQ(ierr);
  }

  ierr = PetscFree(gridsize);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
