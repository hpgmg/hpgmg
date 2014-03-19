#include "fefas.h"
#include "tensor.h"

typedef struct Options_private *Options;
struct Options_private {
  char command[256];
  PetscInt M[3];
  PetscInt p[3];
  PetscInt cmax;
  PetscReal L[3];
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
  o->cmax = 100;
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
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  *opt = o;
  PetscFunctionReturn(0);
}

PetscErrorCode TestGrid()
{
  PetscErrorCode ierr;
  Grid grid;
  Options opt;

  PetscFunctionBegin;
  ierr = OptionsParse("Finite Element FAS Test Grid",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,NULL,opt->cmax,&grid);CHKERRQ(ierr);
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
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,NULL,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,1,1,&dm);CHKERRQ(ierr);
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

  PetscFunctionBegin;
  ierr = OptionsParse("Finite Element FAS Test Element Gradients",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,NULL,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,1,&dm);CHKERRQ(ierr);
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

  ierr = DMFEGetTensorEval(dm,&P,&Q,&B,&D,NULL,NULL);CHKERRQ(ierr);
  ierr = DMFEGetNumElements(dm,&nelems);CHKERRQ(ierr);
  ierr = PetscMalloc2(P*P*P*ne,&ue,3*Q*Q*Q*ne,&du);CHKERRQ(ierr);

  ierr = VecGetArrayRead(L,&l);CHKERRQ(ierr);
  for (PetscInt e=0; e<nelems; e+=ne) {
    ierr = DMFEExtractElements(dm,l,e,ne,ue);CHKERRQ(ierr);
    ierr = PetscMemzero(du,3*Q*Q*Q*ne*sizeof(*du));CHKERRQ(ierr);
    ierr = TensorContract(ne,1,P,Q,D,B,B,TENSOR_EVAL,ue,&du[0*Q*Q*Q*ne]);CHKERRQ(ierr);
    ierr = TensorContract(ne,1,P,Q,B,D,B,TENSOR_EVAL,ue,&du[1*Q*Q*Q*ne]);CHKERRQ(ierr);
    ierr = TensorContract(ne,1,P,Q,B,B,D,TENSOR_EVAL,ue,&du[2*Q*Q*Q*ne]);CHKERRQ(ierr);
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
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,NULL,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,1,&dm);CHKERRQ(ierr);
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
    // ierr = VecView(Gc,NULL);CHKERRQ(ierr);
    ierr = VecDestroy(&Gc);CHKERRQ(ierr);
  } else {
    ierr = DMFEInject(dm,G,NULL);CHKERRQ(ierr);
  }
  ierr = VecDestroy(&G);CHKERRQ(ierr);
  ierr = VecDestroy(&L);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = DMDestroy(&dmcoarse);CHKERRQ(ierr);
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
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,NULL,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,1,&dm);CHKERRQ(ierr);
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
  ierr = DMDestroy(&dmcoarse);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
