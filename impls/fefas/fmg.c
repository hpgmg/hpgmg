#include "fefas.h"
#include "op/fefas-op.h"
#include <petscksp.h>

typedef struct Options_private *Options;
struct Options_private {
  PetscInt M[3];
  PetscInt p[3];
  PetscInt cmax;
  PetscReal L[3];
};

typedef struct MG_private *MG;
struct MG_private {
  DM dm;
  KSP ksp;
  MG coarse;
};

static PetscErrorCode OptionsParse(const char *header,Options *opt)
{
  PetscErrorCode ierr;
  Options o;
  PetscInt three,M_max;

  PetscFunctionBegin;
  *opt = NULL;
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

static PetscErrorCode MGCreate(Op op,DM dm,PetscInt nlevels,MG *newmg) {
  PetscErrorCode ierr;
  MG mg;

  PetscFunctionBegin;
  ierr = PetscNew(&mg);CHKERRQ(ierr);
  mg->dm = dm;
  *newmg = mg;
  for (PetscInt lev=nlevels-1; ; lev--) {
    DM dmcoarse;
    if (mg->dm) { // I have some grid at this level
      Mat A;
      PC pc;
      char prefix[256];
      ierr = KSPCreate(PetscObjectComm((PetscObject)mg->dm),&mg->ksp);CHKERRQ(ierr);
      ierr = KSPGetPC(mg->ksp,&pc);CHKERRQ(ierr);
      ierr = KSPSetConvergenceTest(mg->ksp,KSPConvergedSkip,NULL,NULL);CHKERRQ(ierr);
      ierr = KSPSetNormType(mg->ksp,KSP_NORM_NONE);CHKERRQ(ierr);
      ierr = KSPSetType(mg->ksp,KSPCHEBYSHEV);CHKERRQ(ierr);
      ierr = KSPChebyshevSetEigenvalues(mg->ksp,2,0.2);CHKERRQ(ierr);
      ierr = KSPSetInitialGuessNonzero(mg->ksp,PETSC_TRUE);CHKERRQ(ierr);
      ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);
      ierr = OpGetMat(op,mg->dm,&A);CHKERRQ(ierr);
      ierr = KSPSetOperators(mg->ksp,A,A);CHKERRQ(ierr);
      ierr = MatDestroy(&A);CHKERRQ(ierr);
      ierr = PetscSNPrintf(prefix,sizeof prefix,"mg_%D_",lev);CHKERRQ(ierr);
      ierr = KSPSetOptionsPrefix(mg->ksp,prefix);CHKERRQ(ierr);
      ierr = KSPSetFromOptions(mg->ksp);CHKERRQ(ierr);
    }
    if (lev == 0) break;
    if (mg->dm) {
      ierr = DMFECoarsen(mg->dm,&dmcoarse);CHKERRQ(ierr);
    } else dmcoarse = NULL;
    MPI_Barrier(PETSC_COMM_WORLD);
    ierr = PetscNew(&mg->coarse);CHKERRQ(ierr);
    mg->coarse->dm = dmcoarse;
    mg = mg->coarse;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MGDestroy(MG *mg) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*mg) PetscFunctionReturn(0);
  ierr = MGDestroy(&(*mg)->coarse);CHKERRQ(ierr);
  ierr = KSPDestroy(&(*mg)->ksp);CHKERRQ(ierr);
  ierr = PetscFree(*mg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Written in FAS form
//   Ac uc = R bf + Ac Rhat uf - R Af uf
// Collecting the affine terms:
//   Ac uc = R (bf - Af uf) + Ac Rhat uf
static PetscErrorCode MGVCycle(Op op,MG mg,PetscInt presmooths,PetscInt postsmooths,Vec B,Vec U) {
  PetscErrorCode ierr;
  DM dm = mg->dm,dmcoarse;
  Vec V,Vc,Uc;

  PetscFunctionBegin;
  ierr = KSPSetTolerances(mg->ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,presmooths);CHKERRQ(ierr);
  ierr = KSPSolve(mg->ksp,B,U);CHKERRQ(ierr);
  if (!mg->coarse) PetscFunctionReturn(0);

  dmcoarse = mg->coarse->dm;
  ierr = DMGetGlobalVector(dm,&V);CHKERRQ(ierr);
  ierr = OpApply(op,dm,U,V);CHKERRQ(ierr);
  ierr = VecAYPX(V,-1,B);CHKERRQ(ierr); // b - A u
  if (dmcoarse) {
    ierr = DMGetGlobalVector(dmcoarse,&Uc);CHKERRQ(ierr);
    ierr = DMGetGlobalVector(dmcoarse,&Vc);CHKERRQ(ierr);
  } else {
    Uc = NULL;
    Vc = NULL;
  }
  ierr = OpRestrictState(op,dm,U,Uc);CHKERRQ(ierr);    // Rhat u
  if (dmcoarse) {
    ierr = OpApply(op,dmcoarse,Uc,Vc);CHKERRQ(ierr);
  }
  ierr = OpRestrictResidual(op,dm,V,Vc);CHKERRQ(ierr); // R (b - A u) + Ac (Rhat u)

  if (dmcoarse) {
    Vec Yc;
    ierr = DMFEZeroBoundaries(dmcoarse,Vc);CHKERRQ(ierr);
    ierr = DMGetGlobalVector(dmcoarse,&Yc);CHKERRQ(ierr);
    ierr = VecCopy(Uc,Yc);CHKERRQ(ierr);
    ierr = MGVCycle(op,mg->coarse,presmooths,postsmooths,Vc,Uc);CHKERRQ(ierr);
    ierr = VecAXPY(Uc,-1.,Yc);CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(dmcoarse,&Yc);CHKERRQ(ierr);
  }
  ierr = OpInterpolate(op,dm,Uc,V);CHKERRQ(ierr);
  ierr = VecAXPY(U,1,V);CHKERRQ(ierr);
  if (dmcoarse) {
    ierr = DMRestoreGlobalVector(dmcoarse,&Uc);CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(dmcoarse,&Vc);CHKERRQ(ierr);
  }
  ierr = DMRestoreGlobalVector(dm,&V);CHKERRQ(ierr);

  ierr = KSPSetTolerances(mg->ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,presmooths);CHKERRQ(ierr);
  ierr = KSPSolve(mg->ksp,B,U);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RunMGV()
{
  PetscErrorCode ierr;
  Grid grid;
  Options opt;
  Op op;
  PetscInt fedegree,dof,nlevels;
  DM dm;
  Vec U0,U,F;
  MG mg;
  PetscReal normU0,normError,normF,normResid;

  PetscFunctionBegin;
  ierr = OpCreateFromOptions(PETSC_COMM_WORLD,&op);CHKERRQ(ierr);
  ierr = OpGetFEDegree(op,&fedegree);CHKERRQ(ierr);
  ierr = OpGetDof(op,&dof);CHKERRQ(ierr);
  ierr = OptionsParse("Finite Element FAS FMG solver",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,NULL,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = GridGetNumLevels(grid,&nlevels);CHKERRQ(ierr);

  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,dof,&dm);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&U0);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&U);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&F);CHKERRQ(ierr);
  ierr = OpForcing(op,dm,F);CHKERRQ(ierr);
  ierr = OpSolution(op,dm,U0);CHKERRQ(ierr);
  ierr = VecNorm(U0,NORM_2,&normU0);CHKERRQ(ierr);
  ierr = VecNorm(F,NORM_2,&normF);CHKERRQ(ierr);

  ierr = MGCreate(op,dm,nlevels,&mg);CHKERRQ(ierr);
  for (PetscInt i=0; i<5; i++) {
    PetscInt presmooths = 1,postsmooths = 1;
    Vec Y;
    ierr = MGVCycle(op,mg,presmooths,postsmooths,F,U);CHKERRQ(ierr);
    ierr = DMGetGlobalVector(dm,&Y);CHKERRQ(ierr);
    ierr = VecCopy(U,Y);CHKERRQ(ierr);
    ierr = VecAXPY(Y,-1,U0);CHKERRQ(ierr);
    ierr = VecNorm(Y,NORM_2,&normError);CHKERRQ(ierr);
    ierr = OpApply(op,dm,U,Y);CHKERRQ(ierr);
    ierr = VecAYPX(Y,-1,F);CHKERRQ(ierr);
    ierr = VecNorm(Y,NORM_2,&normResid);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"V(%D,%D) %2D: |e|_2/|u|_2 %8.2e  |r|_2/|f|_2 %8.2e\n",presmooths,postsmooths,i,(double)(normError/normU0),(double)(normResid/normF));CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(dm,&Y);CHKERRQ(ierr);
  }

  ierr = MGDestroy(&mg);CHKERRQ(ierr);
  ierr = VecDestroy(&U0);CHKERRQ(ierr);
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = OpDestroy(&op);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
