#include "fefas.h"
#include "op/fefas-op.h"
#include <petscksp.h>

typedef struct Options_private *Options;
struct Options_private {
  PetscInt M[3];
  PetscInt p[3];
  PetscInt cmax;
  PetscReal L[3];
  PetscInt  smooth[2];
  PetscBool coord_distort;
  PetscInt addquadpts;
};

struct MG_private {
  DM dm;
  KSP ksp;
  MG coarse;
  PetscReal enormInfty,enormL2;
  PetscReal rnorm2,bnorm2;      // rnorm is normalized by bnorm
  PetscBool monitor;
  PetscBool diagnostics;
  PetscLogStage stage;
  PetscLogEvent V_Cycle;
};

static PetscErrorCode OptionsParse(const char *header,Options *opt)
{
  PetscErrorCode ierr;
  Options o;
  PetscInt three,two,M_max;

  PetscFunctionBegin;
  *opt = NULL;
  ierr = PetscNew(&o);CHKERRQ(ierr);
  o->M[0] = 10;
  o->M[1] = 10;
  o->M[2] = 10;
  o->p[0] = 1;
  o->p[1] = 1;
  o->p[2] = 1;
  // cmax is an annoyingly tricky parameter.  The smallest coarse grid we want to deal with is 3x4x4=48 elements, which
  // corresponds to 7x9x9=567 Q2 finite-element nodes.  That coarse grid size is critical because a non-multigrid method
  // needs to be able to solve fast there (and everyone is waiting for it).  Given communication latency in the
  // microsecond range, one refinement (6x8x8 elements) is sufficient to start distributing, so cmax should be smaller
  // than 3x4x4*2^3=384.
  //
  // Meanwhile, we want to support process grids like 1x1x3 for cases where the total number of processes has such a
  // factor.  The coarsening semantic requires that a 2x2x2 grid resides one process, so the grid sizes are 4x4x4 on
  // 1x1x2 and 8x8x8 on 1x1x3.  In the latter case, each process in the 1x1x3 grid will own 8*8*ceil(8/3)=192 elements.
  // We intend to ensure that no process grids be more anisotropic than 1x1x3.
  o->cmax = 192;
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
  o->smooth[0] = 2;
  o->smooth[1] = 3;
  two = 2;
  ierr = PetscOptionsIntArray("-smooth","V- and F-cycle pre,post smoothing","",o->smooth,&two,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-coord_distort","Distort coordinates within unit cube","",o->coord_distort,&o->coord_distort,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-add_quad_pts","Number of additional quadrature points","",o->addquadpts,&o->addquadpts,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  *opt = o;
  PetscFunctionReturn(0);
}

PetscErrorCode MGMonitorSet(MG mg,PetscBool mon) {
  for ( ; mg; mg=mg->coarse) mg->monitor = mon;
  return 0;
}

PetscErrorCode MGCreate(Op op,DM dm,PetscInt nlevels,MG *newmg) {
  PetscErrorCode ierr;
  MG mg;
  PetscInt two;
  PetscReal eig_target[2];
  PetscBool monitor;

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(PetscObjectComm((PetscObject)dm),NULL,"MG Options",NULL);CHKERRQ(ierr);
  two = 2;
  eig_target[0] = 1.4;
  eig_target[1] = 0.4;
  ierr = PetscOptionsRealArray("-mg_eig_target","Target max,min eigenvalues on levels","",eig_target,&two,NULL);CHKERRQ(ierr);
  monitor = PETSC_FALSE;
  ierr = PetscOptionsBool("-mg_monitor","Monitor convergence at the end of each MG cycle","",monitor,&monitor,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
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
      if (lev) {
        ierr = KSPSetConvergenceTest(mg->ksp,KSPConvergedSkip,NULL,NULL);CHKERRQ(ierr);
        ierr = KSPSetNormType(mg->ksp,KSP_NORM_NONE);CHKERRQ(ierr);
        ierr = KSPSetType(mg->ksp,KSPCHEBYSHEV);CHKERRQ(ierr);
      } else {
        ierr = KSPSetNormType(mg->ksp,KSP_NORM_NATURAL);CHKERRQ(ierr);
        ierr = KSPSetType(mg->ksp,KSPCG);CHKERRQ(ierr);
        ierr = KSPSetTolerances(mg->ksp,1e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
      }
      ierr = KSPChebyshevSetEigenvalues(mg->ksp,eig_target[0],eig_target[1]);CHKERRQ(ierr);
      ierr = KSPSetInitialGuessNonzero(mg->ksp,PETSC_TRUE);CHKERRQ(ierr);
      ierr = KSPGetPC(mg->ksp,&pc);CHKERRQ(ierr);
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
  if (monitor) {ierr = MGMonitorSet(*newmg,PETSC_TRUE);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

PetscErrorCode MGDestroy(MG *mg) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*mg) PetscFunctionReturn(0);
  ierr = MGDestroy(&(*mg)->coarse);CHKERRQ(ierr);
  ierr = KSPDestroy(&(*mg)->ksp);CHKERRQ(ierr);
  ierr = PetscFree(*mg);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// PCJacobi does not extract the diagonal until it is applied to a vector.
PetscErrorCode MGSetUpPC(MG mg) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (; mg; mg=mg->coarse) {
    if (mg->dm) {
      PC pc;
      Vec U,V;
      ierr = KSPGetPC(mg->ksp,&pc);CHKERRQ(ierr);
      ierr = DMGetGlobalVector(mg->dm,&U);CHKERRQ(ierr);
      ierr = DMGetGlobalVector(mg->dm,&V);CHKERRQ(ierr);
      ierr = VecZeroEntries(U);CHKERRQ(ierr);
      ierr = PCApply(pc,U,V);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(mg->dm,&U);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(mg->dm,&V);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscReal ConvergenceRate(PetscReal normCoarse,PetscReal normFine) {
  // Try to avoid reporting noisy rates in pre-asymptotic regime
  if (normCoarse < 1e3*PETSC_MACHINE_EPSILON && normFine > 1e3*PETSC_MACHINE_EPSILON) return 0;
  if (normCoarse == 0 || normFine == 0) return 0;
  return log2(normCoarse/normFine);
}

static PetscErrorCode MGSetProfiling(MG mg,PetscBool flg) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (flg) {
    PetscInt lev;
    ierr = DMFEGetInfo(mg->dm,NULL,&lev,NULL,NULL,NULL);CHKERRQ(ierr);
    for ( ; lev >= 0; mg=mg->coarse,lev--) {
      char name[256];
      ierr = PetscSNPrintf(name,sizeof name,"FMG Level %D",lev);CHKERRQ(ierr);
      ierr = PetscLogStageRegister(name,&mg->stage);CHKERRQ(ierr);
      ierr = PetscSNPrintf(name,sizeof name,"MGVcycleLevel%D",lev);CHKERRQ(ierr);
      ierr = PetscLogEventRegister(name,DM_CLASSID,&mg->V_Cycle);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode MGRecordDiagnostics(Op op,MG mg,Vec B,Vec U) {
  PetscErrorCode ierr;
  Vec Y;

  PetscFunctionBegin;
  if (mg->monitor || mg->diagnostics) {
    ierr = DMGetGlobalVector(mg->dm,&Y);CHKERRQ(ierr);
    ierr = OpApply(op,mg->dm,U,Y);CHKERRQ(ierr);
    ierr = VecAYPX(Y,-1.,B);CHKERRQ(ierr);
    ierr = VecNorm(Y,NORM_2,&mg->rnorm2);CHKERRQ(ierr);
    if (mg->bnorm2 > 1e3*PETSC_MACHINE_EPSILON && mg->rnorm2 > 1e3*PETSC_MACHINE_EPSILON) mg->rnorm2 /= mg->bnorm2;
    ierr = DMRestoreGlobalVector(mg->dm,&Y);CHKERRQ(ierr);

    ierr = OpIntegrateNorms(op,mg->dm,U,&mg->enormInfty,&mg->enormL2);CHKERRQ(ierr);
  }
  if (mg->monitor) {
    PetscInt fedegree,level,mlocal[3],Mglobal[3],procs[3];
    PetscReal erateInfty = 0,erateL2 = 0,rrate2 = 0;
    ierr = DMFEGetInfo(mg->dm,&fedegree,&level,mlocal,Mglobal,procs);CHKERRQ(ierr);
    if (mg->coarse) {
      erateInfty = ConvergenceRate(mg->coarse->enormInfty,mg->enormInfty);
      erateL2    = ConvergenceRate(mg->coarse->enormL2,mg->enormL2);
      rrate2     = ConvergenceRate(mg->coarse->rnorm2,mg->rnorm2);
    }
    ierr = PetscPrintf(PetscObjectComm((PetscObject)mg->dm),"Q%D %2D e_max %8.2e(%3.1f) e_L2 %8.2e(%3.1f) r_2 %8.2e(%3.1f) G[%5D%5D%5D] L[%4D%4D%4D] P[%3D%3D%3D]\n",
                       fedegree,level,
                       mg->enormInfty,erateInfty,    // normalized by exact solution
                       mg->enormL2,erateL2,          // normalized by exact solution
                       mg->rnorm2,rrate2,            // algebraic, normalized by (non-tau-corrected) algebraic forcing
                       Mglobal[0],Mglobal[1],Mglobal[2],
                       mlocal[0],mlocal[1],mlocal[2],
                       procs[0],procs[1],procs[2]);CHKERRQ(ierr);
  }
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
  ierr = PetscLogEventBegin(mg->V_Cycle,mg->dm,B,U,0);CHKERRQ(ierr);
  if (presmooths) {
    ierr = KSPSetTolerances(mg->ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,mg->coarse?presmooths-1:20);CHKERRQ(ierr);
    ierr = KSPSolve(mg->ksp,B,U);CHKERRQ(ierr);
  }
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

  if (postsmooths) {
    ierr = KSPSetTolerances(mg->ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,postsmooths-1);CHKERRQ(ierr);
    ierr = KSPSolve(mg->ksp,B,U);CHKERRQ(ierr);
  }
  ierr = PetscLogEventEnd(mg->V_Cycle,mg->dm,B,U,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode MGFCycle(Op op,MG mg,PetscInt presmooths,PetscInt postsmooths,Vec B,Vec U) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscLogStagePush(mg->stage);CHKERRQ(ierr);
  ierr = VecNorm(B,NORM_2,&mg->bnorm2);CHKERRQ(ierr);
  if (mg->coarse) {
    Vec Uc,Bc;
    DM dmcoarse = mg->coarse->dm;
    if (dmcoarse) {
      ierr = DMGetGlobalVector(dmcoarse,&Bc);CHKERRQ(ierr);
      ierr = VecZeroEntries(Bc);CHKERRQ(ierr);
    } else Bc = NULL;
    ierr = OpRestrictResidual(op,mg->dm,B,Bc);CHKERRQ(ierr);
    if (dmcoarse) {
      ierr = DMFEZeroBoundaries(dmcoarse,Bc);CHKERRQ(ierr);
      ierr = DMGetGlobalVector(dmcoarse,&Uc);CHKERRQ(ierr);
      ierr = MGFCycle(op,mg->coarse,presmooths,postsmooths,Bc,Uc);CHKERRQ(ierr);
    } else Uc = NULL;
    ierr = OpInterpolate(op,mg->dm,Uc,U);CHKERRQ(ierr);
    if (dmcoarse) {
      ierr = DMRestoreGlobalVector(dmcoarse,&Uc);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(dmcoarse,&Bc);CHKERRQ(ierr);
    }
  }
  ierr = MGVCycle(op,mg,presmooths,postsmooths,B,U);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);
  ierr = MGRecordDiagnostics(op,mg,B,U);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RunMGV()
{
  PetscErrorCode ierr;
  Grid grid;
  Options opt;
  Op op;
  PetscInt fedegree,dof,addquadpts,nlevels;
  DM dm;
  Vec U0,U,F;
  MG mg;
  PetscReal normU0,normError,normF,normResid;

  PetscFunctionBegin;
  ierr = OpCreateFromOptions(PETSC_COMM_WORLD,&op);CHKERRQ(ierr);
  ierr = OpGetFEDegree(op,&fedegree);CHKERRQ(ierr);
  ierr = OpGetDof(op,&dof);CHKERRQ(ierr);
  ierr = OpGetAddQuadPts(op,&addquadpts);CHKERRQ(ierr);
  ierr = OptionsParse("Finite Element FAS FMG solver",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = GridGetNumLevels(grid,&nlevels);CHKERRQ(ierr);

  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = DMCreateFE(grid,fedegree,dof,addquadpts,&dm);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);

  ierr = DMCreateGlobalVector(dm,&U0);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&U);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&F);CHKERRQ(ierr);
  ierr = OpForcing(op,dm,F);CHKERRQ(ierr);
  ierr = OpSolution(op,dm,U0);CHKERRQ(ierr);
  ierr = VecNorm(U0,NORM_2,&normU0);CHKERRQ(ierr);
  ierr = VecNorm(F,NORM_2,&normF);CHKERRQ(ierr);

  ierr = MGCreate(op,dm,nlevels,&mg);CHKERRQ(ierr);
  ierr = MGSetProfiling(mg,PETSC_TRUE);CHKERRQ(ierr);
  for (PetscInt i=0; i<5; i++) {
    Vec Y;
    ierr = MGVCycle(op,mg,opt->smooth[0],opt->smooth[1],F,U);CHKERRQ(ierr);
    ierr = DMGetGlobalVector(dm,&Y);CHKERRQ(ierr);
    ierr = VecCopy(U,Y);CHKERRQ(ierr);
    ierr = VecAXPY(Y,-1,U0);CHKERRQ(ierr);
    ierr = VecNorm(Y,NORM_2,&normError);CHKERRQ(ierr);
    ierr = OpApply(op,dm,U,Y);CHKERRQ(ierr);
    ierr = VecAYPX(Y,-1,F);CHKERRQ(ierr);
    ierr = VecNorm(Y,NORM_2,&normResid);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"V(%D,%D) %2D: |e|_2/|u|_2 %8.2e  |r|_2/|f|_2 %8.2e\n",opt->smooth[0],opt->smooth[1],i,(double)(normError/normU0),(double)(normResid/normF));CHKERRQ(ierr);
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

// Smooth 5% mesh distortion implemented by "swirling" the region inside |(x,y)|_2 < 1.
PetscErrorCode DMCoordDistort(DM dm,const PetscReal L[]) {
  PetscErrorCode ierr;
  PetscScalar *xx;
  PetscInt m;
  Vec X;

  PetscFunctionBegin;
  ierr = DMGetCoordinatesLocal(dm,&X);CHKERRQ(ierr);
  ierr = VecGetLocalSize(X,&m);CHKERRQ(ierr);
  ierr = VecGetArray(X,&xx);CHKERRQ(ierr);
  for (PetscInt i=0; i<m; i+=3) {
    PetscReal x = xx[i]/L[0],y = xx[i+1]/L[1],z = xx[i+2]/L[2],r2,theta,newx,newy;
    r2 = PetscMin(PetscSqrtReal(PetscSqr(2*x-1) + PetscSqr(2*y-1)),1);
    theta = 0.1*PetscSqr(PetscCosReal(PETSC_PI*r2/2)) * PetscSinReal(PETSC_PI*z);
    newx = (.5 + PetscCosReal(theta)*(x-.5) - PetscSinReal(theta)*(y-.5))*L[0];
    newy = (.5 + PetscSinReal(theta)*(x-.5) + PetscCosReal(theta)*(y-.5))*L[1];
    if (0) printf("%6.3f %6.3f %6.3f:  r2 %6.3f  theta %6.3f  dx %6.3f %6.3f\n",xx[i+0],xx[i+1],xx[i+2],r2,theta,newx-xx[i],newy-xx[i+1]);
    xx[i+0] = newx;
    xx[i+1] = newy;
  }
  ierr = VecRestoreArray(X,&xx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RunFMG()
{
  PetscErrorCode ierr;
  Grid grid;
  Options opt;
  Op op;
  PetscInt fedegree,dof,addquadpts,nlevels;
  DM dm;
  Vec U0,U,F;
  MG mg;
  PetscReal normU0,normError,normF,normResid;

  PetscFunctionBegin;
  ierr = OpCreateFromOptions(PETSC_COMM_WORLD,&op);CHKERRQ(ierr);
  ierr = OpGetFEDegree(op,&fedegree);CHKERRQ(ierr);
  ierr = OpGetDof(op,&dof);CHKERRQ(ierr);
  ierr = OpGetAddQuadPts(op,&addquadpts);CHKERRQ(ierr);
  ierr = OptionsParse("Finite Element FAS FMG solver",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = GridGetNumLevels(grid,&nlevels);CHKERRQ(ierr);

  ierr = DMCreateFE(grid,fedegree,dof,addquadpts,&dm);CHKERRQ(ierr);
  ierr = DMFESetUniformCoordinates(dm,opt->L);CHKERRQ(ierr);
  if (opt->coord_distort) {ierr = DMCoordDistort(dm,opt->L);CHKERRQ(ierr);}

  ierr = DMCreateGlobalVector(dm,&U0);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&U);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&F);CHKERRQ(ierr);
  ierr = OpForcing(op,dm,F);CHKERRQ(ierr);
  ierr = OpSolution(op,dm,U0);CHKERRQ(ierr);
  ierr = VecNorm(U0,NORM_2,&normU0);CHKERRQ(ierr);
  ierr = VecNorm(F,NORM_2,&normF);CHKERRQ(ierr);

  ierr = MGCreate(op,dm,nlevels,&mg);CHKERRQ(ierr);
  ierr = MGSetProfiling(mg,PETSC_TRUE);CHKERRQ(ierr);
  for (PetscInt i=0; i<3; i++) {
    Vec Y;
    if (!i) {ierr = MGFCycle(op,mg,opt->smooth[0],opt->smooth[1],F,U);CHKERRQ(ierr);}
    else    {ierr = MGVCycle(op,mg,opt->smooth[0],opt->smooth[1],F,U);CHKERRQ(ierr);}
    ierr = DMGetGlobalVector(dm,&Y);CHKERRQ(ierr);
    ierr = VecCopy(U,Y);CHKERRQ(ierr);
    ierr = VecAXPY(Y,-1,U0);CHKERRQ(ierr);
    ierr = VecNorm(Y,NORM_2,&normError);CHKERRQ(ierr);
    ierr = OpApply(op,dm,U,Y);CHKERRQ(ierr);
    ierr = VecAYPX(Y,-1,F);CHKERRQ(ierr);
    ierr = VecNorm(Y,NORM_2,&normResid);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"%s(%D,%D) %2D: |e|_2/|u|_2 %8.2e  |r|_2/|f|_2 %8.2e\n",i?"V":"F",opt->smooth[0],opt->smooth[1],i,(double)(normError/normU0),(double)(normResid/normF));CHKERRQ(ierr);
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
