static const char help[] = "Geometric multigrid solver for finite-element elasticity.\n\n";

#include "fefas.h"

PetscErrorCode TestGrid(void);
PetscErrorCode TestFESpace(void);
PetscErrorCode TestFEGrad(void);
PetscErrorCode TestFEInject(void);
PetscErrorCode TestFEInterp(void);
PetscErrorCode TestFERestrict(void);
PetscErrorCode TestOpApply(void);
PetscErrorCode TestOpDiagonal(void);
PetscErrorCode TestKSPSolve(void);
PetscErrorCode TestSampler(void);
PetscErrorCode RunMGV(void);
PetscErrorCode RunFMG(void);
PetscErrorCode RunSample(void);

static PetscErrorCode ActionParse(int argc,char *argv[],PetscErrorCode (**action)(void))
{
  PetscFunctionList actionlist = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *action = NULL;

  ierr = PetscFunctionListAdd(&actionlist,"test-grid",TestGrid);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-fespace",TestFESpace);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-fegrad",TestFEGrad);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-feinject",TestFEInject);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-feinterp",TestFEInterp);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-ferestrict",TestFERestrict);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-opapply",TestOpApply);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-opdiagonal",TestOpDiagonal);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-kspsolve",TestKSPSolve);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-sampler",TestSampler);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"mgv",RunMGV);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"fmg",RunFMG);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"sample",RunSample);CHKERRQ(ierr);

  if (argc < 2 || !argv[1] || argv[1][0] == '-') {
    ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDERR_WORLD,"First argument '%s' must be an action:",argc>=2&&argv[1]?argv[1]:"");CHKERRQ(ierr);
    ierr = PetscFunctionListView(actionlist,PETSC_VIEWER_STDERR_WORLD);CHKERRQ(ierr);
    goto out;
  }
  ierr = PetscFunctionListFind(actionlist,argv[1],action);CHKERRQ(ierr);
  if (!*action) {
    ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDERR_WORLD,"Unknown action '%s':",argc>=2&&argv[1]?argv[1]:"");CHKERRQ(ierr);
    ierr = PetscFunctionListView(actionlist,PETSC_VIEWER_STDERR_WORLD);CHKERRQ(ierr);
    goto out;
  }
  out:
  ierr = PetscFunctionListDestroy(&actionlist);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[])
{
  PetscErrorCode ierr,(*action)(void);

  PetscInitialize(&argc,&argv,NULL,help);
  ierr = ActionParse(argc,argv,&action);CHKERRQ(ierr);
  if (!action) {
    PetscFinalize();
    return 1;
  }
  ierr = (*action)();CHKERRQ(ierr);
  PetscFinalize();
  return 0;
}
