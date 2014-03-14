static const char help[] = "Geometric multigrid solver for finite-element elasticity.\n\n";

#include "fefas.h"

typedef struct Options_private *Options;
struct Options_private {
  char command[256];
  PetscInt M[3];
  PetscInt p[3];
  PetscInt cmax;
};

static PetscErrorCode TestGrid(Options opt)
{
  PetscErrorCode ierr;
  Grid grid;

  PetscFunctionBegin;
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,NULL,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode ActionParse(int argc,char *argv[],PetscErrorCode (**action)(Options))
{
  PetscFunctionList actionlist = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *action = NULL;

  ierr = PetscFunctionListAdd(&actionlist,"test-grid",TestGrid);CHKERRQ(ierr);

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

static PetscErrorCode OptionsParse(Options *opt)
{
  PetscErrorCode ierr;
  Options o;
  PetscInt three;

  PetscFunctionBegin;
  ierr = PetscNew(&o);CHKERRQ(ierr);
  o->M[0] = 10;
  o->M[1] = 10;
  o->M[2] = 10;
  o->p[0] = 1;
  o->p[1] = 1;
  o->p[2] = 1;
  o->cmax = 100;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Finite Element FAS solver",NULL);CHKERRQ(ierr);
  three = 3;
  ierr = PetscOptionsIntArray("-M","Fine grid dimensions","",o->M,&three,NULL);CHKERRQ(ierr);
  three = 3;
  ierr = PetscOptionsIntArray("-p","Process grid dimensions","",o->p,&three,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-cmax","Max coarse grid size","",o->cmax,&o->cmax,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  *opt = o;
  PetscFunctionReturn(0);
}

int main(int argc, char *argv[])
{
  PetscErrorCode ierr,(*action)(Options);
  Options opt;

  PetscInitialize(&argc,&argv,NULL,help);
  ierr = ActionParse(argc,argv,&action);CHKERRQ(ierr);
  if (!action) {
    PetscFinalize();
    return 1;
  }
  ierr = OptionsParse(&opt);CHKERRQ(ierr);
  ierr = (*action)(opt);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFinalize();
  return 0;
}
