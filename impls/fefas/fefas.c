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

static PetscErrorCode TestFESpace(Options opt)
{
  PetscErrorCode ierr;
  Grid grid;
  DM dm;
  Vec G,L;
  PetscInt i,rstart,rend;
  PetscScalar *g;

  PetscFunctionBegin;
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,NULL,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = DMCreateFESpace(grid,1,1,&dm);CHKERRQ(ierr);
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
  ierr = DMDestroyFESpace(&dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode ActionParse(int argc,char *argv[],PetscErrorCode (**action)(Options))
{
  PetscFunctionList actionlist = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  *action = NULL;

  ierr = PetscFunctionListAdd(&actionlist,"test-grid",TestGrid);CHKERRQ(ierr);
  ierr = PetscFunctionListAdd(&actionlist,"test-fespace",TestFESpace);CHKERRQ(ierr);

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
