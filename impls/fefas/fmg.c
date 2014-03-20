#include "fefas.h"
#include "op/fefas-op.h"

typedef struct Options_private *Options;
struct Options_private {
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

#if 0
static PetscErrorCode VCycle(Op op,DM dm,Vec X)
{
  PetscErrorCode ierr;
  DM dmcoarse;

  PetscFunctionBegin;
  ierr = DMFECoarsen(dm,&dmcoarse);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

PetscErrorCode FMG()
{
  PetscErrorCode ierr;
  Grid grid;
  Options opt;
  Op op;

  PetscFunctionBegin;
  ierr = OpCreateFromOptions(PETSC_COMM_WORLD,&op);CHKERRQ(ierr);
  ierr = OptionsParse("Finite Element FAS FMG solver",&opt);CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,opt->M,opt->p,NULL,opt->cmax,&grid);CHKERRQ(ierr);
  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  ierr = OpDestroy(&op);CHKERRQ(ierr);
  ierr = PetscFree(opt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
