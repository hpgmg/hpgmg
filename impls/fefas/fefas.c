static const char help[] = "Geometric multigrid solver for finite-element elasticity.\n\n";

#include "fefas.h"

int main(int argc, char *argv[])
{
  PetscErrorCode ierr;
  Grid grid;
  PetscInt M[3] = {10,10,10},p[3] = {1,1,1},three = 3,cmax=15;

  PetscInitialize(&argc,&argv,NULL,help);
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Finite Element FAS solver",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsIntArray("-M","Fine grid dimensions","",M,&three,NULL);CHKERRQ(ierr);
  three = 3;
  ierr = PetscOptionsIntArray("-p","Process grid dimensions","",p,&three,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-cmax","Max coarse grid size","",cmax,&cmax,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,M,p,NULL,cmax,&grid);CHKERRQ(ierr);
  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  PetscFinalize();
  return 0;
}
