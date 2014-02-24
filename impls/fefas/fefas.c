static const char help[] = "Geometric multigrid solver for finite-element elasticity.\n\n";

#include "fefas.h"

int main(int argc, char *argv[])
{
  PetscErrorCode ierr;
  Grid grid;
  PetscInt M[3] = {10,10,10},three = 3,L=1;

  PetscInitialize(&argc,&argv,NULL,help);
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,"Finite Element FAS solver",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsIntArray("-M","Coarse grid dimensions","",M,&three,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-L","Number of level refinements","",L,&L,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  ierr = GridCreate(PETSC_COMM_WORLD,M[0],M[1],M[2],L,&grid);CHKERRQ(ierr);
  ierr = GridView(grid);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  PetscFinalize();
  return 0;
}
