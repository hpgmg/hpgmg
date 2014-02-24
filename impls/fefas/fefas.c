static const char help[] = "Geometric multigrid solver for finite-element elasticity.\n\n";

#include <petscksp.h>

int main(int argc, char *argv[])
{

  PetscInitialize(&argc,&argv,NULL,help);
  PetscPrintf(PETSC_COMM_WORLD,"Not much here yet.\n");
  PetscFinalize();
  return 0;
}
