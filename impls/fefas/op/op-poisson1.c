#include "fefas-op.h"

PetscErrorCode OpCreate_Poisson1(Op op)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Create poisson1\n");CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
