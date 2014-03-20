#include "fefas-op.h"

PetscErrorCode OpRegisterAll_Generated(void);

static PetscFunctionList OpList;
static PetscBool OpPackageInitialized;

struct Op_private {
  MPI_Comm comm;                /* Finest level comm (only for diagnostics at setup time) */
  PetscErrorCode (*Destroy)(Op);
  void *data;
};

PetscErrorCode OpCreateFromOptions(MPI_Comm comm,Op *op)
{
  PetscErrorCode ierr,(*f)(Op);
  Op o;
  char opname[256] = "poisson1";

  PetscFunctionBegin;
  ierr = OpInitializePackage();CHKERRQ(ierr);
  ierr = PetscNew(&o);CHKERRQ(ierr);
  ierr = PetscCommDuplicate(comm,&o->comm,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBegin(comm,NULL,"Operator options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsFList("-op_type","Operator type","",OpList,opname,opname,sizeof opname,NULL);CHKERRQ(ierr);
  ierr = PetscFunctionListFind(OpList,opname,&f);CHKERRQ(ierr);
  ierr = (*f)(o);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);
  *op = o;
  PetscFunctionReturn(0);
}

PetscErrorCode OpDestroy(Op *op)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*op) PetscFunctionReturn(0);
  if ((*op)->Destroy) {
    ierr = (*op)->Destroy(*op);CHKERRQ(ierr);
  }
  ierr = PetscCommDestroy(&(*op)->comm);CHKERRQ(ierr);
  ierr = PetscFree(*op);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpRegister(const char *name,PetscErrorCode (*f)(Op))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListAdd(&OpList,name,f);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode OpFinalizePackage()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFunctionListDestroy(&OpList);CHKERRQ(ierr);
  OpPackageInitialized = PETSC_FALSE;
  PetscFunctionReturn(0);
}

PetscErrorCode OpInitializePackage()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (OpPackageInitialized) PetscFunctionReturn(0);
  ierr = OpRegisterAll_Generated();CHKERRQ(ierr);
  ierr = PetscRegisterFinalize(OpFinalizePackage);CHKERRQ(ierr);
  OpPackageInitialized = PETSC_TRUE;
  PetscFunctionReturn(0);
}
