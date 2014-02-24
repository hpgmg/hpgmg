#ifndef _fefas_h
#define _fefas_h

#include <petscmat.h>

typedef struct Grid_private *Grid;

PetscErrorCode GridCreate(MPI_Comm comm,PetscInt m,PetscInt n,PetscInt p,PetscInt refinelev,Grid *grid);
PetscErrorCode GridDestroy(Grid *grid);
PetscErrorCode GridView(Grid grid);

#endif
