#ifndef _fefas_h
#define _fefas_h

#include <petscmat.h>

typedef struct Grid_private *Grid;

PetscErrorCode GridCreate(MPI_Comm comm,const PetscInt M[3],const PetscInt p[3],const PetscInt pw[3],PetscInt cmax,Grid *grid);
PetscErrorCode GridDestroy(Grid *grid);
PetscErrorCode GridView(Grid grid);

#endif
