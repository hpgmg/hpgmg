#ifndef _fefas_h
#define _fefas_h

#include <petscdm.h>

typedef struct Grid_private *Grid;

PetscErrorCode GridCreate(MPI_Comm comm,const PetscInt M[3],const PetscInt p[3],const PetscInt pw[3],PetscInt cmax,Grid *grid);
PetscErrorCode GridDestroy(Grid *grid);
PetscErrorCode GridView(Grid grid);
PetscErrorCode DMCreateFE(Grid grid,PetscInt fedegree,PetscInt dof,DM *dmfe);
PetscErrorCode DMDestroyFE(DM *dm);
PetscErrorCode DMFESetUniformCoordinates(DM dm,const PetscReal L[]);
PetscErrorCode DMFEGetTensorEval(DM dm,PetscInt *P,PetscInt *Q,const PetscReal **B,const PetscReal **D,const PetscReal **x,PetscReal **w);
PetscErrorCode DMFEGetNumElements(DM dm,PetscInt *nelems);
PetscErrorCode DMFEExtractElements(DM dm,const PetscScalar *u,PetscInt elem,PetscInt ne,PetscScalar *y);
PetscErrorCode DMFESpaceSetElements(DM dm,PetscScalar *u,PetscInt elem,PetscInt ne,InsertMode imode,const PetscScalar *y);

#endif
