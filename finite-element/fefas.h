#ifndef _fefas_h
#define _fefas_h

#include <petscdm.h>
#include <stdint.h>
#include "op/fefas-op.h"

typedef struct Grid_private *Grid;

typedef enum {DOMAIN_INTERIOR=0x1,DOMAIN_EXTERIOR=0x2,DOMAIN_CLOSURE=0x3} DomainMode;

PetscErrorCode GridCreate(MPI_Comm comm,const PetscInt M[3],const PetscInt p[3],PetscInt cmax,Grid *grid);
PetscErrorCode GridDestroy(Grid *grid);
PetscErrorCode GridView(Grid grid);
PetscErrorCode GridGetNumLevels(Grid grid,PetscInt *nlevels);
PetscInt GridLevelFromM(const PetscInt M[3]);
PetscErrorCode DMCreateFE(Grid grid,PetscInt fedegree,PetscInt dof,DM *dmfe);
PetscErrorCode DMDestroyFE(DM *dm);
PetscErrorCode DMFESetUniformCoordinates(DM dm,const PetscReal L[]);
PetscErrorCode DMFEGetUniformCoordinates(DM dm,PetscReal L[]);
PetscErrorCode DMFEGetInfo(DM dm,PetscInt *fedegree,PetscInt *level,PetscInt mlocal[],PetscInt Mglobal[],PetscInt procs[]);
PetscErrorCode DMFEGetTensorEval(DM dm,PetscInt *P,PetscInt *Q,const PetscReal **B,const PetscReal **D,const PetscReal **x,const PetscReal **w,const PetscReal **w3);
PetscErrorCode DMFEGetNumElements(DM dm,PetscInt *nelems);
PetscErrorCode DMFEExtractElements(DM dm,const PetscScalar *u,PetscInt elem,PetscInt ne,PetscScalar *y);
PetscErrorCode DMFESetElements(DM dm,PetscScalar *u,PetscInt elem,PetscInt ne,InsertMode imode,DomainMode dmode,const PetscScalar *y);
PetscErrorCode DMFECoarsen(DM dm,DM *dmcoarse);
PetscErrorCode DMFEInject(DM dm,Vec Uf,Vec Uc);
PetscErrorCode DMFEInterpolate(DM dm,Vec Uc,Vec Uf);
PetscErrorCode DMFERestrict(DM dm,Vec Uf,Vec Uc);
PetscErrorCode DMFEZeroBoundaries(DM dm,Vec U);
PetscErrorCode DMCoordDistort(DM dm,const PetscReal L[]);

typedef struct MG_private *MG;
PetscErrorCode MGCreate(Op op,DM dm,PetscInt nlevels,MG *newmg);
PetscErrorCode MGDestroy(MG *mg);
PetscErrorCode MGMonitorSet(MG mg,PetscBool mon);
PetscErrorCode MGSetUpPC(MG mg);
PetscErrorCode MGFCycle(Op op,MG mg,PetscInt presmooths,PetscInt postsmooths,Vec B,Vec U);

PetscInt SampleGridNumLevels(const PetscInt p[]);
int64_t SampleGridNumElements(const PetscInt p[]);
PetscErrorCode SampleGridRangeCreate(PetscMPIInt nranks,PetscInt minlocal,PetscInt maxlocal,PetscInt maxsamples,PetscInt *nsamples,PetscInt **gridsizes);
PetscErrorCode ProcessGridFindSquarest(PetscMPIInt nranks,PetscInt squarest[3]);

PetscErrorCode MemoryGetUsage(double *heapused,double *heapavail);

#endif
