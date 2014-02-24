#include "fefas.h"

struct Grid_private {
  MPI_Comm comm;
  Grid coarse;
  PetscInt refinelev;
  PetscInt M[3];
  PetscInt s[3],m[3];
  PetscMPIInt neighborranks[3][3][3];
};

PetscErrorCode GridCreate(MPI_Comm comm,PetscInt m,PetscInt n,PetscInt p,PetscInt refinelev,Grid *grid)
{
  PetscErrorCode ierr;
  Grid g;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,&g);CHKERRQ(ierr);
  ierr = MPI_Comm_dup(comm,&g->comm);CHKERRQ(ierr);
  g->refinelev = refinelev;
  if (refinelev) {
    PetscMPIInt size;
    ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
    if (1 << (3*refinelev) == size) {
      /* Redundant coarse grid
       * We order ranks hierarchically, similar to cache-oblivious orderings.
       */
      PetscMPIInt rank,crank,color,key,i;
      MPI_Comm redcomm;

      ierr = MPI_Comm_rank(g->comm,&rank);CHKERRQ(ierr);
      color = rank%8;
      key = rank/8;
      ierr = MPI_Comm_split(g->comm,color,key,&redcomm);CHKERRQ(ierr);
      ierr = GridCreate(redcomm,m,n,p,refinelev-1,&g->coarse);CHKERRQ(ierr);
      ierr = MPI_Comm_rank(g->coarse->comm,&crank);CHKERRQ(ierr);
      for (i=0; i<3; i++) {
        g->M[i] = 2*g->coarse->M[i];
        g->m[i] = g->coarse->m[i];
        g->s[i] = 2*g->coarse->s[i] + (color & (1 << i) ? g->m[i] : 0);
      }
      ierr = MPI_Comm_free(&redcomm);CHKERRQ(ierr);
    } else {
      /* Coarse grid is on same process set
       * All coarsening can be done in-place.
       */
      PetscInt i;

      if ((1 << 3*refinelev) % size) SETERRQ4(comm,PETSC_ERR_ARG_INCOMP,"Communicator size %d incompatible with refinelevel %D, should divide 8^%D=%D",size,refinelev,refinelev,1 << (3*refinelev));
      ierr = GridCreate(g->comm,m,n,p,refinelev-1,&g->coarse);CHKERRQ(ierr);
      for (i=0; i<3; i++) {
        g->M[i] = 2*g->coarse->M[i];
        g->m[i] = 2*g->coarse->m[i];
        g->s[i] = 2*g->coarse->s[i];
      }
      ierr = PetscMemcpy(g->neighborranks[0][0],g->coarse->neighborranks[0][0],sizeof g->neighborranks);CHKERRQ(ierr);
    }
  } else {
    g->coarse = NULL;
    g->m[0] = g->M[0] = m;
    g->m[1] = g->M[1] = n;
    g->m[2] = g->M[2] = p;
    g->s[0] = g->s[1] = g->s[2] = 0;
    ierr = PetscMemzero(g->neighborranks[0][0],sizeof g->neighborranks);CHKERRQ(ierr);
  }

  *grid = g;
  PetscFunctionReturn(0);
}

PetscErrorCode GridDestroy(Grid *grid)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*grid) PetscFunctionReturn(0);
  ierr = GridDestroy(&(*grid)->coarse);CHKERRQ(ierr);
  ierr = MPI_Comm_free(&(*grid)->comm);CHKERRQ(ierr);
  ierr = PetscFree(*grid);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode GridView(Grid grid)
{
  PetscErrorCode ierr;
  Grid g;
  PetscMPIInt rank;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(grid->comm,&rank);CHKERRQ(ierr);
  for (g=grid; g; g=g->coarse) {
    ierr = PetscSynchronizedPrintf(grid->comm,"[%d] Level %D: [%4D:%4D,%4D:%4D,%4D:%4D] of [%4D,%4D,%4D]\n",rank,g->refinelev,
                                   g->s[0],g->s[0]+g->m[0],
                                   g->s[1],g->s[1]+g->m[1],
                                   g->s[2],g->s[2]+g->m[2],
                                   g->M[0],g->M[1],g->M[2]);CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(grid->comm,PETSC_STDOUT);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
