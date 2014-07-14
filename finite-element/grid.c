#include "fefas.h"
#include <petscsf.h>
#include <petscdmshell.h>
#include <petscdt.h>
#include <stdint.h>

typedef uint64_t zcode;

// Coarse grids should never be smaller than fine grids, so 27 neighbors should always be enough.
#define COARSE_KNOWN_MAX 27

struct Grid_private {
  PetscInt refct;
  MPI_Comm comm;
  Grid coarse;
  PetscInt level;
  PetscInt M[3],p[3];
  PetscInt s[3],m[3];   // Owned cells of global grid
  struct {
    PetscMPIInt rank;
    PetscInt ri[3];
    PetscInt s[3],m[3];
  } coarseknown[COARSE_KNOWN_MAX];
  PetscInt ncoarseknown;
  PetscMPIInt neighborranks[3][3][3];
};

// This type is hidden behind DM (not exposed publicly)
typedef struct FE_private *FE;
struct FE_private {
  Grid grid;
  DM dmcoarse;
  PetscInt degree;      // Finite element polynomial degree
  PetscInt dof;         // Number of degrees of freedom per vertex
  PetscInt om[3];       // Array dimensions of owned part of global vectors
  // Local vectors have sufficient fringe to include C-points needed for interpolation and (soon) overlap for segmental
  // refinement.  Vertices 0 and 0+lM-1 will always correspond to a vertex at the corner of a coarse elements.  In the
  // Q2 example below:
  //
  //   C is a coarse corner point
  //   F is a fine corner point and coarse interior point
  //   f is a fine interior point
  //
  // Suppose that the fine-grid partition had the current process owning vertices [2,a), thus the active region (with a
  // non-overlapping element partition) is [2,a].
  //
  // 0C -- 1f -- 2F -- 3f -- 4C -- 5f -- 6F -- 7f -- 8C -- 9f -- aF -- bf -- cC
  //             ^                                           ^     ^           ^
  //           ls=2                                        om=lm-1           lM=c+1=13
  //                                                             lm=a+1-2=9
  PetscInt lM[3];       // Array dimensions of local vectors (includes enough fringe for interpolation)
  PetscInt ls[3],lm[3]; // Start and extent of active part of local vector
  PetscReal Luniform[3];
  PetscBool hascoordinates;
  MPI_Datatype unit;
  PetscSF sf;
  PetscSF sfinject;
  PetscSF sfinjectLocal;
  struct {
    PetscReal *B;
    PetscReal *D;
    PetscReal *x;
    PetscReal *w;
    PetscReal *interp;
    PetscReal *w3;
  } ref;
};

static PetscInt Idx3(const PetscInt m[],PetscInt i,PetscInt j,PetscInt k) { return (i*m[1]+j)*m[2]+k; }
static PetscInt FEIdxO(FE fe,PetscInt i,PetscInt j,PetscInt k) { return Idx3(fe->om,i,j,k); }
static PetscInt FEIdxL(FE fe,PetscInt i,PetscInt j,PetscInt k) { return Idx3(fe->lM,i,j,k); }
static PetscInt FEIdxLs(FE fe,PetscInt i,PetscInt j,PetscInt k) { return Idx3(fe->lM,fe->ls[0]+i,fe->ls[1]+j,fe->ls[2]+k); }

// Walk through the coarse grids we know about and find the (rank,offset) of local index i,j,k (must be a C-point)
static PetscErrorCode FEGridFindCoarseRankIndex(FE fe,PetscInt i,PetscInt j,PetscInt k,PetscInt *rank,PetscInt *index) {
  Grid grid = fe->grid;
  PetscInt g[3] = {(fe->grid->s[0]/2)*2*fe->degree+fe->ls[0]+i,
                   (fe->grid->s[1]/2)*2*fe->degree+fe->ls[1]+j,
                   (fe->grid->s[2]/2)*2*fe->degree+fe->ls[2]+k}; // Global index in FE node space

  PetscFunctionBegin;
  for (PetscInt l=0; l<3; l++) if (g[l] % 2) SETERRQ6(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Expected C-point %D %D %D (global %D %D %D)",i,j,k,g[0],g[1],g[2]);
  for (PetscInt l=0; l<3; l++) g[l] /= 2; // Now a global coarse FE node location
  for (PetscInt l=0; l<grid->ncoarseknown; l++) { // Try each coarse grid that we know about
    PetscInt Cos[3],Com[3];
    for (PetscInt d=0; d<3; d++) { // Determine ownership in finite-element space
      PetscInt s = grid->coarseknown[l].s[d],m = grid->coarseknown[l].m[d];
      Cos[d] = s*fe->degree;
      Com[d] = m*fe->degree + (s+m == fe->grid->M[d]/2);
    }
    if (   Cos[0] <= g[0] && g[0] < Cos[0]+Com[0]
        && Cos[1] <= g[1] && g[1] < Cos[1]+Com[1]
        && Cos[2] <= g[2] && g[2] < Cos[2]+Com[2]) {
      *rank = grid->coarseknown[l].rank;
      *index = Idx3(Com,g[0]-Cos[0],g[1]-Cos[1],g[2]-Cos[2]);
      PetscFunctionReturn(0);
    }
  }
  SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Could not find owner of %D,%D,%D",g[0],g[1],g[2]);
  PetscFunctionReturn(0);
}

static PetscInt CeilDiv(PetscInt a,PetscInt b) {return a/b + !!(a%b);}

static unsigned ZCodeSplit1(zcode z) {
  zcode t;
  t = (z & 00101010101010101) | ((z>>2 ) & 00202020202020202);
  t = (t & 00003000300030003) | ((t>>4 ) & 00014001400140014);
  t = (t & 00000001700000017) | ((t>>8 ) & 00000036000000360);
  t = (t & 00000000000000377) | ((t>>16) & 00000000000177400);
  return (unsigned)t;
}

static void ZCodeSplit(zcode z,PetscMPIInt i[3]) {
  i[0] = ZCodeSplit1(z >> 2);
  i[1] = ZCodeSplit1(z >> 1);
  i[2] = ZCodeSplit1(z >> 0);
}

static zcode ZCodeFromRank(PetscMPIInt rank,const PetscInt p[3]) {
  zcode z;
  PetscMPIInt r;
  for (z=0,r=-1; r != rank; z++) {
    PetscMPIInt i[3];
    ZCodeSplit(z,i);
    if (i[0] < p[0] && i[1] < p[1] && i[2] < p[2]) r++;
  }
  return z-1;
}

PetscInt GridLevelFromM(const PetscInt M[3]) {
  PetscInt m[3] = {M[0],M[1],M[2]},level;
  for (level=0; m[0]%2==0 && m[1]%2==0 && m[2]%2==0; level++) {
    m[0] /= 2;
    m[1] /= 2;
    m[2] /= 2;
  }
  return level;
}

// The range {0..M-1} is partitioned into p parts.  Find which part i falls in.
static PetscInt PartitionFind(PetscInt M,PetscInt p,PetscInt i) {
  PetscInt t = i/(M/p);  // M/p is a lower bound for actual subdomain size, so t is an upper bound
  while ((M/p)*t + PetscMin(M%p,t) > i) t--;
  return t;
}

// The range {0..M-1} is partitioned into p parts.  Find the start and extent of part t.
static PetscErrorCode PartitionGetRange(PetscInt M,PetscInt p,PetscInt t,PetscInt *s,PetscInt *m) {
  *m = M/p + (M%p > t);
  *s = (M/p)*t + PetscMin(M%p,t);
  return 0;
}

PetscErrorCode GridCreate(MPI_Comm comm,const PetscInt M[3],const PetscInt p[3],PetscInt cmax,Grid *grid)
{
  PetscErrorCode ierr;
  Grid g;
  PetscMPIInt size,rank;
  PetscInt CM[3],Cp[3],j;
  zcode z;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,&g);CHKERRQ(ierr);
  ierr = PetscCommDuplicate(comm,&g->comm,NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (size != p[0]*p[1]*p[2]) SETERRQ4(comm,PETSC_ERR_ARG_INCOMP,"Communicator size %d incompatible with process grid %D,%D,%D",size,p[0],p[1],p[2]);

  for (j=0; j<3; j++) {
    Cp[j] = CeilDiv(p[j],2); // Proposed coarse process set if we restrict communicator
    CM[j] = M[j]/2;          // Coarse grid size (if there is a coarser grid)
    g->M[j] = M[j];
    g->p[j] = p[j];
  }
  g->level = GridLevelFromM(M);

  // Find ownership
  z = ZCodeFromRank(rank,p);

  if (M[0]%2 || M[1]%2 || M[2]%2) { // Already on the coarsest possible grid
    if (size != 1) SETERRQ7(comm,PETSC_ERR_ARG_INCOMP,"Grid %D,%D,%D reached on P[%D %D %D], try increasing cmax %D or global grid to coarsen processes sooner",M[0],M[1],M[2],p[0],p[1],p[2],cmax);
    // Coarsest grid, on a single process
    g->coarse = NULL;
    for (j=0; j<3; j++) {
      g->m[j] = M[j];
      g->s[j] = 0;
    }
    for (j=0; j<27; j++) g->neighborranks[j/9%3][j/3%3][j%3] = -1;
    g->ncoarseknown = 0;
  } else {
    PetscMPIInt ri[3],r;
    zcode t;
    PetscInt cnt,coarsecnt,nneighbors,Cm[3],Cs[3],mask;

    if (CeilDiv(CM[0],Cp[0]) * CeilDiv(CM[1],Cp[1]) * CeilDiv(CM[2],Cp[2]) > cmax || size == 1) {
      for (j=0; j<3; j++) Cp[j] = p[j]; // Coarsen on the same process set
      ierr = GridCreate(comm,CM,Cp,cmax,&g->coarse);CHKERRQ(ierr);
      mask = 00;
    } else { // z&07==0 will participate in coarse grid
      MPI_Comm ccomm;
      ierr = MPI_Comm_split(comm,z&07?MPI_UNDEFINED:0,0,&ccomm);CHKERRQ(ierr);
      if (ccomm != MPI_COMM_NULL) { // I continue to coarse grid
        ierr = GridCreate(ccomm,CM,Cp,cmax,&g->coarse);CHKERRQ(ierr);
        ierr = MPI_Comm_free(&ccomm);CHKERRQ(ierr);
      } else g->coarse = NULL;
      mask = 07;
    }

    ZCodeSplit(z,ri);
    for (j=0; j<3; j++) {
      // Determine the elements that I will own
      ierr = PartitionGetRange(M[j],p[j],ri[j],&g->s[j],&g->m[j]);CHKERRQ(ierr);
      // Range of coarse elements that cover my subdomain
      Cs[j] = g->s[j]/2;
      Cm[j] = PetscMin(CeilDiv(g->s[j]+g->m[j],2)+1,CM[j]) - Cs[j];
    }

    // Determine which coarse processes we need to know about
    g->ncoarseknown = 0;
    for (PetscInt i=Cs[0]; i<Cs[0]+Cm[0]; i++) {
      for (PetscInt j=Cs[1]; j<Cs[1]+Cm[1]; j++) {
        for (PetscInt k=Cs[2]; k<Cs[2]+Cm[2]; k++) {
          PetscInt pi = PartitionFind(CM[0],Cp[0],i);
          PetscInt pj = PartitionFind(CM[1],Cp[1],j);
          PetscInt pk = PartitionFind(CM[2],Cp[2],k);
          PetscInt l;
          for (l=0; l<g->ncoarseknown; l++) {
            if (g->coarseknown[l].ri[0] == pi && g->coarseknown[l].ri[1] == pj && g->coarseknown[l].ri[2] == pk) goto skip; // we already know about this point
          }
          if (l == COARSE_KNOWN_MAX) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Too many coarse known; probably algorithm bug");
          g->coarseknown[l].ri[0] = pi;
          g->coarseknown[l].ri[1] = pj;
          g->coarseknown[l].ri[2] = pk;
          ierr = PartitionGetRange(CM[0],Cp[0],pi,&g->coarseknown[l].s[0],&g->coarseknown[l].m[0]);CHKERRQ(ierr);
          ierr = PartitionGetRange(CM[1],Cp[1],pj,&g->coarseknown[l].s[1],&g->coarseknown[l].m[1]);CHKERRQ(ierr);
          ierr = PartitionGetRange(CM[2],Cp[2],pk,&g->coarseknown[l].s[2],&g->coarseknown[l].m[2]);CHKERRQ(ierr);
          g->coarseknown[l].rank = -1;
          g->ncoarseknown++;
          skip: continue;
        }
      }
    }

    for (cnt=0,nneighbors=0; cnt<27; cnt++) { // Count how many horizontal neighbors exist on grid
      PetscMPIInt i[3] = {ri[0]+cnt/9%3-1,ri[1]+cnt/3%3-1,ri[2]+cnt%3-1};
      if (0 <= i[0] && i[0] < p[0] && 0 <= i[1] && i[1] < p[1] && 0 <= i[2] && i[2] < p[2]) nneighbors++;
      g->neighborranks[cnt/9%3][cnt/3%3][cnt%3] = -1;
    }
    for (t=0,cnt=0,coarsecnt=0,r=0; cnt<nneighbors || coarsecnt < g->ncoarseknown; t++) {
      PetscMPIInt i[3];
      ZCodeSplit(t,i);
      if (i[0] < p[0] && i[1] < p[1] && i[2] < p[2]) {
        for (j=0; j<3; j++) i[j] -= ri[j];
        if (PetscAbs(i[0]) < 2 && PetscAbs(i[1]) < 2 && PetscAbs(i[2]) < 2) {
          g->neighborranks[i[0]+1][i[1]+1][i[2]+1] = r;
          cnt++;
        }
        if ((t&mask) == 0) { // zcode=t participates in the coarse grid; find out if we asked for it
          PetscMPIInt ci[3];
          ZCodeSplit(mask?t>>3:t,ci);
          for (PetscInt l=0; l<g->ncoarseknown; l++) {
            if (ci[0] == g->coarseknown[l].ri[0] && ci[1] == g->coarseknown[l].ri[1] && ci[2] == g->coarseknown[l].ri[2]) {
              g->coarseknown[l].rank = r;
              coarsecnt++;
              break;
            }
          }
        }
        r++;
      }
    }
  }

  if (g->m[0]*g->m[1]*g->m[2] <= 0) SETERRQ7(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Invalid grid L[%D,%D,%D] of G[%D,%D,%D] on level %D",g->m[0],g->m[1],g->m[2],g->M[0],g->M[1],g->M[2],g->level);

  g->refct = 1;
  *grid = g;
  PetscFunctionReturn(0);
}

PetscErrorCode GridDestroy(Grid *grid)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*grid) PetscFunctionReturn(0);
  if (--(*grid)->refct > 0) PetscFunctionReturn(0);
  ierr = GridDestroy(&(*grid)->coarse);CHKERRQ(ierr);
  ierr = PetscCommDestroy(&(*grid)->comm);CHKERRQ(ierr);
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
    PetscInt M[3] = {g->M[0],g->M[1],g->M[2]},j;
    for (j=0; M[0]%2==0 &&  M[1]%2==0 &&  M[2]%2==0; j++) { // Compute current level
      M[0] /= 2; M[1] /= 2; M[2] /= 2;
    }
    ierr = PetscSynchronizedPrintf(grid->comm,"[%d] Level %D: [%4D:%4D,%4D:%4D,%4D:%4D] of [%4D,%4D,%4D] on [%3D,%3D,%3D]\n",rank,j,
                                   g->s[0],g->s[0]+g->m[0],
                                   g->s[1],g->s[1]+g->m[1],
                                   g->s[2],g->s[2]+g->m[2],
                                   g->M[0],g->M[1],g->M[2],
                                   g->p[0],g->p[1],g->p[2]);CHKERRQ(ierr);
  }
  ierr = PetscSynchronizedFlush(grid->comm,PETSC_STDOUT);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode GridGetNumLevels(Grid grid,PetscInt *nlevels) {
  PetscErrorCode ierr;
  Grid g;

  PetscFunctionBegin;
  *nlevels = 0;
  for (g=grid; g; g=g->coarse,(*nlevels)++);
  ierr = MPI_Allreduce(MPI_IN_PLACE,nlevels,1,MPIU_INT,MPI_MAX,grid->comm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMCreateGlobalVector_FE(DM dm,Vec *G)
{
  PetscErrorCode ierr;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  ierr = VecCreate(PetscObjectComm((PetscObject)dm),G);CHKERRQ(ierr);
  ierr = VecSetBlockSize(*G,fe->dof);CHKERRQ(ierr);
  ierr = VecSetSizes(*G,fe->om[0]*fe->om[1]*fe->om[2]*fe->dof,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetUp(*G);CHKERRQ(ierr);
  ierr = VecSetDM(*G,dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMCreateLocalVector_FE(DM dm,Vec *G)
{
  PetscErrorCode ierr;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_SELF,G);CHKERRQ(ierr);
  ierr = VecSetBlockSize(*G,fe->dof);CHKERRQ(ierr);
  ierr = VecSetSizes(*G,fe->lM[0]*fe->lM[1]*fe->lM[2]*fe->dof,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetUp(*G);CHKERRQ(ierr);
  ierr = VecSetDM(*G,dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMGlobalToLocalBegin_FE(DM dm,Vec G,InsertMode imode,Vec L)
{
  PetscErrorCode ierr;
  FE fe;
  const PetscScalar *g;
  PetscScalar *l;

  PetscFunctionBegin;
  if (imode != INSERT_VALUES) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"InsertMode");
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  ierr = VecGetArrayRead(G,&g);CHKERRQ(ierr);
  ierr = VecGetArray(L,&l);CHKERRQ(ierr);
  ierr = PetscSFBcastBegin(fe->sf,fe->unit,g,l);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(G,&g);CHKERRQ(ierr);
  ierr = VecRestoreArray(L,&l);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMGlobalToLocalEnd_FE(DM dm,Vec G,InsertMode imode,Vec L)
{
  PetscErrorCode ierr;
  FE fe;
  PetscInt i,j,k,d;
  const PetscScalar *g;
  PetscScalar *l;

  PetscFunctionBegin;
  if (imode != INSERT_VALUES) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"InsertMode");
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  ierr = VecGetArrayRead(G,&g);CHKERRQ(ierr);
  ierr = VecGetArray(L,&l);CHKERRQ(ierr);

  // Copy over local part
  for (i=0; i<fe->om[0]; i++) {
    for (j=0; j<fe->om[1]; j++) {
      for (k=0; k<fe->om[2]; k++) {
        for (d=0; d<fe->dof; d++) {
          l[FEIdxLs(fe,i,j,k)*fe->dof+d] = g[FEIdxO(fe,i,j,k)*fe->dof+d];
        }
      }
    }
  }
  ierr = PetscSFBcastEnd(fe->sf,fe->unit,g,l);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(G,&g);CHKERRQ(ierr);
  ierr = VecRestoreArray(L,&l);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMLocalToGlobalBegin_FE(DM dm,Vec L,InsertMode imode,Vec G)
{
  PetscErrorCode ierr;
  FE fe;
  const PetscScalar *l;
  PetscScalar *g;
  MPI_Op op;

  PetscFunctionBegin;
  switch (imode) {
  case ADD_VALUES: op = MPIU_SUM; break;
  case INSERT_VALUES: PetscFunctionReturn(0); // Local copy will be done in DMLocalToGlobalEnd_FE
  default: SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"InsertMode");
  }
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  ierr = VecGetArrayRead(L,&l);CHKERRQ(ierr);
  ierr = VecGetArray(G,&g);CHKERRQ(ierr);
  ierr = PetscSFReduceBegin(fe->sf,fe->unit,l,g,op);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(L,&l);CHKERRQ(ierr);
  ierr = VecRestoreArray(G,&g);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMLocalToGlobalEnd_FE(DM dm,Vec L,InsertMode imode,Vec G)
{
  PetscErrorCode ierr;
  FE fe;
  PetscInt i,j,k,d;
  const PetscScalar *l;
  PetscScalar *g;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  ierr = VecGetArrayRead(L,&l);CHKERRQ(ierr);
  ierr = VecGetArray(G,&g);CHKERRQ(ierr);

  // Add local part
  for (i=0; i<fe->om[0]; i++) {
    for (j=0; j<fe->om[1]; j++) {
      for (k=0; k<fe->om[2]; k++) {
        for (d=0; d<fe->dof; d++) {
          PetscInt src = FEIdxLs(fe,i,j,k)*fe->dof+d;
          PetscInt dst = FEIdxO(fe,i,j,k)*fe->dof+d;
          if (imode == ADD_VALUES) g[dst] += l[src];
          else                     g[dst]  = l[src];
        }
      }
    }
  }
  if (imode == ADD_VALUES) {
    PetscLogFlops(fe->om[0]*fe->om[1]*fe->om[2]*fe->dof);
    ierr = PetscSFReduceEnd(fe->sf,fe->unit,l,g,MPIU_SUM);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(L,&l);CHKERRQ(ierr);
  ierr = VecRestoreArray(G,&g);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMFESetUniformCoordinates(DM dm,const PetscReal L[])
{
  PetscErrorCode ierr;
  DM dmc;
  FE fe;
  Vec X;
  PetscScalar *x;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  if (fe->degree > 2) SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"fe->degree %D > 2",fe->degree);
  ierr = DMCreateFE(fe->grid,fe->degree,3,&dmc);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(dmc,&X);CHKERRQ(ierr);
  ierr = VecGetArray(X,&x);CHKERRQ(ierr);
  for (PetscInt i=0; i<fe->lM[0]; i++) {
    for (PetscInt j=0; j<fe->lM[1]; j++) {
      for (PetscInt k=0; k<fe->lM[2]; k++) {
        x[FEIdxL(fe,i,j,k)*3+0] = L[0]*(fe->grid->s[0]/2*2 + (PetscReal)i/fe->degree)/fe->grid->M[0];
        x[FEIdxL(fe,i,j,k)*3+1] = L[1]*(fe->grid->s[1]/2*2 + (PetscReal)j/fe->degree)/fe->grid->M[1];
        x[FEIdxL(fe,i,j,k)*3+2] = L[2]*(fe->grid->s[2]/2*2 + (PetscReal)k/fe->degree)/fe->grid->M[2];
      }
    }
  }
  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = DMSetCoordinateDM(dm,dmc);CHKERRQ(ierr);
  ierr = DMSetCoordinatesLocal(dm,X);CHKERRQ(ierr);
  ierr = PetscMemcpy(fe->Luniform,L,sizeof fe->Luniform);CHKERRQ(ierr);
  fe->hascoordinates = PETSC_TRUE;
  ierr = DMDestroy(&dmc);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMFEGetUniformCoordinates(DM dm,PetscReal L[]) {
  PetscErrorCode ierr;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  if (!fe->hascoordinates) SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_ARG_WRONGSTATE,"Uniform coordinates not set");
  ierr = PetscMemcpy(L,fe->Luniform,sizeof fe->Luniform);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMFEGetInfo(DM dm,PetscInt *fedegree,PetscInt *level,PetscInt mlocal[],PetscInt Mglobal[],PetscInt procs[]) {
  PetscErrorCode ierr;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  if (fedegree) *fedegree = fe->degree;
  if (level)    *level = fe->grid->level;
  for (PetscInt i=0; i<3; i++) {
    if (mlocal)  mlocal[i] = fe->grid->m[i];
    if (Mglobal) Mglobal[i] = fe->grid->M[i];
    if (procs)   procs[i] = fe->grid->p[i];
  }
  PetscFunctionReturn(0);
}

// Injection: state restriction
//   primary order zero: high-frequency harmonics on fine grid alias at full amplitude
//   secondary order infinite: low-frequencies are exact on coarse grid
//
// All fine-grid processes must call this function.  Those not participating in coarse grid should pass Uc=NULL.
PetscErrorCode DMFEInject(DM dm,Vec Uf,Vec Uc)
{
  PetscErrorCode ierr;
  const PetscScalar *uf;
  PetscScalar *uc = NULL;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  if (fe->degree > 2) SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"fe->degree %D > 2",fe->degree);
  ierr = VecGetArrayRead(Uf,&uf);CHKERRQ(ierr);
  if (Uc) {ierr = VecGetArray(Uc,&uc);CHKERRQ(ierr);}

  ierr = PetscSFReduceBegin(fe->sfinject,fe->unit,uf,uc,MPI_REPLACE);CHKERRQ(ierr);
  ierr = PetscSFReduceEnd(fe->sfinject,fe->unit,uf,uc,MPI_REPLACE);CHKERRQ(ierr);

  ierr = VecRestoreArrayRead(Uf,&uf);CHKERRQ(ierr);
  if (Uc) {ierr = VecRestoreArray(Uc,&uc);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

// Interpolation: embedding of coarse space in fine space
//
// Processers not involved in coarse grid should pass Uc=NULL
PetscErrorCode DMFEInterpolate(DM dm,Vec Uc,Vec Uf)
{
  PetscErrorCode ierr;
  PetscScalar *uf;
  const PetscScalar *uc = NULL;
  const PetscInt *lM;
  FE fe;
  Vec Ufl;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  if (fe->degree > 2) SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"fe->degree %D > 2",fe->degree);

  ierr = DMGetLocalVector(dm,&Ufl);CHKERRQ(ierr);
  ierr = VecZeroEntries(Ufl);CHKERRQ(ierr);
  ierr = VecGetArray(Ufl,&uf);CHKERRQ(ierr);
  if (Uc) {
    ierr = VecGetArrayRead(Uc,&uc);CHKERRQ(ierr);
  }

  // Transpose of injection to populate C-points in fine grid
  ierr = PetscSFBcastBegin(fe->sfinjectLocal,fe->unit,uc,uf);CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(fe->sfinjectLocal,fe->unit,uc,uf);CHKERRQ(ierr);
  if (Uc) {
    ierr = VecRestoreArrayRead(Uc,&uc);CHKERRQ(ierr);
  }

  // Fill in missing entries in k
  lM = fe->lM;
  for (PetscInt i=0; i<lM[0]; i+=2) {
    for (PetscInt j=0; j<lM[1]; j+=2) {
      for (PetscInt k=1; k<lM[2]; k+=2*fe->degree) {
        for (PetscInt kk=0; kk<fe->degree; kk++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              uf[FEIdxL(fe,i,j,k+2*kk)*fe->dof+d] += fe->ref.interp[kk*(fe->degree+1)+l] * uf[FEIdxL(fe,i,j,k-1+2*l)*fe->dof+d];
            }
          }
        }
      }
    }
    for (PetscInt j=1; j<lM[1]; j+=2*fe->degree) {
      for (PetscInt jj=0; jj<fe->degree; jj++) {
        for (PetscInt k=0; k<lM[2]; k++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              uf[FEIdxL(fe,i,j+2*jj,k)*fe->dof+d] += fe->ref.interp[jj*(fe->degree+1)+l] * uf[FEIdxL(fe,i,j-1+2*l,k)*fe->dof+d];
            }
          }
        }
      }
    }
  }
  for (PetscInt i=1; i<lM[0]; i+=2*fe->degree) {
    for (PetscInt ii=0; ii<fe->degree; ii++) {
      for (PetscInt j=0; j<lM[1]; j++) {
        for (PetscInt k=0; k<lM[2]; k++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              uf[FEIdxL(fe,i+2*ii,j,k)*fe->dof+d] += fe->ref.interp[ii*(fe->degree+1)+l] * uf[FEIdxL(fe,i-1+2*l,j,k)*fe->dof+d];
            }
          }
        }
      }
    }
  }
  PetscLogFlops(lM[0]*lM[1]*lM[2]*fe->dof*(fe->degree+1)/8*2 + lM[0]*lM[1]*lM[2]*fe->dof*(fe->degree+1)/4*2 + lM[0]*lM[1]*lM[2]*fe->dof*(fe->degree+1)/2*2);
  ierr = VecRestoreArray(Ufl,&uf);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(dm,Ufl,INSERT_VALUES,Uf);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(dm,Ufl,INSERT_VALUES,Uf);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Ufl);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Residual Restriction: integration of fine points over coarse elements
//
// Processers not involved in coarse grid should pass Uc=NULL
PetscErrorCode DMFERestrict(DM dm,Vec Uf,Vec Uc)
{
  PetscErrorCode ierr;
  PetscScalar *ufl,*uc;
  const PetscScalar *uf;
  const PetscInt *lM;
  FE fe;
  Vec Ufl;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  if (fe->degree > 2) SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"fe->degree %D > 2",fe->degree);

  // Embed global space in local space (with zeros in ghost regions)
  ierr = DMGetLocalVector(dm,&Ufl);CHKERRQ(ierr);
  ierr = VecZeroEntries(Ufl);CHKERRQ(ierr);
  ierr = VecGetArray(Ufl,&ufl);CHKERRQ(ierr);
  ierr = VecGetArrayRead(Uf,&uf);CHKERRQ(ierr);
  for (PetscInt i=0; i<fe->om[0]; i++) {
    for (PetscInt j=0; j<fe->om[1]; j++) {
      for (PetscInt k=0; k<fe->om[2]; k++) {
        for (PetscInt d=0; d<fe->dof; d++) {
          ufl[FEIdxLs(fe,i,j,k)*fe->dof+d] = uf[FEIdxO(fe,i,j,k)*fe->dof+d];
        }
      }
    }
  }
  ierr = VecRestoreArrayRead(Uf,&uf);CHKERRQ(ierr);
  // Integrate using transpose of interpolation in each direction

  // Sum entries to C-points in i, then j, then k
  lM = fe->lM;
  for (PetscInt i=1; i<lM[0]; i+=2*fe->degree) {
    for (PetscInt ii=0; ii<fe->degree; ii++) {
      for (PetscInt j=0; j<lM[1]; j++) {
        for (PetscInt k=0; k<lM[2]; k++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              ufl[FEIdxL(fe,i-1+2*l,j,k)*fe->dof+d] += ufl[FEIdxL(fe,i+2*ii,j,k)*fe->dof+d] * fe->ref.interp[ii*(fe->degree+1)+l];
            }
          }
        }
      }
    }
  }
  for (PetscInt i=0; i<lM[0]; i+=2) {
    for (PetscInt j=1; j<lM[1]; j+=2*fe->degree) {
      for (PetscInt jj=0; jj<fe->degree; jj++) {
        for (PetscInt k=0; k<lM[2]; k++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              ufl[FEIdxL(fe,i,j-1+2*l,k)*fe->dof+d] += ufl[FEIdxL(fe,i,j+2*jj,k)*fe->dof+d] * fe->ref.interp[jj*(fe->degree+1)+l];
            }
          }
        }
      }
    }
    for (PetscInt j=0; j<lM[1]; j+=2) {
      for (PetscInt k=1; k<lM[2]; k+=2*fe->degree) {
        for (PetscInt kk=0; kk<fe->degree; kk++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              ufl[FEIdxL(fe,i,j,k-1+2*l)*fe->dof+d] += ufl[FEIdxL(fe,i,j,k+2*kk)*fe->dof+d] * fe->ref.interp[kk*(fe->degree+1)+l];
            }
          }
        }
      }
    }
  }
  PetscLogFlops(lM[0]*lM[1]*lM[2]*fe->dof*(fe->degree+1)/8*2 + lM[0]*lM[1]*lM[2]*fe->dof*(fe->degree+1)/4*2 + lM[0]*lM[1]*lM[2]*fe->dof*(fe->degree+1)/2*2);

  if (Uc) {
    ierr = VecGetArray(Uc,&uc);CHKERRQ(ierr);
  } else uc = NULL;

  // Inject from C-points in fine grid to coarse grid
  ierr = PetscSFReduceBegin(fe->sfinjectLocal,fe->unit,ufl,uc,MPIU_SUM);CHKERRQ(ierr);
  ierr = PetscSFReduceEnd(fe->sfinjectLocal,fe->unit,ufl,uc,MPIU_SUM);CHKERRQ(ierr);

  ierr = VecRestoreArray(Ufl,&ufl);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Ufl);CHKERRQ(ierr);

  if (Uc) {
    ierr = VecRestoreArray(Uc,&uc);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMFEZeroBoundaries(DM dm,Vec U) {
  PetscErrorCode ierr;
  PetscScalar *u;
  PetscInt gs[3],gM[3],*om;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  om = fe->om;
  for (PetscInt i=0; i<3; i++) {
    gs[i] = fe->grid->s[i]*fe->degree;
    gM[i] = fe->grid->M[i]*fe->degree+1; // global boundaries
  }
  ierr = VecGetArray(U,&u);CHKERRQ(ierr);
  for (PetscInt i=gs[0]; i<gs[0]+om[0]; i++) {
    for (PetscInt j=gs[1]; j<gs[1]+om[1]; j++) {
      for (PetscInt k=gs[2]; k<gs[2]+om[2]; k++) {
        if ((0<i && i<gM[0]-1) && (0<j && j<gM[1]-1) && (0<k && k<gM[2]-1)) continue;
        for (PetscInt d=0; d<fe->dof; d++) {
          u[FEIdxO(fe,i-gs[0],j-gs[1],k-gs[2])*fe->dof+d] = 0;
        }
      }
    }
  }
  ierr = VecRestoreArray(U,&u);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMFECoarsen(DM dm,DM *dmcoarse)
{
  PetscErrorCode ierr;
  FE fe,fecoarse = NULL;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  if (fe->dmcoarse) {
    *dmcoarse = fe->dmcoarse;
    PetscFunctionReturn(0);
  }
  if (fe->grid->coarse) {       /* My process participates in the coarse grid */
    ierr = DMCreateFE(fe->grid->coarse,fe->degree,fe->dof,&fe->dmcoarse);CHKERRQ(ierr);
    ierr = DMGetApplicationContext(fe->dmcoarse,&fecoarse);CHKERRQ(ierr);
  }

  if (1) {
    PetscInt nleaves,nroots,leaf,*ilocal,i,j,k;
    PetscSFNode *iremote;

    // Build injection for global spaces (only uses owned values)
    // ls[0]%2==0 when low corner of subdomain is even (C-point)
    nleaves = CeilDiv(fe->om[0]-(fe->ls[0]%2),2) * CeilDiv(fe->om[1]-(fe->ls[1]%2),2) * CeilDiv(fe->om[2]-(fe->ls[2]%2),2);
    ierr = PetscMalloc1(nleaves,&ilocal);CHKERRQ(ierr);
    ierr = PetscMalloc1(nleaves,&iremote);CHKERRQ(ierr);
    for (i=(fe->ls[0]%2),leaf=0; i<fe->om[0]; i+=2) {
      for (j=(fe->ls[1]%2); j<fe->om[1]; j+=2) {
        for (k=(fe->ls[2]%2); k<fe->om[2]; k+=2) {
          ilocal[leaf] = FEIdxO(fe,i,j,k);
          ierr = FEGridFindCoarseRankIndex(fe,i,j,k,&iremote[leaf].rank,&iremote[leaf].index);CHKERRQ(ierr);
          leaf++;
        }
      }
    }
    if (nleaves != leaf) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_PLIB,"nleaves %D != leaf %D",nleaves,leaf);
    ierr = PetscSFCreate(fe->grid->comm,&fe->sfinject);CHKERRQ(ierr);
    nroots = fecoarse ? fecoarse->om[0]*fecoarse->om[1]*fecoarse->om[2] : 0;
    ierr = PetscSFSetGraph(fe->sfinject,nroots,nleaves,ilocal,PETSC_OWN_POINTER,iremote,PETSC_OWN_POINTER);CHKERRQ(ierr);

    // Injection from fine local space to coarse global space
    nleaves = CeilDiv(fe->lM[0],2) * CeilDiv(fe->lM[1],2) * CeilDiv(fe->lM[2],2);
    ierr = PetscMalloc1(nleaves,&ilocal);CHKERRQ(ierr);
    ierr = PetscMalloc1(nleaves,&iremote);CHKERRQ(ierr);
    for (i=0,leaf=0; i<fe->lM[0]; i+=2) {
      for (j=0; j<fe->lM[1]; j+=2) {
        for (k=0; k<fe->lM[2]; k+=2) {
          ilocal[leaf] = Idx3(fe->lM,i,j,k);
          ierr = FEGridFindCoarseRankIndex(fe,i-fe->ls[0],j-fe->ls[1],k-fe->ls[2],&iremote[leaf].rank,&iremote[leaf].index);CHKERRQ(ierr);
          leaf++;
        }
      }
    }
    ierr = PetscSFCreate(fe->grid->comm,&fe->sfinjectLocal);CHKERRQ(ierr);
    nroots = fecoarse ? fecoarse->om[0]*fecoarse->om[1]*fecoarse->om[2] : 0;
    ierr = PetscSFSetGraph(fe->sfinjectLocal,nroots,nleaves,ilocal,PETSC_OWN_POINTER,iremote,PETSC_OWN_POINTER);CHKERRQ(ierr);
  }

  if (fe->hascoordinates) {
    DM dmc,dmc_coarse;
    Vec Xf,Xc;
    ierr = DMGetCoordinateDM(dm,&dmc);CHKERRQ(ierr);
    ierr = DMGetCoordinates(dm,&Xf);CHKERRQ(ierr);
    ierr = DMFECoarsen(dmc,&dmc_coarse);CHKERRQ(ierr);
    if (dmc_coarse) {
      ierr = DMCreateGlobalVector(dmc_coarse,&Xc);CHKERRQ(ierr);
    } else Xc = NULL;
    ierr = DMFEInject(dmc,Xf,Xc);CHKERRQ(ierr);
    if (dmc_coarse) {
      ierr = DMSetCoordinateDM(fe->dmcoarse,dmc_coarse);CHKERRQ(ierr);
      ierr = DMSetCoordinates(fe->dmcoarse,Xc);CHKERRQ(ierr);
      ierr = PetscMemcpy(fecoarse->Luniform,fe->Luniform,sizeof fe->Luniform);CHKERRQ(ierr);
      fecoarse->hascoordinates = PETSC_TRUE;
    }
    ierr = VecDestroy(&Xc);CHKERRQ(ierr);
  }

  *dmcoarse = fe->dmcoarse;
  PetscFunctionReturn(0);
}

static PetscErrorCode FEBasisEval(FE fe,PetscReal q,PetscReal B[],PetscReal D[])
{

  PetscFunctionBegin;
  switch (fe->degree) {
  case 1:
    B[0] = (1 - q)/2;
    B[1] = (1 + q)/2;
    if (D) {
      D[0] = -1./2;
      D[1] = 1./2;
    }
    break;
  case 2:
    B[0] = .5*(PetscSqr(q) - q);
    B[1] = 1 - PetscSqr(q);
    B[2] = .5*(PetscSqr(q) + q);
    if (D) {
      D[0] = q - .5;
      D[1] = -2*q;
      D[2] = q + .5;
    }
    break;
  default: SETERRQ1(fe->grid->comm,PETSC_ERR_SUP,"fe->degree %D",fe->degree);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode FESetUp(FE fe)
{
  PetscErrorCode ierr;
  PetscInt i,P,Q;

  PetscFunctionBegin;
  // Create reference element evaluation
  P = fe->degree+1;
  Q = fe->degree+1;
  ierr = PetscMalloc6(P*Q,&fe->ref.B,P*Q,&fe->ref.D,Q,&fe->ref.x,Q,&fe->ref.w,fe->degree*(fe->degree+1),&fe->ref.interp,Q*Q*Q,&fe->ref.w3);CHKERRQ(ierr);
  ierr = PetscDTGaussQuadrature(Q,-1,1,fe->ref.x,fe->ref.w);CHKERRQ(ierr);
  for (PetscInt i=0; i<Q; i++) {
    for (PetscInt j=0; j<Q; j++) {
      for (PetscInt k=0; k<Q; k++) {
        fe->ref.w3[(i*Q+j)*Q+k] = fe->ref.w[i] * fe->ref.w[j] * fe->ref.w[k];
      }
    }
  }
  for (i=0; i<Q; i++) {
    const PetscReal q = fe->ref.x[i];
    ierr = FEBasisEval(fe,q,&fe->ref.B[i*P],&fe->ref.D[i*P]);CHKERRQ(ierr);
  }
  // Interpolation to fill in nodes not nested in coarse grid
  // p=1: C0 -- f1 -- C2
  // p=2: C0 -- f1 -- C2 -- f3 -- C4
  for (i=0; i<fe->degree; i++) {
    const PetscReal floc[][2] = {{0.,0.},{0.,0.},{-.5,.5}};
    PetscReal q = floc[fe->degree][i];
    ierr = FEBasisEval(fe,q,&fe->ref.interp[i*(fe->degree+1)],NULL);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMFEGetTensorEval(DM dm,PetscInt *P,PetscInt *Q,const PetscReal **B,const PetscReal **D,const PetscReal **x,const PetscReal **w,const PetscReal **w3)
{
  PetscErrorCode ierr;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);

  *P = fe->degree + 1;
  *Q = fe->degree + 1;
  if (B) *B = fe->ref.B;
  if (D) *D = fe->ref.D;
  if (x) *x = fe->ref.x;
  if (w) *w = fe->ref.w;
  if (w3) *w3 = fe->ref.w3;
  PetscFunctionReturn(0);
}

PetscErrorCode DMFEGetNumElements(DM dm,PetscInt *nelems)
{
  PetscErrorCode ierr;
  FE fe;
  const PetscInt *m;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  m = fe->grid->m;
  *nelems = m[0]*m[1]*m[2];
  PetscFunctionReturn(0);
}

// Extract elements elem to elem+ne from local array u[grid i,j,k][0:dof], returning the result in y.
// y is padded by replicating last element in case of irregular ending
// vectorization-friendly ordering: y [0:dof] [0:(2*fedegree+1)^3] [0:ne]
PetscErrorCode DMFEExtractElements(DM dm,const PetscScalar *u,PetscInt elem,PetscInt ne,PetscScalar *y)
{
  PetscErrorCode ierr;
  FE fe;
  PetscInt P,fedegree,e;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  fedegree = fe->degree;
  P = fedegree + 1;

  for (e=elem; e<elem+ne; e++) {
    const PetscInt *m = fe->grid->m;
    PetscInt E = PetscMin(e,m[0]*m[1]*m[2]-1); // Last element replicated if we spill out of owned subdomain
    PetscInt i = E / (m[1]*m[2]);
    PetscInt j = (E - i*m[1]*m[2]) / m[2];
    PetscInt k = E - (i*m[1] + j)*m[2];
    PetscInt ii,jj,kk,d;
    for (d=0; d<fe->dof; d++) {
      for (ii=0; ii<P; ii++) {
        for (jj=0; jj<P; jj++) {
          for (kk=0; kk<P; kk++) {
            y[(((d*P+ii)*P+jj)*P+kk)*ne+(e-elem)] = u[FEIdxLs(fe,i*fedegree+ii,j*fedegree+jj,k*fedegree+kk)*fe->dof+d];
          }
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

// Sum/insert into elements elem:elem+ne in local vector u, using element contributions from y
// Any "elements" beyond the locally-owned part are ignored
PetscErrorCode DMFESetElements(DM dm,PetscScalar *u,PetscInt elem,PetscInt ne,InsertMode imode,DomainMode dmode,const PetscScalar *y)
{
  PetscErrorCode ierr;
  FE fe;
  PetscInt P,fedegree,e;
  const PetscInt *m;
  PetscInt gs[3],gM[3];

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);

  fedegree = fe->degree;
  P = fedegree + 1;
  m = fe->grid->m;
  for (PetscInt i=0; i<3; i++) {
    gs[i] = fe->grid->s[i]*fedegree;
    gM[i] = fe->grid->M[i]*fedegree+1; // global boundaries
  }

  for (e=elem; e<PetscMin(elem+ne,m[0]*m[1]*m[2]); e++) {
    PetscInt E = e;
    PetscInt i = E / (m[1]*m[2]);
    PetscInt j = (E - i*m[1]*m[2]) / m[2];
    PetscInt k = E - (i*m[1] + j)*m[2];
    PetscInt ii,jj,kk,d;
    for (d=0; d<fe->dof; d++) {
      for (ii=0; ii<P; ii++) {
        for (jj=0; jj<P; jj++) {
          for (kk=0; kk<P; kk++) {
            PetscInt iu = i*fedegree+ii,ju = j*fedegree+jj,ku = k*fedegree+kk;
            PetscInt src = (((d*P+ii)*P+jj)*P+kk)*ne + e-elem;
            PetscInt dst = FEIdxLs(fe,iu,ju,ku)*fe->dof+d;
            if ((0<gs[0]+iu && gs[0]+iu<gM[0]-1) && (0<gs[1]+ju && gs[1]+ju<gM[1]-1) && (0<gs[2]+ku && gs[2]+ku<gM[2]-1)) {
              if (PetscUnlikely((dmode & DOMAIN_INTERIOR) == 0)) continue;
            } else {
              if ((dmode & DOMAIN_EXTERIOR) == 0) continue;
            }
            if (imode == ADD_VALUES) u[dst] += y[src];
            else                     u[dst]  = y[src];
          }
        }
      }
    }
  }
  if (imode == ADD_VALUES) PetscLogFlops(fe->dof*P*P*P*(PetscMin(elem+ne,m[0]*m[1]*m[2])-elem));
  PetscFunctionReturn(0);
}

static PetscErrorCode FEDestroy(void **ctx)
{
  PetscErrorCode ierr;
  FE fe = (FE)*ctx;

  PetscFunctionBegin;
  if (!fe) PetscFunctionReturn(0);
  ierr = GridDestroy(&fe->grid);CHKERRQ(ierr);
  ierr = MPI_Type_free(&fe->unit);CHKERRQ(ierr);
  ierr = PetscSFDestroy(&fe->sf);CHKERRQ(ierr);
  ierr = PetscSFDestroy(&fe->sfinject);CHKERRQ(ierr);
  ierr = PetscSFDestroy(&fe->sfinjectLocal);CHKERRQ(ierr);
  ierr = DMDestroy(&fe->dmcoarse);CHKERRQ(ierr);
  ierr = PetscFree6(fe->ref.B,fe->ref.D,fe->ref.x,fe->ref.w,fe->ref.interp,fe->ref.w3);CHKERRQ(ierr);
  ierr = PetscFree(*ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Each process always owns the nodes in the negative direction (lower left).  The last process in each direction also
// owns its outer boundary.
PetscErrorCode DMCreateFE(Grid grid,PetscInt fedegree,PetscInt dof,DM *dmfe)
{
  PetscErrorCode ierr;
  PetscInt    *ilocal,i,j,k,nleaves,leaf,their_om[3],rneighbor_om[3];
  PetscSFNode *iremote;
  FE     fe;
  DM          dm;

  PetscFunctionBegin;
  ierr = PetscNew(&fe);CHKERRQ(ierr);
  grid->refct++;
  fe->grid = grid;
  fe->degree = fedegree;
  fe->dof = dof;
  for (i=0; i<3; i++) {
    PetscInt gs = 2*(grid->s[i]/2)*fedegree;                       // coordinate of start
    PetscInt ge = 2*CeilDiv(grid->s[i]+grid->m[i],2)*fedegree + 1; // one past coordinate of end
    fe->lM[i] = ge - gs;                                           // extent of local array
    fe->ls[i] = grid->s[i]*fedegree - gs;                          // first active part of local array
    fe->lm[i] = (grid->s[i]+grid->m[i])*fedegree + 1 - (gs + fe->ls[i]); // extent of active part in local array
    if (grid->neighborranks[1+(i==0)][1+(i==1)][1+(i==2)] >= 0) { // My neighbor exists so I don't own that fringe
      fe->om[i] = fe->lm[i]-1;
    } else {                    // I own my high boundary
      fe->om[i] = fe->lm[i];
    }
  }
  for (i=0; i<3; i++) {
    if (grid->s[i]+grid->m[i] == grid->M[i]) rneighbor_om[i] = -1;
    else {
      PetscInt ri = PartitionFind(grid->M[i],grid->p[i],grid->s[i]+grid->m[i]),s,m;
      ierr = PartitionGetRange(grid->M[i],grid->p[i],ri,&s,&m);CHKERRQ(ierr);
      rneighbor_om[i] = m*fedegree + (s+m == grid->M[i]);
    }
  }

  // Create neighbor scatter (roots=global, leaves=local)
  nleaves = fe->lm[0]*fe->lm[1]*fe->lm[2] - fe->om[0]*fe->om[1]*fe->om[2];
  ierr = PetscMalloc1(nleaves,&ilocal);CHKERRQ(ierr);
  ierr = PetscMalloc1(nleaves,&iremote);CHKERRQ(ierr);
  leaf = 0;
  for (i=0; i<fe->lm[0]; i++) {
    their_om[0] = i >= fe->om[0] ? rneighbor_om[0] : fe->om[0];
    for (j=0; j<fe->lm[1]; j++) {
      their_om[1] = j >= fe->om[1] ? rneighbor_om[1] : fe->om[1];
      for (k=0; k<fe->lm[2]; k++) {
        their_om[2] = k >= fe->om[2] ? rneighbor_om[2] : fe->om[2];
        if (i >= fe->om[0] || j >= fe->om[1] || k >= fe->om[2]) { // Someone else owns this vertex
          ilocal[leaf] = FEIdxLs(fe,i,j,k);
          iremote[leaf].rank = grid->neighborranks[1+(i>=fe->om[0])][1+(j>=fe->om[1])][1+(k>=fe->om[2])];
          iremote[leaf].index = ((i%fe->om[0])*their_om[1] + (j%fe->om[1]))*their_om[2] + (k%fe->om[2]);
          leaf++;
        }
      }
    }
  }
  ierr = PetscSFCreate(grid->comm,&fe->sf);CHKERRQ(ierr);
  ierr = PetscSFSetGraph(fe->sf,fe->om[0]*fe->om[1]*fe->om[2],nleaves,ilocal,PETSC_OWN_POINTER,iremote,PETSC_OWN_POINTER);CHKERRQ(ierr);
  ierr = MPI_Type_contiguous(dof,MPIU_SCALAR,&fe->unit);CHKERRQ(ierr);
  ierr = MPI_Type_commit(&fe->unit);CHKERRQ(ierr);

  ierr = FESetUp(fe);CHKERRQ(ierr);

  ierr = DMShellCreate(grid->comm,&dm);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(dm,fe);CHKERRQ(ierr);
  ierr = DMSetApplicationContextDestroy(dm,FEDestroy);CHKERRQ(ierr);
  ierr = DMShellSetLocalToGlobal(dm,DMLocalToGlobalBegin_FE,DMLocalToGlobalEnd_FE);CHKERRQ(ierr);
  ierr = DMShellSetGlobalToLocal(dm,DMGlobalToLocalBegin_FE,DMGlobalToLocalEnd_FE);CHKERRQ(ierr);
  ierr = DMShellSetCreateGlobalVector(dm,DMCreateGlobalVector_FE);CHKERRQ(ierr);
  ierr = DMShellSetCreateLocalVector(dm,DMCreateLocalVector_FE);CHKERRQ(ierr);
  *dmfe = dm;
  PetscFunctionReturn(0);
}
