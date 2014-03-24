#include "fefas.h"
#include <petscsf.h>
#include <petscdmshell.h>
#include <petscdt.h>
#include <stdint.h>

typedef uint64_t zcode;

struct Grid_private {
  PetscInt refct;
  MPI_Comm comm;
  MPI_Comm pcomm; // Communicator to talk to parent.
  Grid coarse;
  PetscInt M[3],p[3];
  PetscInt s[3],m[3];   // Owned cells of global grid
  PetscInt cs[3],cm[3]; // My cells of global coarse grid
  PetscInt Cs[3],Cm[3]; // Parent's ownership on global coarse grid
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
  PetscInt lm[3];       // Array dimensions of local vectors
  PetscInt Clm[3];      // Array dimensions of local coarse grid that we contribute to
  PetscInt Com[3];      // Array dimensions of owned coarse grid that we contribute to
  PetscInt cls[3];      // Start of array part that we contribute to coarse grid
  PetscInt rneighbor_om[3][3];
  PetscBool hascoordinates;
  MPI_Datatype unit;
  PetscSF sf;
  PetscSF sfinject;
  PetscSF sfinjectLocal;
  PetscSegBuffer seg;
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
static PetscInt FEIdxL(FE fe,PetscInt i,PetscInt j,PetscInt k) { return Idx3(fe->lm,i,j,k); }

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

PetscErrorCode GridCreate(MPI_Comm comm,const PetscInt M[3],const PetscInt p[3],const PetscInt pw[3],PetscInt cmax,Grid *grid)
{
  PetscErrorCode ierr;
  Grid g;
  PetscMPIInt size,rank;
  PetscInt Mc[3],j;
  zcode z;

  PetscFunctionBegin;
  ierr = PetscMalloc1(1,&g);CHKERRQ(ierr);
  ierr = PetscCommDuplicate(comm,&g->comm,NULL);CHKERRQ(ierr);
  g->pcomm = MPI_COMM_NULL;
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  if (size != p[0]*p[1]*p[2]) SETERRQ4(comm,PETSC_ERR_ARG_INCOMP,"Communicator size %d incompatible with process grid %D,%D,%D",size,p[0],p[1],p[2]);

  for (j=0; j<3; j++) {
    Mc[j] = M[j]/2;
    g->M[j] = M[j];
    g->p[j] = p[j];
  }

  // Find ownership
  z = ZCodeFromRank(rank,p);

  if (CeilDiv(M[0],p[0])*CeilDiv(M[1],p[1])*CeilDiv(M[2],p[2]) > cmax || size == 1) {
    // Coarse grid is on same process set, in-place
    if (M[0]%2 || M[1]%2 || M[2]%2) {
      if (size != 1) SETERRQ4(comm,PETSC_ERR_ARG_INCOMP,"Grid %D,%D,%D exceeds cmax %D, but cannot be coarsened",M[0],M[1],M[2],cmax);
      // Coarsest grid, on a single process
      g->coarse = NULL;
      for (j=0; j<3; j++) {
        g->m[j] = M[j];
        g->s[j] = 0;
        g->cm[j] = -1;
        g->cs[j] = -1;
      }
      for (j=0; j<27; j++) g->neighborranks[j/9%3][j/3%3][j%3] = -1;
    } else {
      ierr = GridCreate(comm,Mc,p,NULL,cmax,&g->coarse);CHKERRQ(ierr);
      for (j=0; j<3; j++) {
        g->m[j] = 2*g->coarse->m[j];
        g->s[j] = 2*g->coarse->s[j];
        g->cm[j] = g->coarse->m[j];
        g->cs[j] = g->coarse->s[j];
      }
      ierr = PetscMemcpy(g->neighborranks[0][0],g->coarse->neighborranks[0][0],sizeof g->neighborranks);CHKERRQ(ierr);
    }
  } else if (size > 1) { // z&07==0 will participate in coarse grid
    PetscMPIInt ri[3],prank,r;
    zcode t;
    MPI_Comm ccomm;
    PetscInt mpw[3],cpw[3],csm[6],cnt,nneighbors;

    if (M[0]%2 || M[1]%2 || M[2]%2) SETERRQ6(comm,PETSC_ERR_ARG_INCOMP,"Grid %D,%D,%D cannot be coarsened on process grid %D,%D,%D",M[0],M[1],M[2],p[0],p[1],p[2]);
    ierr = MPI_Comm_split(comm,z>>3,0,&g->pcomm);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(g->pcomm,&prank);CHKERRQ(ierr);
    ierr = MPI_Comm_split(comm,prank?MPI_UNDEFINED:0,0,&ccomm);CHKERRQ(ierr);

    for (j=0; j<3; j++) { // My contribution to coarse grid weight, only along axes
      zcode mask[] = {03,05,06};
      mpw[j] = z&mask[j] ? 0 : pw?pw[j]:1;
    }
    // Total coarse-grid share for this cluster.
    ierr = MPI_Allreduce(mpw,cpw,3,MPIU_INT,MPI_SUM,g->pcomm);CHKERRQ(ierr);
    if (!prank) { // rank 0 continues to coarse grid
      PetscInt cp[3];
      for (j=0; j<3; j++) cp[j] = CeilDiv(p[j],2);
      ierr = GridCreate(ccomm,Mc,cp,cpw,cmax,&g->coarse);CHKERRQ(ierr);
      ierr = MPI_Comm_free(&ccomm);CHKERRQ(ierr);
      for (j=0; j<3; j++) {
        csm[2*j+0] = g->coarse->s[j];
        csm[2*j+1] = g->coarse->m[j];
      }
    } else {
      g->coarse = NULL;
    }
    ierr = MPI_Bcast(csm,6,MPIU_INT,0,g->pcomm);CHKERRQ(ierr);
    for (j=0; j<3; j++) {
      int jbit = !!(z & (04>>j)); // Am I low or high within this dimension?
      PetscInt low_share = jbit ? cpw[j]-(pw?pw[j]:1) : (pw?pw[j]:1);
      PetscInt ms[2];
      ms[0] = CeilDiv(csm[2*j+1] * low_share,cpw[j]);
      ms[1] = csm[2*j+1] - ms[0];
      g->cm[j] = ms[jbit];
      g->cs[j] = csm[2*j+0] + jbit*ms[0];
      g->s[j] = 2*g->cs[j];
      g->m[j] = 2*g->cm[j];
      g->Cs[j] = csm[2*j+0];
      g->Cm[j] = csm[2*j+1];
    }
    ZCodeSplit(z,ri);

    for (cnt=0,nneighbors=0; cnt<27; cnt++) { // Count how many neighbors exist on grid
      PetscMPIInt i[3] = {ri[0]+cnt/9%3-1,ri[1]+cnt/3%3-1,ri[2]+cnt%3-1};
      if (0 <= i[0] && i[0] < p[0] && 0 <= i[1] && i[1] < p[1] && 0 <= i[2] && i[2] < p[2]) nneighbors++;
      g->neighborranks[cnt/9%3][cnt/3%3][cnt%3] = -1;
    }
    for (t=0,cnt=0,r=0; cnt<nneighbors; t++) {
      PetscMPIInt i[3];
      ZCodeSplit(t,i);
      if (i[0] < p[0] && i[1] < p[1] && i[2] < p[2]) {
        for (j=0; j<3; j++) i[j] -= ri[j];
        if (PetscAbs(i[0]) < 2 && PetscAbs(i[1]) < 2 && PetscAbs(i[2]) < 2) {
          g->neighborranks[i[0]+1][i[1]+1][i[2]+1] = r;
          cnt++;
        }
        r++;
      }
    }
  } else SETERRQ3(comm,PETSC_ERR_SUP,"Multiprocess coarse grid %D,%D,%D",M[0],M[1],M[2]);

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
  if ((*grid)->pcomm != MPI_COMM_NULL) {ierr = MPI_Comm_free(&(*grid)->pcomm);CHKERRQ(ierr);}
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
  ierr = VecSetSizes(*G,fe->lm[0]*fe->lm[1]*fe->lm[2]*fe->dof,PETSC_DETERMINE);CHKERRQ(ierr);
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
          l[FEIdxL(fe,i,j,k)*fe->dof+d] = g[FEIdxO(fe,i,j,k)*fe->dof+d];
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
  case INSERT_VALUES: op = MPI_REPLACE; break;
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
  MPI_Op op;

  PetscFunctionBegin;
  switch (imode) {
  case ADD_VALUES: op = MPIU_SUM; break;
  case INSERT_VALUES: op = MPI_REPLACE; break;
  default: SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"InsertMode");
  }
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  ierr = VecGetArrayRead(L,&l);CHKERRQ(ierr);
  ierr = VecGetArray(G,&g);CHKERRQ(ierr);

  // Add local part
  for (i=0; i<fe->om[0]; i++) {
    for (j=0; j<fe->om[1]; j++) {
      for (k=0; k<fe->om[2]; k++) {
        for (d=0; d<fe->dof; d++) {
          PetscInt src = FEIdxL(fe,i,j,k)*fe->dof+d;
          PetscInt dst = FEIdxO(fe,i,j,k)*fe->dof+d;
          if (imode == ADD_VALUES) g[dst] += l[src];
          else                     g[dst]  = l[src];
        }
      }
    }
  }
  ierr = PetscSFReduceEnd(fe->sf,fe->unit,l,g,op);CHKERRQ(ierr);
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
  for (PetscInt i=0; i<fe->lm[0]; i++) {
    for (PetscInt j=0; j<fe->lm[1]; j++) {
      for (PetscInt k=0; k<fe->lm[2]; k++) {
        x[FEIdxL(fe,i,j,k)*3+0] = L[0]*(fe->grid->s[0] + (PetscReal)i/fe->degree)/fe->grid->M[0];
        x[FEIdxL(fe,i,j,k)*3+1] = L[1]*(fe->grid->s[1] + (PetscReal)j/fe->degree)/fe->grid->M[1];
        x[FEIdxL(fe,i,j,k)*3+2] = L[2]*(fe->grid->s[2] + (PetscReal)k/fe->degree)/fe->grid->M[2];
      }
    }
  }
  ierr = VecRestoreArray(X,&x);CHKERRQ(ierr);
  ierr = DMSetCoordinateDM(dm,dmc);CHKERRQ(ierr);
  ierr = DMSetCoordinatesLocal(dm,X);CHKERRQ(ierr);
  fe->hascoordinates = PETSC_TRUE;
  ierr = DMDestroy(&dmc);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
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

  if (fe->sfinject) { // Communicate to parent process (data size expected to be relatively small)
    ierr = PetscSFReduceBegin(fe->sfinject,fe->unit,uf,uc,MPI_REPLACE);CHKERRQ(ierr);
    ierr = PetscSFReduceEnd(fe->sfinject,fe->unit,uf,uc,MPI_REPLACE);CHKERRQ(ierr);
  } else { // in-place coarsening on same process set
    FE fecoarse;
    ierr = DMGetApplicationContext(fe->dmcoarse,&fecoarse);CHKERRQ(ierr);
    for (PetscInt i=0; i<fecoarse->om[0]; i++) {
      for (PetscInt j=0; j<fecoarse->om[1]; j++) {
        for (PetscInt k=0; k<fecoarse->om[2]; k++) {
          for (PetscInt d=0; d<fe->dof; d++) {
            uc[FEIdxO(fecoarse,i,j,k)*fe->dof+d] = uf[FEIdxO(fe,i*2,j*2,k*2)*fe->dof+d];
          }
        }
      }
    }
  }

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
  const PetscInt *lm;
  FE fe;
  Vec Ucl,Ufl;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  if (fe->degree > 2) SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"fe->degree %D > 2",fe->degree);

  ierr = DMGetLocalVector(dm,&Ufl);CHKERRQ(ierr);
  ierr = VecZeroEntries(Ufl);CHKERRQ(ierr);
  ierr = VecGetArray(Ufl,&uf);CHKERRQ(ierr);
  if (Uc) {
    ierr = DMGetLocalVector(fe->dmcoarse,&Ucl);CHKERRQ(ierr);
    ierr = DMGlobalToLocalBegin(fe->dmcoarse,Uc,INSERT_VALUES,Ucl);CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(fe->dmcoarse,Uc,INSERT_VALUES,Ucl);CHKERRQ(ierr);
    ierr = VecGetArrayRead(Ucl,&uc);CHKERRQ(ierr);
  }

  // Transpose of injection to populate C-points in fine grid
  if (fe->sfinject) { // Communicate from parent process (data size expected to be relatively small)
    ierr = PetscSFBcastBegin(fe->sfinjectLocal,fe->unit,uc,uf);CHKERRQ(ierr);
    ierr = PetscSFBcastEnd(fe->sfinjectLocal,fe->unit,uc,uf);CHKERRQ(ierr);
  } else { // in-place on same process set
    FE fecoarse;
    ierr = DMGetApplicationContext(fe->dmcoarse,&fecoarse);CHKERRQ(ierr);
    for (PetscInt i=0; i<fecoarse->lm[0]; i++) {
      for (PetscInt j=0; j<fecoarse->lm[1]; j++) {
        for (PetscInt k=0; k<fecoarse->lm[2]; k++) {
          for (PetscInt d=0; d<fe->dof; d++) {
            uf[FEIdxL(fe,i*2,j*2,k*2)*fe->dof+d] = uc[FEIdxL(fecoarse,i,j,k)*fe->dof+d];
          }
        }
      }
    }
  }
  if (Uc) {
    ierr = VecRestoreArrayRead(Ucl,&uc);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(fe->dmcoarse,&Ucl);CHKERRQ(ierr);
  }

  // Fill in missing entries in k
  lm = fe->lm;
  for (PetscInt i=0; i<lm[0]; i+=2) {
    for (PetscInt j=0; j<lm[1]; j+=2) {
      for (PetscInt k=1; k<lm[2]; k+=2*fe->degree) {
        for (PetscInt kk=0; kk<fe->degree; kk++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              uf[FEIdxL(fe,i,j,k+2*kk)*fe->dof+d] += fe->ref.interp[kk*(fe->degree+1)+l] * uf[FEIdxL(fe,i,j,k-1+2*l)*fe->dof+d];
            }
          }
        }
      }
    }
    for (PetscInt j=1; j<lm[1]; j+=2*fe->degree) {
      for (PetscInt jj=0; jj<fe->degree; jj++) {
        for (PetscInt k=0; k<lm[2]; k++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              uf[FEIdxL(fe,i,j+2*jj,k)*fe->dof+d] += fe->ref.interp[jj*(fe->degree+1)+l] * uf[FEIdxL(fe,i,j-1+2*l,k)*fe->dof+d];
            }
          }
        }
      }
    }
  }
  for (PetscInt i=1; i<lm[0]; i+=2*fe->degree) {
    for (PetscInt ii=0; ii<fe->degree; ii++) {
      for (PetscInt j=0; j<lm[1]; j++) {
        for (PetscInt k=0; k<lm[2]; k++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              uf[FEIdxL(fe,i+2*ii,j,k)*fe->dof+d] += fe->ref.interp[ii*(fe->degree+1)+l] * uf[FEIdxL(fe,i-1+2*l,j,k)*fe->dof+d];
            }
          }
        }
      }
    }
  }
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
  const PetscInt *lm;
  FE fe;
  Vec Ucl,Ufl;

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
          ufl[FEIdxL(fe,i,j,k)*fe->dof+d] = uf[FEIdxO(fe,i,j,k)*fe->dof+d];
        }
      }
    }
  }
  ierr = VecRestoreArrayRead(Uf,&uf);CHKERRQ(ierr);

  // Integrate using transpose of interpolation in each direction

  // Fill in missing entries in k
  lm = fe->lm;
  for (PetscInt i=1; i<lm[0]; i+=2*fe->degree) {
    for (PetscInt ii=0; ii<fe->degree; ii++) {
      for (PetscInt j=0; j<lm[1]; j++) {
        for (PetscInt k=0; k<lm[2]; k++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              ufl[FEIdxL(fe,i-1+2*l,j,k)*fe->dof+d] += ufl[FEIdxL(fe,i+2*ii,j,k)*fe->dof+d] * fe->ref.interp[ii*(fe->degree+1)+l];
            }
          }
        }
      }
    }
  }
  for (PetscInt i=0; i<lm[0]; i+=2) {
    for (PetscInt j=1; j<lm[1]; j+=2*fe->degree) {
      for (PetscInt jj=0; jj<fe->degree; jj++) {
        for (PetscInt k=0; k<lm[2]; k++) {
          for (PetscInt l=0; l<fe->degree+1; l++) {
            for (PetscInt d=0; d<fe->dof; d++) {
              ufl[FEIdxL(fe,i,j-1+2*l,k)*fe->dof+d] += ufl[FEIdxL(fe,i,j+2*jj,k)*fe->dof+d] * fe->ref.interp[jj*(fe->degree+1)+l];
            }
          }
        }
      }
    }
    for (PetscInt j=0; j<lm[1]; j+=2) {
      for (PetscInt k=1; k<lm[2]; k+=2*fe->degree) {
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

  if (Uc) {
    ierr = DMGetLocalVector(fe->dmcoarse,&Ucl);CHKERRQ(ierr);
    ierr = VecZeroEntries(Ucl);CHKERRQ(ierr);
    ierr = VecGetArray(Ucl,&uc);CHKERRQ(ierr);
  } else uc = NULL;

  // Transpose of injection to populate C-points in fine grid
  if (fe->sfinject) { // Communicate from parent process (data size expected to be relatively small)
    ierr = PetscSFReduceBegin(fe->sfinjectLocal,fe->unit,ufl,uc,MPIU_SUM);CHKERRQ(ierr);
    ierr = PetscSFReduceEnd(fe->sfinjectLocal,fe->unit,ufl,uc,MPIU_SUM);CHKERRQ(ierr);
  } else { // in-place on same process set
    FE fecoarse;
    ierr = DMGetApplicationContext(fe->dmcoarse,&fecoarse);CHKERRQ(ierr);
    for (PetscInt i=0; i<fecoarse->lm[0]; i++) {
      for (PetscInt j=0; j<fecoarse->lm[1]; j++) {
        for (PetscInt k=0; k<fecoarse->lm[2]; k++) {
          for (PetscInt d=0; d<fe->dof; d++) {
            uc[FEIdxL(fecoarse,i,j,k)*fe->dof+d] += ufl[FEIdxL(fe,i*2,j*2,k*2)*fe->dof+d];
          }
        }
      }
    }
  }
  ierr = VecRestoreArray(Ufl,&ufl);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm,&Ufl);CHKERRQ(ierr);

  if (Uc) {
    ierr = VecRestoreArray(Ucl,&uc);CHKERRQ(ierr);
    ierr = DMLocalToGlobalBegin(fe->dmcoarse,Ucl,ADD_VALUES,Uc);CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(fe->dmcoarse,Ucl,ADD_VALUES,Uc);CHKERRQ(ierr);
    ierr = DMRestoreLocalVector(fe->dmcoarse,&Ucl);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

PetscErrorCode DMFECoarsen(DM dm,DM *dmcoarse)
{
  PetscErrorCode ierr;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);
  if (fe->dmcoarse) {
    *dmcoarse = fe->dmcoarse;
    PetscFunctionReturn(0);
  }
  if (fe->grid->coarse) {       /* My process participates in the coarse grid */
    FE fecoarse;
    ierr = DMCreateFE(fe->grid->coarse,fe->degree,fe->dof,&fe->dmcoarse);CHKERRQ(ierr);
    ierr = DMGetApplicationContext(fe->dmcoarse,&fecoarse);CHKERRQ(ierr);
    for (PetscInt i=0; i<3; i++) {
      fe->Clm[i] = fecoarse->lm[i];
      fe->Com[i] = fecoarse->om[i];
    }
  }
  if (fe->grid->pcomm != MPI_COMM_NULL) {
    ierr = MPI_Bcast(fe->Clm,3,MPIU_INT,0,fe->grid->pcomm);CHKERRQ(ierr);
    ierr = MPI_Bcast(fe->Com,3,MPIU_INT,0,fe->grid->pcomm);CHKERRQ(ierr);
  }
  for (PetscInt i=0; i<3; i++) fe->cls[i] = (fe->grid->s[i]/2 - fe->grid->Cs[i])*fe->degree;

  if (fe->grid->pcomm != MPI_COMM_NULL) { // Coarsening is not in-place so we need to build sfinject
    PetscInt nleaves,leaf,*ilocal,i,j,k;
    PetscSFNode *iremote;

    // Build injection for global spaces (only uses owned values)
    nleaves = CeilDiv(fe->om[0],2) * CeilDiv(fe->om[1],2) * CeilDiv(fe->om[2],2);
    ierr = PetscMalloc1(nleaves,&ilocal);CHKERRQ(ierr);
    ierr = PetscMalloc1(nleaves,&iremote);CHKERRQ(ierr);
    for (i=0,leaf=0; 2*i<fe->om[0]; i++) {
      for (j=0; 2*j<fe->om[1]; j++) {
        for (k=0; 2*k<fe->om[2]; k++) {
          ilocal[leaf] = FEIdxO(fe,2*i,2*j,2*k);
          iremote[leaf].rank = 0;  // rank 0 of pcomm always owns all of my coarse nodes
          iremote[leaf].index = Idx3(fe->Com,fe->cls[0]+i,fe->cls[1]+j,fe->cls[2]+k);
          leaf++;
        }
      }
    }
    ierr = PetscSFCreate(fe->grid->pcomm,&fe->sfinject);CHKERRQ(ierr);
    ierr = PetscSFSetGraph(fe->sfinject,fe->grid->coarse?fe->Com[0]*fe->Com[1]*fe->Com[2]:0,nleaves,ilocal,PETSC_OWN_POINTER,iremote,PETSC_OWN_POINTER);CHKERRQ(ierr);

    // Injection between local spaces (could/should be fused with local neighbor update)
    nleaves = CeilDiv(fe->lm[0],2) * CeilDiv(fe->lm[1],2) * CeilDiv(fe->lm[2],2);
    ierr = PetscMalloc1(nleaves,&ilocal);CHKERRQ(ierr);
    ierr = PetscMalloc1(nleaves,&iremote);CHKERRQ(ierr);
    for (i=0,leaf=0; 2*i<fe->lm[0]; i++) {
      for (j=0; 2*j<fe->lm[1]; j++) {
        for (k=0; 2*k<fe->lm[2]; k++) {
          ilocal[leaf] = FEIdxL(fe,2*i,2*j,2*k);
          iremote[leaf].rank = 0;  // rank 0 of pcomm always owns all of my coarse nodes
          iremote[leaf].index = Idx3(fe->Clm,fe->cls[0]+i,fe->cls[1]+j,fe->cls[2]+k);
          leaf++;
        }
      }
    }
    ierr = PetscSFCreate(fe->grid->pcomm,&fe->sfinjectLocal);CHKERRQ(ierr);
    ierr = PetscSFSetGraph(fe->sfinjectLocal,fe->grid->coarse?fe->Clm[0]*fe->Clm[1]*fe->Clm[2]:0,nleaves,ilocal,PETSC_OWN_POINTER,iremote,PETSC_OWN_POINTER);CHKERRQ(ierr);
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
  if (!fe->seg) {ierr = PetscSegBufferCreate(sizeof(PetscScalar),ne*fe->dof*P*P*P,&fe->seg);CHKERRQ(ierr);}

  for (e=elem; e<elem+ne; e++) {
    const PetscInt *m = fe->grid->m,*lm = fe->lm;
    PetscInt E = PetscMin(e,m[0]*m[1]*m[2]-1); // Last element replicated if we spill out of owned subdomain
    PetscInt i = E / (m[1]*m[2]);
    PetscInt j = (E - i*m[1]*m[2]) / m[2];
    PetscInt k = E - (i*m[1] + j)*m[2];
    PetscInt ii,jj,kk,d;
    for (d=0; d<fe->dof; d++) {
      for (ii=0; ii<P; ii++) {
        for (jj=0; jj<P; jj++) {
          for (kk=0; kk<P; kk++) {
            y[(((d*P+ii)*P+jj)*P+kk)*ne+(e-elem)] = u[(((i*fedegree+ii)*lm[1]+j*fedegree+jj)*lm[2]+k*fedegree+kk)*fe->dof+d];
          }
        }
      }
    }
  }
  PetscFunctionReturn(0);
}

// Sum/insert into elements elem:elem+ne in local vector u, using element contributions from y
// Any "elements" beyond the locally-owned part are ignored
PetscErrorCode DMFESetElements(DM dm,PetscScalar *u,PetscInt elem,PetscInt ne,InsertMode imode,const PetscScalar *y)
{
  PetscErrorCode ierr;
  FE fe;
  PetscInt P,fedegree,e;
  const PetscInt *m;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);

  fedegree = fe->degree;
  P = fedegree + 1;
  m = fe->grid->m;

  for (e=elem; e<PetscMin(elem+ne,m[0]*m[1]*m[2]); e++) {
    const PetscInt *lm = fe->lm;
    PetscInt E = e;
    PetscInt i = E / (m[1]*m[2]);
    PetscInt j = (E - i*m[1]*m[2]) / m[2];
    PetscInt k = E - (i*m[1] + j)*m[2];
    PetscInt ii,jj,kk,d;
    for (d=0; d<fe->dof; d++) {
      for (ii=0; ii<P; ii++) {
        for (jj=0; jj<P; jj++) {
          for (kk=0; kk<P; kk++) {
            PetscInt src = (((d*P+ii)*P+jj)*P+kk)*ne + e-elem;
            PetscInt dst = (((i*fedegree+ii)*lm[1]+j*fedegree+jj)*lm[2]+k*fedegree+kk)*fe->dof+d;
            if (imode == ADD_VALUES) u[dst] += y[src];
            else                     u[dst]  = y[src];
          }
        }
      }
    }
  }
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
  ierr = PetscSegBufferDestroy(&fe->seg);CHKERRQ(ierr);
  ierr = PetscFree(*ctx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

// Each process always owns the nodes in the negative direction (lower left).  The last process in each direction also
// owns its outer boundary.
PetscErrorCode DMCreateFE(Grid grid,PetscInt fedegree,PetscInt dof,DM *dmfe)
{
  PetscErrorCode ierr;
  PetscInt    *ilocal,i,j,k,nleaves,leaf,their_om[3];
  PetscSFNode *iremote;
  FE     fe;
  DM          dm;
  PetscMPIInt tag;

  PetscFunctionBegin;
  ierr = PetscNew(&fe);CHKERRQ(ierr);
  grid->refct++;
  fe->grid = grid;
  fe->degree = fedegree;
  fe->dof = dof;
  for (i=0; i<3; i++) {
    fe->lm[i] = fedegree*grid->m[i]+1;
    if (grid->neighborranks[1+(i==0)][1+(i==1)][1+(i==2)] >= 0) { // My neighbor exists so I don't own that fringe
      fe->om[i] = fe->lm[i]-1;
    } else {                    // I own my high boundary
      fe->om[i] = fe->lm[i];
    }
    fe->Com[i] = -1;
    fe->Clm[i] = -1;
  }
  ierr = PetscCommGetNewTag(grid->comm,&tag);CHKERRQ(ierr);
  for (i=0; i<3; i++) {
    MPI_Request req = MPI_REQUEST_NULL;
    PetscMPIInt rankR = grid->neighborranks[1+(i==0)][1+(i==1)][1+(i==2)];
    PetscMPIInt rankL = grid->neighborranks[1-(i==0)][1-(i==1)][1-(i==2)];
    if (rankR >= 0) {
      ierr = MPI_Irecv(fe->rneighbor_om[i],3,MPIU_INT,rankR,tag,grid->comm,&req);CHKERRQ(ierr);
    }
    if (rankL >= 0) {
      ierr = MPI_Send(fe->om,3,MPIU_INT,rankL,tag,grid->comm);CHKERRQ(ierr);
    }
    ierr = MPI_Wait(&req,MPI_STATUS_IGNORE);CHKERRQ(ierr);
  }

  // Create neighbor scatter (roots=global, leaves=local)
  nleaves = fe->lm[0]*fe->lm[1]*fe->lm[2] - fe->om[0]*fe->om[1]*fe->om[2];
  ierr = PetscMalloc1(nleaves,&ilocal);CHKERRQ(ierr);
  ierr = PetscMalloc1(nleaves,&iremote);CHKERRQ(ierr);
  leaf = 0;
  for (i=0; i<fe->lm[0]; i++) {
    their_om[0] = i >= fe->om[0] ? fe->rneighbor_om[0][0] : fe->om[0];
    for (j=0; j<fe->lm[1]; j++) {
      their_om[1] = j >= fe->om[1] ? fe->rneighbor_om[1][1] : fe->om[1];
      for (k=0; k<fe->lm[2]; k++) {
        their_om[2] = k >= fe->om[2] ? fe->rneighbor_om[2][2] : fe->om[2];
        if (i >= fe->om[0] || j >= fe->om[1] || k >= fe->om[2]) { // Someone else owns this vertex
          ilocal[leaf] = FEIdxL(fe,i,j,k);
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
