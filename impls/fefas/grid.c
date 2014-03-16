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
  PetscInt cs[3],cm[3]; // My cells of coarse grid
  PetscMPIInt neighborranks[3][3][3];
};

// This type is hidden behind DM (not exposed publicly)
typedef struct FE_private *FE;
struct FE_private {
  Grid grid;
  FE fecoords;     // Dirty hack to work around application contexts not being destroyed
  PetscInt degree;      // Finite element polynomial degree
  PetscInt dof;         // Number of degrees of freedom per vertex
  PetscInt om[3];       // Array dimensions of owned part of global vectors
  PetscInt lm[3];       // Array dimensions of local vectors
  MPI_Datatype unit;
  PetscSF sf;
  PetscSegBuffer seg;
  struct {
    PetscReal *B;
    PetscReal *D;
    PetscReal *x;
    PetscReal *w;
  } ref;
};

static PetscInt FEIdxO(FE fe,PetscInt i,PetscInt j,PetscInt k) { return (i*fe->om[1] + j)*fe->om[2] + k; }
static PetscInt FEIdxL(FE fe,PetscInt i,PetscInt j,PetscInt k) { return (i*fe->lm[1] + j)*fe->lm[2] + k; }
static PetscInt FENeighborDim(FE fe,PetscInt l) {
  return fe->grid->s[l] + 2*fe->grid->m[l] < fe->grid->M[l] // Can we fit another domain the same size as mine without reaching boundary?
    ? fe->om[l]
    : fe->degree*(fe->grid->M[l] - (fe->grid->s[l]+fe->grid->m[l])) + 1;
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
  ierr = DMGetApplicationContext(dmc,&fe->fecoords);CHKERRQ(ierr); // Dirty hack so we can free it
  ierr = DMDestroy(&dmc);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMFEGetTensorEval(DM dm,PetscInt *P,PetscInt *Q,const PetscReal **B,const PetscReal **D,const PetscReal **x,PetscReal **w)
{
  PetscErrorCode ierr;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(dm,&fe);CHKERRQ(ierr);

  *P = fe->degree + 1;
  *Q = fe->degree + 1;
  if (!fe->ref.B) {
    PetscInt i;
    ierr = PetscMalloc4(*P * *Q,&fe->ref.B,*P * *Q,&fe->ref.D,*Q,&fe->ref.x,*Q,&fe->ref.w);CHKERRQ(ierr);
    ierr = PetscDTGaussQuadrature(*Q,-1,1,fe->ref.x,fe->ref.w);CHKERRQ(ierr);
    for (i=0; i<*Q; i++) {
      const PetscReal q = fe->ref.x[i];
      switch (fe->degree) {
      case 1:
        fe->ref.B[i*(*P)+0] = (1 - q)/2;
        fe->ref.D[i*(*P)+0] = -1./2;
        fe->ref.B[i*(*P)+1] = (1 + q)/2;
        fe->ref.D[i*(*P)+1] = 1./2;
        break;
      case 2:
        fe->ref.B[i*(*P)+0] = .5*(PetscSqr(q) - q);
	fe->ref.D[i*(*P)+0] = q - .5;
	fe->ref.B[i*(*P)+1] = 1 - PetscSqr(q);
	fe->ref.D[i*(*P)+1] = -2*q;
	fe->ref.B[i*(*P)+2] = .5*(PetscSqr(q) + q);
	fe->ref.D[i*(*P)+2] = q + .5;
        break;
      default: SETERRQ1(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"fe->degree %D",fe->degree);
      }
    }
  }
  if (B) *B = fe->ref.B;
  if (D) *D = fe->ref.D;
  if (x) *x = fe->ref.x;
  if (w) *w = fe->ref.w;
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
  P = 2*fedegree + 1;
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
            PetscInt src = (((d*P+ii)*P+jj)*P+kk) + e-elem;
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

// Each process always owns the nodes in the negative direction (lower left).  The last process in each direction also
// owns its outer boundary.
PetscErrorCode DMCreateFE(Grid grid,PetscInt fedegree,PetscInt dof,DM *dmfe)
{
  PetscErrorCode ierr;
  PetscInt    *ilocal,i,j,k,nleaves,leaf,their_om[3];
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
    fe->lm[i] = fedegree*grid->m[i]+1;
    if (grid->neighborranks[1+(i==0)][1+(i==1)][1+(i==2)] >= 0) { // My neighbor exists so I don't own that fringe
      fe->om[i] = fe->lm[i]-1;
    } else {                    // I own my high boundary
      fe->om[i] = fe->lm[i];
    }
  }
  nleaves = fe->lm[0]*fe->lm[1]*fe->lm[2] - fe->om[0]*fe->om[1]*fe->om[2];
  ierr = PetscMalloc1(nleaves,&ilocal);CHKERRQ(ierr);
  ierr = PetscMalloc1(nleaves,&iremote);CHKERRQ(ierr);
  leaf = 0;
  for (i=0; i<fe->lm[0]; i++) {
    their_om[0] = i >= fe->om[0] ? FENeighborDim(fe,0) : fe->om[0];
    for (j=0; j<fe->lm[1]; j++) {
      their_om[1] = j >= fe->om[1] ? FENeighborDim(fe,1) : fe->om[1];
      for (k=0; k<fe->lm[2]; k++) {
        their_om[2] = k >= fe->om[2] ? FENeighborDim(fe,2) : fe->om[2];
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

  ierr = DMShellCreate(grid->comm,&dm);CHKERRQ(ierr);
  ierr = DMSetApplicationContext(dm,fe);CHKERRQ(ierr);
  ierr = DMShellSetLocalToGlobal(dm,DMLocalToGlobalBegin_FE,DMLocalToGlobalEnd_FE);CHKERRQ(ierr);
  ierr = DMShellSetGlobalToLocal(dm,DMGlobalToLocalBegin_FE,DMGlobalToLocalEnd_FE);CHKERRQ(ierr);
  ierr = DMShellSetCreateGlobalVector(dm,DMCreateGlobalVector_FE);CHKERRQ(ierr);
  ierr = DMShellSetCreateLocalVector(dm,DMCreateLocalVector_FE);CHKERRQ(ierr);
  *dmfe = dm;
  PetscFunctionReturn(0);
}

PetscErrorCode DMDestroyFE(DM *dm)
{
  PetscErrorCode ierr;
  FE fe;

  PetscFunctionBegin;
  ierr = DMGetApplicationContext(*dm,&fe);CHKERRQ(ierr);
  while (fe) {
    FE fetmp;
    ierr = GridDestroy(&fe->grid);CHKERRQ(ierr);
    ierr = MPI_Type_free(&fe->unit);CHKERRQ(ierr);
    ierr = PetscSFDestroy(&fe->sf);CHKERRQ(ierr);
    ierr = PetscFree4(fe->ref.B,fe->ref.D,fe->ref.x,fe->ref.w);CHKERRQ(ierr);
    ierr = PetscSegBufferDestroy(&fe->seg);CHKERRQ(ierr);
    fetmp = fe;
    fe = fe->fecoords;
    ierr = PetscFree(fetmp);CHKERRQ(ierr);
  }
  ierr = DMDestroy(dm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
