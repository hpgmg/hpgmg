#include "fefas.h"
#include <stdint.h>
#include <inttypes.h>

static const PetscInt _coarse_grids[][3] = {{2,2,2}, {2,2,3}, {2,2,4}, {2,3,3}, {2,3,4}, {3,3,3}, {3,3,4}, {3,4,4}};
static const size_t _n_coarse_grids = sizeof _coarse_grids / sizeof _coarse_grids[0];

static PetscInt CeilDiv(PetscInt a,PetscInt b) {return a/b + !!(a%b);}

struct Sampler_private {
  PetscInt cursample;
  PetscInt nsamples;
  PetscInt *samples;
};

// Find the squarest grid, with best[0] <= best[1] <= best[2].
PetscErrorCode ProcessGridFindSquarest(PetscMPIInt nranks,PetscInt best[3]) {
  PetscMPIInt target,s,a,b,c;

  PetscFunctionBegin;
  target = PetscCeilReal(PetscPowReal(nranks,1./3));
  if (target*target*target > nranks) target--; // if ceil was overzealous
  for (a=target; a>=1; a--) {
    if (nranks%a) continue; // Not a candidate factor
    s = PetscCeilReal(PetscSqrtReal(nranks/a));
    if (s*s > nranks/a) s--; // If our ceil was overzealous
    for (b=s; b>=a; b--) {
      if (nranks/a % b == 0) {  // The first proper divisor is the one I want
        c = nranks/a/b;
        best[0] = a;
        best[1] = b;
        best[2] = c;
        PetscFunctionReturn(0);
      }
    }
  }
  PetscFunctionReturn(0);
}

PetscInt ProcessGridNumLevels(const PetscInt p[]) {
  PetscMPIInt pmax = p[2],plev;
  if (pmax < 1) return -1;
  for (plev=0; pmax > 1; pmax = CeilDiv(pmax,2)) plev++;
  return plev;
}

int64_t SampleGridNumElements(const PetscInt M[]) {
  return (int64_t)M[0]*M[1]*M[2];
}

static PetscErrorCode FindCompatibleProblemSize(PetscMPIInt nranks,PetscInt targetlocal,PetscInt size[3]) {
  PetscErrorCode ierr;
  PetscInt plev,pgrid[3];
  int64_t best = -1;

  PetscFunctionBegin;
  for (PetscInt i=0; i<3; i++) size[i] = -1;

  // We need our computational grids to have at least as many levels as the process grids
  ierr = ProcessGridFindSquarest(nranks,pgrid);CHKERRQ(ierr);
  plev = ProcessGridNumLevels(pgrid);

  for (PetscInt i=0; i<_n_coarse_grids; i++) {
    const PetscInt *c = _coarse_grids[i];
    for (int64_t gsize = (int64_t)c[0]*c[1]*c[2]*PetscPowInt(8,plev),lev=plev; gsize < (int64_t)targetlocal*nranks; gsize *= 8,lev++) {
      if (gsize > best) {
        best = gsize;
        for (PetscInt j=0; j<3; j++) size[j] = c[j] * PetscPowInt(2,lev);
      }
    }
  }
  PetscFunctionReturn(0);
}

// gridsizes is array of length nsamples, to be PetscFree'd by caller
PetscErrorCode SampleGridRangeCreate(PetscMPIInt nranks,PetscInt minlocal,PetscInt maxlocal,PetscInt maxsamples,PetscInt *nsamples,PetscInt **gridsizes) {
  PetscErrorCode ierr;
  int64_t target;
  PetscInt gsize[100][3],n;

  PetscFunctionBegin;
  if (maxsamples < 2) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_INCOMP,"The max number of samples must be at least 2");

  // Build a list of compatible grid sizes in descending order
  for (target=maxlocal,n=0; target>=minlocal; n++) {
    ierr = FindCompatibleProblemSize(nranks,target,gsize[n]);CHKERRQ(ierr);
    if (SampleGridNumElements(gsize[n]) < minlocal) {
      n--;
      break;
    }
    target = (SampleGridNumElements(gsize[n]) - 1)/nranks;
  }

  // Filter the list by greedily removing interior sample locations whose removal would leave behind the smallest
  // possible ratio between successive sizes
  while (n > maxsamples) {
    PetscInt loc = -1;
    double   ratio = 1e10;
    for (PetscInt i=1; i<n-1; i++) {
      double r = (double)SampleGridNumElements(gsize[i-1]) / SampleGridNumElements(gsize[i+1]);
      if (r < ratio) {
        loc = i;
        ratio = r;
      }
    }
    ierr = PetscMemmove(gsize[loc],gsize[loc+1],(char*)gsize[n]-(char*)gsize[loc+1]);CHKERRQ(ierr);
    n--;
  }

  *nsamples = n;
  ierr = PetscMalloc(n*sizeof gsize[0],gridsizes);CHKERRQ(ierr);
  ierr = PetscMemcpy(*gridsizes,gsize[0],n*sizeof gsize[0]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
