#include "fefas.h"
#include "petsctime.h"
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

static PetscErrorCode ReportMemoryUsage(MPI_Comm comm,PetscLogDouble memused,PetscLogDouble memavail) {
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Allreduce(MPI_IN_PLACE,&memused,1,MPI_DOUBLE,MPI_MAX,comm);CHKERRQ(ierr);
  ierr = MPI_Allreduce(MPI_IN_PLACE,&memavail,1,MPI_DOUBLE,MPI_MAX,comm);CHKERRQ(ierr);
  if (memavail >= 0) {
    ierr = PetscPrintf(comm,"Max memory per MPI rank: %f GB, %f GB available\n",memused*1e-9,memavail*1e-9);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(comm,"Max memory per MPI rank: %f GB\n",memused*1e-9);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode SampleOnGrid(MPI_Comm comm,Op op,const PetscInt M[3],const PetscInt smooth[2],PetscInt nrepeat,PetscLogDouble mintime,PetscLogDouble *memused,PetscLogDouble *memavail,PetscBool monitor) {
  PetscErrorCode ierr;
  PetscInt pgrid[3],cmax,fedegree,dof,nlevels,M_max;
  PetscMPIInt nranks;
  Grid grid;
  DM dm;
  Vec U,F;
  MG mg;
  PetscReal L[3];
  PetscBool affine;
#ifdef USE_HPM
  char eventname[256];
#endif

  PetscFunctionBegin;
  ierr = OpGetFEDegree(op,&fedegree);CHKERRQ(ierr);
  ierr = OpGetDof(op,&dof);CHKERRQ(ierr);

  ierr = MPI_Comm_size(comm,&nranks);CHKERRQ(ierr);
  ierr = ProcessGridFindSquarest(nranks,pgrid);CHKERRQ(ierr);

  // It would make sense to either use a different coarsening criteria (perhaps even specified by the sampler).  On
  // large numbers of processes, the coarse grids should be square enough that 192 is a good threshold size.
  cmax = 192;

  ierr = GridCreate(comm,M,pgrid,cmax,&grid);CHKERRQ(ierr);
  ierr = GridGetNumLevels(grid,&nlevels);CHKERRQ(ierr);

  ierr = DMCreateFE(grid,fedegree,dof,&dm);CHKERRQ(ierr);
  M_max = PetscMax(M[0],PetscMax(M[1],M[2]));
  L[0] = M[0]*1./M_max;
  L[1] = M[1]*1./M_max;
  L[2] = M[2]*1./M_max;
  ierr = DMFESetUniformCoordinates(dm,L);CHKERRQ(ierr);
  ierr = OpGetAffineOnly(op,&affine);CHKERRQ(ierr);
  if (!affine) {ierr = DMCoordDistort(dm,L);CHKERRQ(ierr);}

  ierr = DMCreateGlobalVector(dm,&U);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&F);CHKERRQ(ierr);
  ierr = OpForcing(op,dm,F);CHKERRQ(ierr);

  ierr = MGCreate(op,dm,nlevels,&mg);CHKERRQ(ierr);
  ierr = MGMonitorSet(mg,monitor);CHKERRQ(ierr);
  ierr = MGSetUpPC(mg);CHKERRQ(ierr);

#ifdef USE_HPM
  ierr = PetscSNPrintf(eventname,sizeof eventname,"Solve G[%D %D %D]",M[0],M[1],M[2]);CHKERRQ(ierr);
  HPM_Start(eventname);
#endif
  PetscInt i = 0;
  PetscLogDouble sampletime = 0;
  while ( (i<nrepeat) || (sampletime < mintime) ) {
    PetscLogDouble t0,t1,elapsed,flops,eqs;
    ierr = VecZeroEntries(U);CHKERRQ(ierr);
    ierr = MPI_Barrier(comm);CHKERRQ(ierr);
    ierr = PetscTime(&t0);CHKERRQ(ierr);
    flops = petsc_TotalFlops;
    ierr = MGFCycle(op,mg,smooth[0],smooth[1],F,U);CHKERRQ(ierr);
    ierr = PetscTime(&t1);CHKERRQ(ierr);
    flops = petsc_TotalFlops - flops;
    elapsed = t1 - t0;
    ierr = MPI_Allreduce(MPI_IN_PLACE,&elapsed,1,MPI_DOUBLE,MPI_MAX,comm);CHKERRQ(ierr);
    ierr = MPI_Allreduce(MPI_IN_PLACE,&flops,1,MPI_DOUBLE,MPI_SUM,comm);CHKERRQ(ierr);
    eqs = (double)(M[0]*fedegree+1)*(M[1]*fedegree+1)*(M[2]*fedegree+1)*dof;
    ierr = PetscPrintf(comm,"Q%D G[%5D%5D%5D] P[%3D%3D%3D] %10.3e s  %10f GF  %10f MEq/s\n",fedegree,M[0],M[1],M[2],pgrid[0],pgrid[1],pgrid[2],t1-t0,flops/elapsed*1e-9,eqs/elapsed*1e-6);CHKERRQ(ierr);
    i++;
    sampletime += elapsed;
  }
#ifdef USE_HPM
  HPM_Stop(eventname);
#endif

  if (memused) {ierr = MemoryGetUsage(memused,memavail);CHKERRQ(ierr);}
  ierr = MGDestroy(&mg);CHKERRQ(ierr);
  ierr = VecDestroy(&U);CHKERRQ(ierr);
  ierr = VecDestroy(&F);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = GridDestroy(&grid);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode RunSample() {
  PetscErrorCode ierr;
  Op op;
  PetscInt pgrid[3],smooth[2] = {3,1},two = 2,maxsamples = 6,repeat = 5,nsamples,(*gridsize)[3];
  PetscReal local[2] = {100,10000};
  PetscReal mintime = 1;
  PetscLogDouble memused,memavail;
  PetscMPIInt nranks;
  MPI_Comm comm = PETSC_COMM_WORLD;

  PetscFunctionBegin;
  ierr = PetscOptionsBegin(comm,NULL,"FMG Performance Sampler options",NULL);CHKERRQ(ierr);
  ierr = PetscOptionsRealArray("-local","range of local problem sizes","",local,&two,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-maxsamples","maximum number of samples across range","",maxsamples,&maxsamples,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-repeat","Minimum number of repetitions for each problem size","",repeat,&repeat,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-mintime","Minimum interval (in seconds) for repeatedly solving each problem size","",mintime,&mintime,NULL);CHKERRQ(ierr);
  two = 2;
  ierr = PetscOptionsIntArray("-smooth","V- and F-cycle pre,post smoothing","",smooth,&two,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  ierr = OpCreateFromOptions(comm,&op);CHKERRQ(ierr);

  ierr = MPI_Comm_size(comm,&nranks);CHKERRQ(ierr);
  ierr = ProcessGridFindSquarest(nranks,pgrid);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"Finite Element FAS Performance Sampler on process grid [%D %D %D] = %d\n",pgrid[0],pgrid[1],pgrid[2],nranks);CHKERRQ(ierr);

  ierr = SampleGridRangeCreate(nranks,(PetscReal)local[0],(PetscReal)local[1],maxsamples,&nsamples,(PetscInt**)&gridsize);CHKERRQ(ierr);

  ierr = MemoryGetUsage(&memused,&memavail);CHKERRQ(ierr);
  ierr = ReportMemoryUsage(comm,memused,memavail);CHKERRQ(ierr);

  ierr = PetscPrintf(comm,"Small Test G[%5D%5D%5D]\n",gridsize[nsamples-1][0],gridsize[nsamples-1][1],gridsize[nsamples-1][2]);CHKERRQ(ierr);
  ierr = SampleOnGrid(comm,op,gridsize[nsamples-1],smooth,1,0,NULL,NULL,PETSC_FALSE);CHKERRQ(ierr);
  ierr = PetscPrintf(comm,"Large Test G[%5D%5D%5D]\n",gridsize[0][0],gridsize[0][1],gridsize[0][2]);CHKERRQ(ierr);
  ierr = SampleOnGrid(comm,op,gridsize[0],smooth,1,0,&memused,&memavail,PETSC_TRUE);CHKERRQ(ierr);

  ierr = ReportMemoryUsage(comm,memused,memavail);CHKERRQ(ierr);

  ierr = PetscPrintf(comm,"Starting performance sampling\n");CHKERRQ(ierr);
  for (PetscInt i=nsamples-1; i>=0; i--) {
    ierr = SampleOnGrid(comm,op,gridsize[i],smooth,repeat,mintime,NULL,NULL,PETSC_FALSE);CHKERRQ(ierr);
  }

  ierr = PetscFree(gridsize);CHKERRQ(ierr);
  ierr = OpDestroy(&op);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
