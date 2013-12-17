//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include "../timer.h"
//------------------------------------------------------------------------------------------------------------------------------
void residual(domain_type * domain, int level,  int res_id, int phi_id, int rhs_id, double a, double b){
  // form coefficients for 13pt
  double coef0  = -90.0/12.0;
  double coef1  =  16.0/12.0;
  double coef2  =  -1.0/12.0;

  // exchange the boundary for x in prep for Ax...
  // for 13-point stencil, needs a 2-deep ghost zone for faces only
  exchange_boundary(domain,level,phi_id,1,0,0);

  // now do residual/restriction proper...
  uint64_t _timeStart = CycleTime();
  int CollaborativeThreadingBoxSize = 100000; // i.e. never
  #ifdef __COLLABORATIVE_THREADING
    CollaborativeThreadingBoxSize = 1 << __COLLABORATIVE_THREADING;
  #endif
  int omp_across_boxes = (domain->subdomains[0].levels[level].dim.i <  CollaborativeThreadingBoxSize);
  int omp_within_a_box = (domain->subdomains[0].levels[level].dim.i >= CollaborativeThreadingBoxSize);
  int box;

  #pragma omp parallel for private(box) if(omp_across_boxes)
  for(box=0;box<domain->subdomains_per_rank;box++){
    int i,j,k;
    int pencil = domain->subdomains[box].levels[level].pencil;int pencil2 = pencil*2;
    int  plane = domain->subdomains[box].levels[level].plane; int  plane2 =  plane*2;
    int ghosts = domain->subdomains[box].levels[level].ghosts;
    int  dim_k = domain->subdomains[box].levels[level].dim.k;
    int  dim_j = domain->subdomains[box].levels[level].dim.j;
    int  dim_i = domain->subdomains[box].levels[level].dim.i;
    double h2inv = 1.0/(domain->h[level]*domain->h[level]);
    double * __restrict__ phi    = domain->subdomains[box].levels[level].grids[  phi_id] + ghosts*(1+pencil+plane); // i.e. [0] = first non ghost zone point
    double * __restrict__ rhs    = domain->subdomains[box].levels[level].grids[  rhs_id] + ghosts*(1+pencil+plane);
    double * __restrict__ alpha  = domain->subdomains[box].levels[level].grids[__alpha ] + ghosts*(1+pencil+plane);
    double * __restrict__ beta_i = domain->subdomains[box].levels[level].grids[__beta_i] + ghosts*(1+pencil+plane);
    double * __restrict__ beta_j = domain->subdomains[box].levels[level].grids[__beta_j] + ghosts*(1+pencil+plane);
    double * __restrict__ beta_k = domain->subdomains[box].levels[level].grids[__beta_k] + ghosts*(1+pencil+plane);
    double * __restrict__ res    = domain->subdomains[box].levels[level].grids[  res_id] + ghosts*(1+pencil+plane);

    #pragma omp parallel for private(k,j,i) if(omp_within_a_box) collapse(2)
    for(k=0;k<dim_k;k++){
    for(j=0;j<dim_j;j++){
    for(i=0;i<dim_i;i++){
      int ijk = i + j*pencil + k*plane;
      double helmholtz =  a*phi[ijk]
                         -b*h2inv*(
                              coef2*(phi[ijk-plane2          ] +
                                     phi[ijk       -pencil2  ] +
                                     phi[ijk               -2] +
                                     phi[ijk               +2] +
                                     phi[ijk       +pencil2  ] +
                                     phi[ijk+plane2          ] ) +
                              coef1*(phi[ijk-plane           ] +
                                     phi[ijk       -pencil   ] +
                                     phi[ijk               -1] +
                                     phi[ijk               +1] +
                                     phi[ijk       +pencil   ] +
                                     phi[ijk+plane           ] ) +
                              coef0*(phi[ijk                 ] )
                         );
      res[ijk] = rhs[ijk]-helmholtz;
    }}}
  }
  domain->cycles.residual[level] += (uint64_t)(CycleTime()-_timeStart);
}

void residual_and_restriction(domain_type *domain, int level_f, int phi_id, int rhs_id, int level_c, int res_id, double a, double b){
  printf("Error, residual_and_restriction has not been implemented...\n");
  fflush(stdout);
  exit(0);
}
//------------------------------------------------------------------------------------------------------------------------------
