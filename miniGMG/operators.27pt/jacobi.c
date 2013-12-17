//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include "../timer.h"
//------------------------------------------------------------------------------------------------------------------------------
void smooth(domain_type * domain, int level, int phi_id, int rhs_id, double a, double b){
  if(numSmooths&1){
    printf("error - numSmooths must be even...\n");
    exit(0);
  }

  int CollaborativeThreadingBoxSize = 100000; // i.e. never
  #ifdef __COLLABORATIVE_THREADING
    CollaborativeThreadingBoxSize = 1 << __COLLABORATIVE_THREADING;
  #endif
  int omp_across_boxes = (domain->subdomains[0].levels[level].dim.i <  CollaborativeThreadingBoxSize);
  int omp_within_a_box = (domain->subdomains[0].levels[level].dim.i >= CollaborativeThreadingBoxSize);

  int box,s;
  int ghosts = domain->ghosts;
  double TwoThirds = 2.0/3.0;
  // form coefficients for 27pt
  double center  = -128.0/30.0;
  double faces   =   14.0/30.0;
  double edges   =    3.0/30.0;
  double corners =    1.0/30.0;

  // if communication-avoiding, need RHS for stencils in ghost zones
  if(ghosts>1)exchange_boundary(domain,level,rhs_id,1,1,1);

  for(s=0;s<numSmooths;s+=ghosts){
    // Jacobi ping pongs between phi and __temp
    if((s&1)==0)exchange_boundary(domain,level,phi_id,1,1,1);  // 27-point stencil, needs 1-deep ghost zone for faces, edges, and corners
           else exchange_boundary(domain,level,__temp,1,1,1);  // 27-point stencil, needs 1-deep ghost zone for faces, edges, and corners

    // now do ghosts communication-avoiding smooths on each box...
    uint64_t _timeStart = CycleTime();

    #pragma omp parallel for private(box) if(omp_across_boxes)
    for(box=0;box<domain->subdomains_per_rank;box++){
      int i,j,k,ss;
      int pencil = domain->subdomains[box].levels[level].pencil;
      int  plane = domain->subdomains[box].levels[level].plane;
      int ghosts = domain->subdomains[box].levels[level].ghosts;
      int  dim_k = domain->subdomains[box].levels[level].dim.k;
      int  dim_j = domain->subdomains[box].levels[level].dim.j;
      int  dim_i = domain->subdomains[box].levels[level].dim.i;
      double h2inv = 1.0/(domain->h[level]*domain->h[level]);
      double * __restrict__ rhs    = domain->subdomains[box].levels[level].grids[  rhs_id] + ghosts*(1+pencil+plane);
      double * __restrict__ alpha  = domain->subdomains[box].levels[level].grids[__alpha ] + ghosts*(1+pencil+plane);
      double * __restrict__ beta_i = domain->subdomains[box].levels[level].grids[__beta_i] + ghosts*(1+pencil+plane);
      double * __restrict__ beta_j = domain->subdomains[box].levels[level].grids[__beta_j] + ghosts*(1+pencil+plane);
      double * __restrict__ beta_k = domain->subdomains[box].levels[level].grids[__beta_k] + ghosts*(1+pencil+plane);
      double * __restrict__ lambda = domain->subdomains[box].levels[level].grids[__lambda] + ghosts*(1+pencil+plane);
  
      int ghostsToOperateOn=ghosts-1;
      for(ss=s;ss<s+ghosts;ss++,ghostsToOperateOn--){
        double * __restrict__ phi;
        double * __restrict__ phi_new;
              if((ss&1)==0){phi    = domain->subdomains[box].levels[level].grids[  phi_id] + ghosts*(1+pencil+plane);
                            phi_new= domain->subdomains[box].levels[level].grids[  __temp] + ghosts*(1+pencil+plane);}
                       else{phi    = domain->subdomains[box].levels[level].grids[  __temp] + ghosts*(1+pencil+plane);
                            phi_new= domain->subdomains[box].levels[level].grids[  phi_id] + ghosts*(1+pencil+plane);}
        #pragma omp parallel for private(k,j,i) if(omp_within_a_box) collapse(2)
        for(k=0-ghostsToOperateOn;k<dim_k+ghostsToOperateOn;k++){
        for(j=0-ghostsToOperateOn;j<dim_j+ghostsToOperateOn;j++){
        for(i=0-ghostsToOperateOn;i<dim_i+ghostsToOperateOn;i++){
          int ijk = i + j*pencil + k*plane;
          double helmholtz =  a*phi[ijk]
                             -b*h2inv*(
                                corners*(phi[ijk-plane-pencil-1] +
                                         phi[ijk-plane-pencil+1] +
                                         phi[ijk-plane+pencil-1] +
                                         phi[ijk-plane+pencil+1] +
                                         phi[ijk+plane-pencil-1] +
                                         phi[ijk+plane-pencil+1] +
                                         phi[ijk+plane+pencil-1] +
                                         phi[ijk+plane+pencil+1] ) +
                                  edges*(phi[ijk-plane-pencil  ] +
                                         phi[ijk-plane       -1] +
                                         phi[ijk-plane       +1] +
                                         phi[ijk-plane+pencil  ] +
                                         phi[ijk      -pencil-1] +
                                         phi[ijk      -pencil+1] +
                                         phi[ijk      +pencil-1] +
                                         phi[ijk      +pencil+1] +
                                         phi[ijk+plane-pencil  ] +
                                         phi[ijk+plane       -1] +
                                         phi[ijk+plane       +1] +
                                         phi[ijk+plane+pencil  ] ) +
                                  faces*(phi[ijk-plane         ] +
                                         phi[ijk      -pencil  ] +
                                         phi[ijk             -1] +
                                         phi[ijk             +1] +
                                         phi[ijk      +pencil  ] +
                                         phi[ijk+plane         ] ) +
                                 center*(phi[ijk               ] )
                              );
          phi_new[ijk] = phi[ijk] - TwoThirds*lambda[ijk]*(helmholtz-rhs[ijk]);
        }}}
      } // ss-loop
    } // box-loop
    domain->cycles.smooth[level] += (uint64_t)(CycleTime()-_timeStart);
  } // s-loop
}

//------------------------------------------------------------------------------------------------------------------------------
