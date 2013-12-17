//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include "../timer.h"
//------------------------------------------------------------------------------------------------------------------------------
void smooth(domain_type * domain, int level, int phi_id, int rhs_id, double a, double b){
  int CollaborativeThreadingBoxSize = 100000; // i.e. never
  #ifdef __COLLABORATIVE_THREADING
    CollaborativeThreadingBoxSize = 1 << __COLLABORATIVE_THREADING;
  #endif
  int omp_across_boxes = (domain->subdomains[0].levels[level].dim.i <  CollaborativeThreadingBoxSize);
  int omp_within_a_box = (domain->subdomains[0].levels[level].dim.i >= CollaborativeThreadingBoxSize);
  int box,s;
  int ghosts = domain->ghosts;

  // if communication-avoiding, need RHS for stencils in ghost zones
  if(ghosts>1)exchange_boundary(domain,level,rhs_id,1,1,1); 

  for(s=0;s<numSmooths;s+=ghosts){
    exchange_boundary(domain,level,phi_id,1,ghosts>1,ghosts>1);  // corners/edges if doing communication-avoiding...

    // now do ghosts communication-avoiding smooths on each box...
    uint64_t _timeStart = CycleTime();
    #pragma omp parallel for private(box) if(omp_across_boxes)
    for(box=0;box<domain->subdomains_per_rank;box++){
      int i,j,k,ss;
      int pencil = domain->subdomains[box].levels[level].pencil;
      int  plane = domain->subdomains[box].levels[level].plane;
    //int ghosts = domain->subdomains[box].levels[level].ghosts;
      int  dim_k = domain->subdomains[box].levels[level].dim.k;
      int  dim_j = domain->subdomains[box].levels[level].dim.j;
      int  dim_i = domain->subdomains[box].levels[level].dim.i;
      double h2inv = 1.0/(domain->h[level]*domain->h[level]);
      double * __restrict__ phi      = domain->subdomains[box].levels[level].grids[  phi_id] + ghosts*(1+pencil+plane); // i.e. [0] = first non ghost zone point
      double * __restrict__ phi_new  = domain->subdomains[box].levels[level].grids[  phi_id] + ghosts*(1+pencil+plane); // i.e. [0] = first non ghost zone point
      double * __restrict__ rhs      = domain->subdomains[box].levels[level].grids[  rhs_id] + ghosts*(1+pencil+plane);
      double * __restrict__ alpha    = domain->subdomains[box].levels[level].grids[__alpha ] + ghosts*(1+pencil+plane);
      double * __restrict__ beta_i   = domain->subdomains[box].levels[level].grids[__beta_i] + ghosts*(1+pencil+plane);
      double * __restrict__ beta_j   = domain->subdomains[box].levels[level].grids[__beta_j] + ghosts*(1+pencil+plane);
      double * __restrict__ beta_k   = domain->subdomains[box].levels[level].grids[__beta_k] + ghosts*(1+pencil+plane);
      double * __restrict__ lambda   = domain->subdomains[box].levels[level].grids[__lambda] + ghosts*(1+pencil+plane);
      double * __restrict__ RedBlack[2] = {domain->subdomains[box].levels[level].RedBlack_FP[0] + ghosts*(1+pencil), 
                                           domain->subdomains[box].levels[level].RedBlack_FP[1] + ghosts*(1+pencil)};
          

      int ghostsToOperateOn=ghosts-1;
      for(ss=s;ss<s+ghosts;ss++,ghostsToOperateOn--){
        #if   defined(__GSRB_CONDITIONAL)
        #warning GSRB on every point with conditional assignment for Red-Black
        #pragma omp parallel for private(k,j,i) if(omp_within_a_box) collapse(2)
        for(k=0-ghostsToOperateOn;k<dim_k+ghostsToOperateOn;k++){
        for(j=0-ghostsToOperateOn;j<dim_j+ghostsToOperateOn;j++){
        #pragma simd always
        for(i=0-ghostsToOperateOn;i<dim_i+ghostsToOperateOn;i++){
              int ijk = i + j*pencil + k*plane;
              int doit = ((i^j^k^ss^1)&1);
              double helmholtz =  a*alpha[ijk]*phi[ijk]
                                 -b*h2inv*(
                                    beta_i[ijk+1     ]*( phi[ijk+1     ]-phi[ijk       ] )
                                   -beta_i[ijk       ]*( phi[ijk       ]-phi[ijk-1     ] )
                                   +beta_j[ijk+pencil]*( phi[ijk+pencil]-phi[ijk       ] )
                                   -beta_j[ijk       ]*( phi[ijk       ]-phi[ijk-pencil] )
                                   +beta_k[ijk+plane ]*( phi[ijk+plane ]-phi[ijk       ] )
                                   -beta_k[ijk       ]*( phi[ijk       ]-phi[ijk-plane ] )
                                  );
              phi_new[ijk] = (doit) ? phi[ijk] - lambda[ijk]*(helmholtz-rhs[ijk]) : phi[ijk];
        }}}
        #elif defined(__GSRB_STRIDE2)
        #warning GSRB using stride-2 accesses
        #pragma omp parallel for private(k,j,i) if(omp_within_a_box) collapse(2)
        for(k=0-ghostsToOperateOn;k<dim_k+ghostsToOperateOn;k++){
        for(j=0-ghostsToOperateOn;j<dim_j+ghostsToOperateOn;j++){
        for(i=((j^k^ss^1)&1)+1-ghosts;i<dim_i+ghostsToOperateOn;i+=2){ // stride-2 GSRB
              int ijk = i + j*pencil + k*plane;
              double helmholtz =  a*alpha[ijk]*phi[ijk]
                                 -b*h2inv*(
                                    beta_i[ijk+1     ]*( phi[ijk+1     ]-phi[ijk       ] )
                                   -beta_i[ijk       ]*( phi[ijk       ]-phi[ijk-1     ] )
                                   +beta_j[ijk+pencil]*( phi[ijk+pencil]-phi[ijk       ] )
                                   -beta_j[ijk       ]*( phi[ijk       ]-phi[ijk-pencil] )
                                   +beta_k[ijk+plane ]*( phi[ijk+plane ]-phi[ijk       ] )
                                   -beta_k[ijk       ]*( phi[ijk       ]-phi[ijk-plane ] )
                                  );
              phi_new[ijk] = phi[ijk] - lambda[ijk]*(helmholtz-rhs[ijk]);
        }}}
        #elif defined(__GSRB_FP)
        #warning GSRB using pre-computed 1.0/0.0 FP array for Red-Black
        #pragma omp parallel for private(k,j,i) if(omp_within_a_box) collapse(2)
        for(k=0-ghostsToOperateOn;k<dim_k+ghostsToOperateOn;k++){
        for(j=0-ghostsToOperateOn;j<dim_j+ghostsToOperateOn;j++){int EvenOdd = (k^ss)&1;
        for(i=0-ghostsToOperateOn;i<dim_i+ghostsToOperateOn;i++){
              int ij  = i + j*pencil;
              int ijk = i + j*pencil + k*plane;
              double helmholtz =  a*alpha[ijk]*phi[ijk]
                                 -b*h2inv*(
                                    beta_i[ijk+1     ]*( phi[ijk+1     ]-phi[ijk       ] )
                                   -beta_i[ijk       ]*( phi[ijk       ]-phi[ijk-1     ] )
                                   +beta_j[ijk+pencil]*( phi[ijk+pencil]-phi[ijk       ] )
                                   -beta_j[ijk       ]*( phi[ijk       ]-phi[ijk-pencil] )
                                   +beta_k[ijk+plane ]*( phi[ijk+plane ]-phi[ijk       ] )
                                   -beta_k[ijk       ]*( phi[ijk       ]-phi[ijk-plane ] )
                                  );
              phi_new[ijk] = phi[ijk] - RedBlack[EvenOdd][ij]*lambda[ijk]*(helmholtz-rhs[ijk]); // compiler seems to get confused unless there are disjoint read/write pointers
        }}}
        #else
        #warning GSRB using if-then-else on loop indices for Red-Black
        #pragma omp parallel for private(k,j,i) if(omp_within_a_box) collapse(2)
        for(k=0-ghostsToOperateOn;k<dim_k+ghostsToOperateOn;k++){
        for(j=0-ghostsToOperateOn;j<dim_j+ghostsToOperateOn;j++){
        for(i=0-ghostsToOperateOn;i<dim_i+ghostsToOperateOn;i++){
        if((i^j^k^ss^1)&1){ // looks very clean when [0] is i,j,k=0,0,0 
              int ijk = i + j*pencil + k*plane;
              double helmholtz =  a*alpha[ijk]*phi[ijk]
                                 -b*h2inv*(
                                    beta_i[ijk+1     ]*( phi[ijk+1     ]-phi[ijk       ] )
                                   -beta_i[ijk       ]*( phi[ijk       ]-phi[ijk-1     ] )
                                   +beta_j[ijk+pencil]*( phi[ijk+pencil]-phi[ijk       ] )
                                   -beta_j[ijk       ]*( phi[ijk       ]-phi[ijk-pencil] )
                                   +beta_k[ijk+plane ]*( phi[ijk+plane ]-phi[ijk       ] )
                                   -beta_k[ijk       ]*( phi[ijk       ]-phi[ijk-plane ] )
                                  );
              phi_new[ijk] = phi[ijk] - lambda[ijk]*(helmholtz-rhs[ijk]);
        }}}}
        #endif
      } // ss-loop
    } // boxes
    domain->cycles.smooth[level] += (uint64_t)(CycleTime()-_timeStart);
  } // s-loop
}


//------------------------------------------------------------------------------------------------------------------------------
