//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include "../timer.h"
//------------------------------------------------------------------------------------------------------------------------------
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
//------------------------------------------------------------------------------------------------------------------------------
// better solution would be to adapt the box size as the problem shrinks...
// i.e. fix unit stride at 4KB and calculate BlockJ = ((STANZA+dim.i-1)/dim.i)
// similarly, fix BlockK to get some reuse and have enough tasks...
//------------------------------------------------------------------------------------------------------------------------------
// Kludge for now...  
#define BlockJ 16
#define BlockK 4
//------------------------------------------------------------------------------------------------------------------------------
void __box_smooth_GSRB_multiple(box_type *box, int phi_id, int rhs_id, double a, double b, int s){
  int jj,kk;
  int pencil = box->pencil;
  int plane = box->plane;
  int ghosts = box->ghosts;
  double h2inv = 1.0/(box->h*box->h);
  double * __restrict__ phi    = box->grids[  phi_id] + ghosts*plane + ghosts*pencil + ghosts; // i.e. [0] = first non ghost zone point
  double * __restrict__ phi_new= box->grids[  phi_id] + ghosts*plane + ghosts*pencil + ghosts;
  double * __restrict__ rhs    = box->grids[  rhs_id] + ghosts*plane + ghosts*pencil + ghosts;
  double * __restrict__ alpha  = box->grids[__alpha ] + ghosts*plane + ghosts*pencil + ghosts;
  double * __restrict__ beta_i = box->grids[__beta_i] + ghosts*plane + ghosts*pencil + ghosts;
  double * __restrict__ beta_j = box->grids[__beta_j] + ghosts*plane + ghosts*pencil + ghosts;
  double * __restrict__ beta_k = box->grids[__beta_k] + ghosts*plane + ghosts*pencil + ghosts;
  double * __restrict__ lambda = box->grids[__lambda] + ghosts*plane + ghosts*pencil + ghosts;
  int ghostsToOperateOn=ghosts-1;
  int ss;
  int big_box=0;

  // don't subdivide small boxes into tasks (too much overhead from omp task...)
  if(box->dim.k>8)big_box=1;
  if(box->dim.j>8)big_box=1;

  // do ghosts iterations on this list of tasks...
  for(ss=s;ss<s+ghosts;ss++,ghostsToOperateOn--){

    // iterate through all cache blocks within this box and queue a task...
    for(kk=0-ghostsToOperateOn;kk<box->dim.k+ghostsToOperateOn;kk+=BlockK){
    for(jj=0-ghostsToOperateOn;jj<box->dim.j+ghostsToOperateOn;jj+=BlockJ){
    #pragma omp task if(big_box)
    {
      int i,j,k;
      int highK,highJ;
      highK = MIN(kk+BlockK,box->dim.k+ghostsToOperateOn);
      highJ = MIN(jj+BlockJ,box->dim.j+ghostsToOperateOn);
      #if   defined(__GSRB_CONDITIONAL)
      #warning GSRB on every point with conditional assignment for Red-Black
      for(k=kk;k<highK;k++){
      for(j=jj;j<highJ;j++){
      #pragma simd always
      for(i=0-ghostsToOperateOn;i<box->dim.i+ghostsToOperateOn;i++){
        int ijk = i + j*pencil + k*plane;
        int doit = ((i^(j^k^ss^1))&1);
        double helmholtz =  a*alpha[ijk]*phi[ijk]
                           -b*h2inv*(
                              beta_i[ijk+1     ]*( phi[ijk+1     ]-phi[ijk       ] )
                             -beta_i[ijk       ]*( phi[ijk       ]-phi[ijk-1     ] )
                             +beta_j[ijk+pencil]*( phi[ijk+pencil]-phi[ijk       ] )
                             -beta_j[ijk       ]*( phi[ijk       ]-phi[ijk-pencil] )
                             +beta_k[ijk+plane ]*( phi[ijk+plane ]-phi[ijk       ] )
                             -beta_k[ijk       ]*( phi[ijk       ]-phi[ijk-plane ] )
                            );
        //double delta = doit ? lambda[ijk]*(helmholtz-rhs[ijk]) : 0.0;
        //phi_new[ijk] = phi[ijk] - delta;
        phi_new[ijk] = doit ? phi[ijk] - lambda[ijk]*(helmholtz-rhs[ijk]) : phi[ijk];
      }}}
      #elif defined(__GSRB_STRIDE2)
      #warning GSRB using stride-2 accesses
      for(k=kk;k<highK;k++){
      for(j=jj;j<highJ;j++){
      for(i=((j^k^ss^1)&1)+1-ghosts;i<box->dim.i+ghostsToOperateOn;i+=2){ // stride-2 GSRB
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
      for(k=kk;k<highK;k++){int EvenOdd = (k^ss)&1;
      for(j=jj;j<highJ;j++){
      for(i=0-ghostsToOperateOn;i<box->dim.i+ghostsToOperateOn;i++){
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
      for(k=kk;k<highK;k++){
      for(j=jj;j<highJ;j++){
      for(i=0-ghostsToOperateOn;i<box->dim.i+ghostsToOperateOn;i++){
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
    }}}

    // If doing communication avoiding, we dependent tasks cannot get too far ahead.  
    // As I have no idea how to perform p2p synchronization among omp tasks, I'll just barrier...
    if(ghostsToOperateOn>0){
    #pragma omp taskwait
    }
  } // ss
}


//------------------------------------------------------------------------------------------------------------------------------
void smooth(domain_type * domain, int level, int phi_id, int rhs_id, double a, double b){
  int box,s;
  int ghosts = domain->ghosts;
  // if communication-avoiding, need RHS for stencils in ghost zones
  if(ghosts>1)exchange_boundary(domain,level,rhs_id,1,1,1);

  for(s=0;s<numSmooths;s+=ghosts){
    exchange_boundary(domain,level,phi_id,1,ghosts>1,ghosts>1);  // corners/edges if doing communication-avoiding...
    uint64_t _timeStart = CycleTime();
    #pragma omp parallel
    {
      int box;
      #pragma omp for private(box) nowait // <<< needs to be omp for rather than single in order to get enough task injection.  <<< needs to be no wait to ensure idle cores can grab tasks asap
      for(box=0;box<domain->subdomains_per_rank;box++){
        __box_smooth_GSRB_multiple(&domain->subdomains[box].levels[level],phi_id,rhs_id,a,b,s);
      }
    }
    domain->cycles.smooth[level] += (uint64_t)(CycleTime()-_timeStart);
  }
}
