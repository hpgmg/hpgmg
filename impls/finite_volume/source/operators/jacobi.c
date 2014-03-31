//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include "../timer.h"
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int x_id, int rhs_id, double a, double b){
  if(__NUM_SMOOTHS&1){
    printf("error - __NUM_SMOOTHS must be even...\n");
    exit(0);
  }


  int box,s;
  int ghosts = level->box_ghosts;
  int radius     = __STENCIL_RADIUS;
  int starShaped = __STENCIL_STAR_SHAPED;
  int communicationAvoiding = ghosts > radius; 
 
  #ifdef __USE_L1JACOBI
  double weight = 1.0;
  #else
  double weight = 2.0/3.0;
  #endif
 
 
  // if communication-avoiding, need updated RHS for stencils in ghost zones
  if(communicationAvoiding)exchange_boundary(level,rhs_id,0); 

  for(s=0;s<__NUM_SMOOTHS;s+=ghosts){
    // Jacobi ping pongs between x_id and __temp
    if((s&1)==0){exchange_boundary(level,  x_id,starShaped && !communicationAvoiding);apply_BCs(level,  x_id);}
            else{exchange_boundary(level,__temp,starShaped && !communicationAvoiding);apply_BCs(level,__temp);}

    // now do ghosts communication-avoiding smooths on each box...
    uint64_t _timeStart = CycleTime();

    #pragma omp parallel for private(box) num_threads(level->concurrent_boxes)
    for(box=0;box<level->num_my_boxes;box++){
      int i,j,k,ss;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int     dim = level->my_boxes[box].dim;
      const double h2inv = 1.0/(level->h*level->h);
      const double * __restrict__ rhs    = level->my_boxes[box].components[  rhs_id] + ghosts*(1+jStride+kStride);
      const double * __restrict__ alpha  = level->my_boxes[box].components[__alpha ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_i = level->my_boxes[box].components[__beta_i] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_j = level->my_boxes[box].components[__beta_j] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_k = level->my_boxes[box].components[__beta_k] + ghosts*(1+jStride+kStride);
      const double * __restrict__ valid  = level->my_boxes[box].components[__valid ] + ghosts*(1+jStride+kStride); // cell is inside the domain
      #ifdef __USE_L1JACOBI
      const double * __restrict__ lambda = level->my_boxes[box].components[__L1inv ] + ghosts*(1+jStride+kStride);
      #else
      const double * __restrict__ lambda = level->my_boxes[box].components[__Dinv  ] + ghosts*(1+jStride+kStride);
      #endif
      int ghostsToOperateOn=ghosts-1;
      for(ss=s;ss<s+ghosts;ss++,ghostsToOperateOn--){
        const double * __restrict__ x_n;
              double * __restrict__ x_np1;
              if((ss&1)==0){x_n   = level->my_boxes[box].components[  x_id] + ghosts*(1+jStride+kStride);
                            x_np1 = level->my_boxes[box].components[__temp] + ghosts*(1+jStride+kStride);}
                       else{x_n   = level->my_boxes[box].components[__temp] + ghosts*(1+jStride+kStride);
                            x_np1 = level->my_boxes[box].components[  x_id] + ghosts*(1+jStride+kStride);}
        #pragma omp parallel for private(k,j,i) num_threads(level->threads_per_box) __OMP_COLLAPSE
        for(k=0-ghostsToOperateOn;k<dim+ghostsToOperateOn;k++){
        for(j=0-ghostsToOperateOn;j<dim+ghostsToOperateOn;j++){
        for(i=0-ghostsToOperateOn;i<dim+ghostsToOperateOn;i++){
          int ijk = i + j*jStride + k*kStride;
          double Ax_n = __apply_op(x_n);
          x_np1[ijk] = x_n[ijk] + weight*lambda[ijk]*(rhs[ijk]-Ax_n);
        }}}
      } // ss-loop
    } // box-loop
    level->cycles.smooth += (uint64_t)(CycleTime()-_timeStart);
  } // s-loop
}

//------------------------------------------------------------------------------------------------------------------------------
