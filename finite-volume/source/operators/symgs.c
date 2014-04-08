//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int phi_id, int rhs_id, double a, double b){
  int box,s;
  int ghosts = level->box_ghosts;
  int radius     = STENCIL_RADIUS;
  int starShaped = STENCIL_STAR_SHAPED;
  int communicationAvoiding = ghosts > radius; 

  for(s=0;s<2*NUM_SMOOTHS;s++){ // there are two sweeps (forward/backward) per GS smooth
    exchange_boundary(level,phi_id,starShaped);
            apply_BCs(level,phi_id);

    // now do ghosts communication-avoiding smooths on each box...
    uint64_t _timeStart = CycleTime();
    #pragma omp parallel for private(box) num_threads(level->concurrent_boxes)
    for(box=0;box<level->num_my_boxes;box++){
      int i,j,k;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int     dim = level->my_boxes[box].dim;
      const double h2inv = 1.0/(level->h*level->h);
            double * __restrict__ phi      = level->my_boxes[box].components[        phi_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
      const double * __restrict__ rhs      = level->my_boxes[box].components[        rhs_id] + ghosts*(1+jStride+kStride);
      const double * __restrict__ alpha    = level->my_boxes[box].components[STENCIL_ALPHA ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_i   = level->my_boxes[box].components[STENCIL_BETA_I] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_j   = level->my_boxes[box].components[STENCIL_BETA_J] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_k   = level->my_boxes[box].components[STENCIL_BETA_K] + ghosts*(1+jStride+kStride);
      const double * __restrict__ Dinv     = level->my_boxes[box].components[STENCIL_DINV  ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ valid    = level->my_boxes[box].components[STENCIL_VALID ] + ghosts*(1+jStride+kStride); // cell is inside the domain
          

      if( (s&0x1)==0 ){ // forward sweep
        for(k=0;k<dim;k++){
        for(j=0;j<dim;j++){
        for(i=0;i<dim;i++){
          int ijk = i + j*jStride + k*kStride;
          double Ax = apply_op_ijk(phi);
          phi[ijk] = phi[ijk] + Dinv[ijk]*(rhs[ijk]-Ax);
        }}}
      }else{ // backward sweep
        for(k=dim-1;k>=0;k--){
        for(j=dim-1;j>=0;j--){
        for(i=dim-1;i>=0;i--){
          int ijk = i + j*jStride + k*kStride;
          double Ax = apply_op_ijk(phi);
          phi[ijk] = phi[ijk] + Dinv[ijk]*(rhs[ijk]-Ax);
        }}}
      }

    } // boxes
    level->cycles.smooth += (uint64_t)(CycleTime()-_timeStart);
  } // s-loop
}


//------------------------------------------------------------------------------------------------------------------------------
