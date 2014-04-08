//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// calculate res_id = rhs_id - A(x_id)

void residual(level_type * level, int res_id, int x_id, int rhs_id, double a, double b){
  int starShaped = STENCIL_STAR_SHAPED;

  // exchange the boundary for x in prep for Ax...
  exchange_boundary(level,x_id,starShaped);
          apply_BCs(level,x_id);

  // now do residual/restriction proper...
  uint64_t _timeStart = CycleTime();
  int box;

  #pragma omp parallel for private(box) num_threads(level->concurrent_boxes)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    int jStride = level->my_boxes[box].jStride;
    int kStride = level->my_boxes[box].kStride;
    int  ghosts = level->my_boxes[box].ghosts;
    int     dim = level->my_boxes[box].dim;
    double h2inv = 1.0/(level->h*level->h);
    const double * __restrict__ x      = level->my_boxes[box].components[          x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    const double * __restrict__ rhs    = level->my_boxes[box].components[        rhs_id] + ghosts*(1+jStride+kStride);
    const double * __restrict__ alpha  = level->my_boxes[box].components[STENCIL_ALPHA ] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_i = level->my_boxes[box].components[STENCIL_BETA_I] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_j = level->my_boxes[box].components[STENCIL_BETA_J] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_k = level->my_boxes[box].components[STENCIL_BETA_K] + ghosts*(1+jStride+kStride);
    const double * __restrict__ valid  = level->my_boxes[box].components[STENCIL_VALID ] + ghosts*(1+jStride+kStride); // cell is inside the domain
          double * __restrict__ res    = level->my_boxes[box].components[        res_id] + ghosts*(1+jStride+kStride);

    #pragma omp parallel for private(k,j,i) num_threads(level->threads_per_box) OMP_COLLAPSE
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
      int ijk = i + j*jStride + k*kStride;
      double Ax = apply_op_ijk(x);
      res[ijk] = rhs[ijk]-Ax;
    }}}
  }
  level->cycles.residual += (uint64_t)(CycleTime()-_timeStart);
}

