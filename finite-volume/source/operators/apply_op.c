//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void apply_op(level_type * level, int Ax_id, int x_id, double a, double b){  // y=Ax
  int starShaped = __STENCIL_STAR_SHAPED;

  // exchange the boundary of x in preparation for Ax
  exchange_boundary(level,x_id,starShaped);
          apply_BCs(level,x_id);

  // now do Ax proper...
  uint64_t _timeStart = CycleTime();
  int box;

  #pragma omp parallel for private(box) num_threads(level->concurrent_boxes)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k,s;
    int jStride = level->my_boxes[box].jStride;
    int kStride = level->my_boxes[box].kStride;
    int  ghosts = level->my_boxes[box].ghosts;
    int     dim = level->my_boxes[box].dim;
    double h2inv = 1.0/(level->h*level->h);
    const double * __restrict__ x      = level->my_boxes[box].components[     x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
          double * __restrict__ Ax     = level->my_boxes[box].components[    Ax_id] + ghosts*(1+jStride+kStride); 
    const double * __restrict__ alpha  = level->my_boxes[box].components[ __alpha ] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_i = level->my_boxes[box].components[ __beta_i] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_j = level->my_boxes[box].components[ __beta_j] + ghosts*(1+jStride+kStride);
    const double * __restrict__ beta_k = level->my_boxes[box].components[ __beta_k] + ghosts*(1+jStride+kStride);
    const double * __restrict__  valid = level->my_boxes[box].components[  __valid] + ghosts*(1+jStride+kStride);

    #pragma omp parallel for private(k,j,i) num_threads(level->threads_per_box) __OMP_COLLAPSE
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
      int ijk = i + j*jStride + k*kStride;
      Ax[ijk] = __apply_op(x);
    }}}
  }
  level->cycles.apply_op += (uint64_t)(CycleTime()-_timeStart);
}
//------------------------------------------------------------------------------------------------------------------------------
