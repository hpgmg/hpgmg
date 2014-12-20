//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void apply_BCs_linear(level_type * level, int x_id){
  if(level->boundary_condition.type == BC_PERIODIC)return; // no BC's to apply !

  // for cell-centered, we need to fill in the ghost zones to apply any BC's
  // this code does a simple linear interpolation for homogeneous dirichlet
  //
  //   . . . . . . . . .          . . . . . . . . .
  //   .       .       .          .       .       .
  //   .   ?   .   ?   .          .+x(0,0).-x(0,0).
  //   .       .       .          .       .       .
  //   . . . . +-------+          . . . . +-------+
  //   .       |       |          .       |       |
  //   .   ?   | x(0,0)|          .-x(0,0)| x(0,0)|
  //   .       |       |          .       |       |
  //   . . . . +-------+          . . . . +-------+
  //           ^
  //           domain boundary is the face... i.e. between two array indices !!! 
  //

  uint64_t _timeStart = CycleTime();
  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ x      = level->my_boxes[box].vectors[        x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
  //double * __restrict__  valid = level->my_boxes[box].vectors[VECTOR_VALID] + ghosts*(1+jStride+kStride);

    int box_on_low_i  = (level->my_boxes[box].low.i     ==            0);
    int box_on_low_j  = (level->my_boxes[box].low.j     ==            0);
    int box_on_low_k  = (level->my_boxes[box].low.k     ==            0);
    int box_on_high_i = (level->my_boxes[box].low.i+dim == level->dim.i);
    int box_on_high_j = (level->my_boxes[box].low.j+dim == level->dim.j);
    int box_on_high_k = (level->my_boxes[box].low.k+dim == level->dim.k);

    if(level->boundary_condition.type == BC_DIRICHLET){
      int i,j,k,normal;
      double s;

      // note, just because you are in a corner ghost zone, doesn't mean you are on the corner of the domain.
      // thus, one needs to calculate the normal to the domain (not normal to box) in each ghost zone region
      // depending on whether this normal is on a domain face, edge, or corner, one needs to choose 's' appropriately

      // calculate a normal vector for this face                                                              // if face is on a domain boundary, impose the boundary condition using the calculated normal
      s=1;if(box_on_low_i ){normal= 1+      0+      0;s*=-1;}                                                           if(box_on_low_i ){i= -1;j  =0;k  =0;for(j=0;j<dim;j++)for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;if(box_on_low_j ){normal= 0+jStride+      0;s*=-1;}                                                           if(box_on_low_j ){i=  0;j= -1;k  =0;for(k=0;k<dim;k++)for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;if(box_on_low_k ){normal= 0+      0+kStride;s*=-1;}                                                           if(box_on_low_k ){i=  0;j  =0;k= -1;for(j=0;j<dim;j++)for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;if(box_on_high_i){normal=-1+      0+      0;s*=-1;}                                                           if(box_on_high_i){i=dim;j  =0;k  =0;for(j=0;j<dim;j++)for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;if(box_on_high_j){normal= 0-jStride+      0;s*=-1;}                                                           if(box_on_high_j){i=  0;j=dim;k  =0;for(k=0;k<dim;k++)for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;if(box_on_high_k){normal= 0+      0-kStride;s*=-1;}                                                           if(box_on_high_k){i=  0;j  =0;k=dim;for(j=0;j<dim;j++)for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}

      // calculate a normal vector for this edge                                                                                              // if edge is on a domain boundary, impose the boundary condition using the calculated normal
      s=1;normal=0;if(box_on_low_j ){normal+=jStride;s*=-1;}if(box_on_low_k ){normal+=kStride;s*=-1;}                                         if(box_on_low_j ||box_on_low_k ){i=  0;j= -1;k= -1;for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_j){normal-=jStride;s*=-1;}if(box_on_low_k ){normal+=kStride;s*=-1;}                                         if(box_on_high_j||box_on_low_k ){i=  0;j=dim;k= -1;for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_low_j ){normal+=jStride;s*=-1;}if(box_on_high_k){normal-=kStride;s*=-1;}                                         if(box_on_low_j ||box_on_high_k){i=  0;j= -1;k=dim;for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_j){normal-=jStride;s*=-1;}if(box_on_high_k){normal-=kStride;s*=-1;}                                         if(box_on_high_j||box_on_high_k){i=  0;j=dim;k=dim;for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_low_i ){normal+=      1;s*=-1;}if(box_on_low_k ){normal+=kStride;s*=-1;}                                         if(box_on_low_i ||box_on_low_k ){i= -1;j=  0;k= -1;for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_i){normal-=      1;s*=-1;}if(box_on_low_k ){normal+=kStride;s*=-1;}                                         if(box_on_high_i||box_on_low_k ){i=dim;j=  0;k= -1;for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_low_i ){normal+=      1;s*=-1;}if(box_on_high_k){normal-=kStride;s*=-1;}                                         if(box_on_low_i ||box_on_high_k){i= -1;j=  0;k=dim;for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_i){normal-=      1;s*=-1;}if(box_on_high_k){normal-=kStride;s*=-1;}                                         if(box_on_high_i||box_on_high_k){i=dim;j=  0;k=dim;for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_low_i ){normal+=      1;s*=-1;}if(box_on_low_j ){normal+=jStride;s*=-1;}                                         if(box_on_low_i ||box_on_low_j ){i= -1;j= -1;k=  0;for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_i){normal-=      1;s*=-1;}if(box_on_low_j ){normal+=jStride;s*=-1;}                                         if(box_on_high_i||box_on_low_j ){i=dim;j= -1;k=  0;for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_low_i ){normal+=      1;s*=-1;}if(box_on_high_j){normal-=jStride;s*=-1;}                                         if(box_on_low_i ||box_on_high_j){i= -1;j=dim;k=  0;for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_i){normal-=      1;s*=-1;}if(box_on_high_j){normal-=jStride;s*=-1;}                                         if(box_on_high_i||box_on_high_j){i=dim;j=dim;k=  0;for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      
      // calculate a normal vector for this corner                                                                                            // if corner is on a domain boundary, impose the boundary condition using the calculated normal
      s=1;normal=0;if(box_on_low_i ){normal+=      1;s*=-1;}if(box_on_low_j ){normal+=jStride;s*=-1;}if(box_on_low_k ){normal+=kStride;s*=-1;}if(box_on_low_i || box_on_low_j || box_on_low_k ){i= -1;j= -1;k= -1;{int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_i){normal-=      1;s*=-1;}if(box_on_low_j ){normal+=jStride;s*=-1;}if(box_on_low_k ){normal+=kStride;s*=-1;}if(box_on_high_i|| box_on_low_j || box_on_low_k ){i=dim;j= -1;k= -1;{int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_low_i ){normal+=      1;s*=-1;}if(box_on_high_j){normal-=jStride;s*=-1;}if(box_on_low_k ){normal+=kStride;s*=-1;}if(box_on_low_i || box_on_high_j|| box_on_low_k ){i= -1;j=dim;k= -1;{int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_i){normal-=      1;s*=-1;}if(box_on_high_j){normal-=jStride;s*=-1;}if(box_on_low_k ){normal+=kStride;s*=-1;}if(box_on_high_i|| box_on_high_j|| box_on_low_k ){i=dim;j=dim;k= -1;{int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_low_i ){normal+=      1;s*=-1;}if(box_on_low_j ){normal+=jStride;s*=-1;}if(box_on_high_k){normal-=kStride;s*=-1;}if(box_on_low_i || box_on_low_j || box_on_high_k){i= -1;j= -1;k=dim;{int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_i){normal-=      1;s*=-1;}if(box_on_low_j ){normal+=jStride;s*=-1;}if(box_on_high_k){normal-=kStride;s*=-1;}if(box_on_high_i|| box_on_low_j || box_on_high_k){i=dim;j= -1;k=dim;{int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_low_i ){normal+=      1;s*=-1;}if(box_on_high_j){normal-=jStride;s*=-1;}if(box_on_high_k){normal-=kStride;s*=-1;}if(box_on_low_i || box_on_high_j|| box_on_high_k){i= -1;j=dim;k=dim;{int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}
      s=1;normal=0;if(box_on_high_i){normal-=      1;s*=-1;}if(box_on_high_j){normal-=jStride;s*=-1;}if(box_on_high_k){normal-=kStride;s*=-1;}if(box_on_high_i|| box_on_high_j|| box_on_high_k){i=dim;j=dim;k=dim;{int ijk=i+j*jStride+k*kStride;x[ijk]=s*x[ijk+normal];}}

    }
  }
  level->cycles.boundary_conditions += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
  // for cell-centered, we need to fill in the ghost zones to apply any BC's
  // this code does a 4th order scheme for homogeneous dirichlet
  //
  //   .       |       |       |       |       |
  //   . . . . +-------+-------+-------+-------+--
  //   .       |       |       |       |       |
  //   .       | x(0,3)| x(1,3)| x(2,3)| x(3,3)|
  //   .       |       |       |       |       |
  //   . . . . +-------+-------+-------+-------+--
  //   .       |       |       |       |       |
  //   .       | x(0,2)| x(1,2)| x(2,2)| x(3,2)|
  //   .       |       |       |       |       |
  //   . . . . +-------+-------+-------+-------+--
  //   .       |       |       |       |       |
  //   . (-1,1)| x(0,1)| x(1,1)| x(2,1)| x(3,1)|
  //   .       |       |       |       |       |
  //   . . . . +-------+-------+-------+-------+--
  //   .       |       |       |       |       |
  //   .       | x(0,0)| x(1,0)| x(2,0)| x(3,0)|
  //   .       |       |       |       |       |
  //   . . . . +-------+-------+-------+-------+-- <<< domain boundary is the face... i.e. between two array indices !!! 
  //   .       .       .       .       .       .
  //   .(-1,-1).       . (1,-1).       .       .
  //   .       .       .       .       .       .
  //   . . . . . . . . . . . . . . . . . . . . . .
  //

//------------------------------------------------------------------------------------------------------------------------------

  // for cell-centered, we need to fill in the ghost zones to apply any BC's
  // this code does a 2nd order scheme for homogeneous dirichlet
  //
  //   .       |       |       |
  //   . . . . +-------+-------+--
  //   .       |       |       | 
  //   . (-1,1)| x(0,1)| x(1,1)| 
  //   .       |       |       | 
  //   . . . . +-------+-------+--
  //   .       |       |       | 
  //   .       | x(0,0)| x(1,0)| 
  //   .       |       |       | 
  //   . . . . +-------+-------+-- <<< domain boundary is the face... i.e. between two array indices !!! 
  //   .       .       .       .
  //   .(-1,-1).       . (1,-1).
  //   .       .       .       .
  //   . . . . . . . . . . . . . .
  //

//------------------------------------------------------------------------------------------------------------------------------
