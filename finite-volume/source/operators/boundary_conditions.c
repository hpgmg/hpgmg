//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void apply_BCs_linear(level_type * level, int x_id){
  if(level->domain_boundary_condition == BC_PERIODIC)return; // no BC's to apply !

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
  int omp_across_boxes = 1;
  int omp_within_a_box = 0;
  int box;

  #pragma omp parallel for private(box) OMP_THREAD_ACROSS_BOXES(level->concurrent_boxes)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k,s;
    int jStride = level->my_boxes[box].jStride;
    int kStride = level->my_boxes[box].kStride;
    int  ghosts = level->my_boxes[box].ghosts;
    int     dim = level->my_boxes[box].dim;
    double h2inv = 1.0/(level->h*level->h);
    double * __restrict__ x      = level->my_boxes[box].vectors[        x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double * __restrict__  valid = level->my_boxes[box].vectors[VECTOR_VALID] + ghosts*(1+jStride+kStride);

    if(level->domain_boundary_condition == BC_DIRICHLET){
      // why these and not -1, -5, +77 ???
                  k= -1;if((level->my_boxes[box].low.k     ==            0))                                                                for(j=0;j<dim;j++)for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk          +kStride];} // face
            j= -1;      if((level->my_boxes[box].low.j     ==            0))                                                                for(k=0;k<dim;k++)for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk  +jStride        ];} // face
      i= -1;            if((level->my_boxes[box].low.i     ==            0))                                                                for(k=0;k<dim;k++)for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk+1                ];} // face
      i=dim;            if((level->my_boxes[box].low.i+dim == level->dim.i))                                                                for(k=0;k<dim;k++)for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk-1                ];} // face
            j=dim;      if((level->my_boxes[box].low.j+dim == level->dim.j))                                                                for(k=0;k<dim;k++)for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk  -jStride        ];} // face
                  k=dim;if((level->my_boxes[box].low.k+dim == level->dim.k))                                                                for(j=0;j<dim;j++)for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk          -kStride];} // face

            j= -1;k= -1;if((level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0))                                for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk  +jStride+kStride];} // edge
      i= -1;      k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.k     ==            0))                                for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk+1        +kStride];} // edge
      i=dim;      k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.k     ==            0))                                for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk-1        +kStride];} // edge
            j=dim;k= -1;if((level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0))                                for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk  -jStride+kStride];} // edge
      i= -1;j= -1;      if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0))                                for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk+1+jStride        ];} // edge
      i=dim;j= -1;      if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0))                                for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk-1+jStride        ];} // edge
      i= -1;j=dim;      if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j))                                for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk+1-jStride        ];} // edge
      i=dim;j=dim;      if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j))                                for(k=0;k<dim;k++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk-1-jStride        ];} // edge
            j= -1;k=dim;if((level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k))                                for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk  +jStride-kStride];} // edge
      i= -1;      k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k))                                for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk+1        -kStride];} // edge
      i=dim;      k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.k+dim == level->dim.k))                                for(j=0;j<dim;j++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk-1        -kStride];} // edge
            j=dim;k=dim;if((level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k))                                for(i=0;i<dim;i++){int ijk=i+j*jStride+k*kStride;x[ijk]=x[ijk  -jStride-kStride];} // edge

      i= -1;j= -1;k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0)){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk+1+jStride+kStride];} // corner
      i=dim;j= -1;k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0)){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk-1+jStride+kStride];} // corner
      i= -1;j=dim;k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0)){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk+1-jStride+kStride];} // corner
      i=dim;j=dim;k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0)){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk-1-jStride+kStride];} // corner
      i= -1;j= -1;k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k)){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk+1+jStride-kStride];} // corner
      i=dim;j= -1;k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k)){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk-1+jStride-kStride];} // corner
      i= -1;j=dim;k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k)){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk+1-jStride-kStride];} // corner
      i=dim;j=dim;k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k)){int ijk=i+j*jStride+k*kStride;x[ijk]=-x[ijk-1-jStride-kStride];} // corner
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
