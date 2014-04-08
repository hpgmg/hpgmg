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

  #pragma omp parallel for private(box) num_threads(level->concurrent_boxes)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k,s;
    int jStride = level->my_boxes[box].jStride;
    int kStride = level->my_boxes[box].kStride;
    int  ghosts = level->my_boxes[box].ghosts;
    int     dim = level->my_boxes[box].dim;
    double h2inv = 1.0/(level->h*level->h);
    double * __restrict__ x      = level->my_boxes[box].components[     x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double * __restrict__  valid = level->my_boxes[box].components[  STENCIL_VALID] + ghosts*(1+jStride+kStride);

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
static inline void BC_4thOrder_face(double * __restrict__ x, int ijk, int aStride){
  x[ijk] = ( -140.0*x[ijk+aStride] + 70.0*x[ijk+2*aStride] - 28.0*x[ijk+3*aStride] + 5.0*x[ijk+4*aStride] ) / 35.0;
}

static inline void BC_4thOrder_edge(double * __restrict__ x, int ijk, int aStride, int bStride){
  x[ijk] = ( 5.0*x[ijk+aStride] - 10.0*x[ijk+2*aStride] + 10.0*x[ijk+3*aStride] - 5.0*x[ijk+4*aStride] + x[ijk+5*aStride] +
             5.0*x[ijk+bStride] - 10.0*x[ijk+2*bStride] + 10.0*x[ijk+3*bStride] - 5.0*x[ijk+4*bStride] + x[ijk+5*bStride] ) * 0.5;
}

static inline void BC_4thOrder_corner(double * __restrict__ x, int ijk, int aStride, int bStride, int cStride){
  x[ijk] = ( 5.0*x[ijk+aStride] - 10.0*x[ijk+2*aStride] + 10.0*x[ijk+3*aStride] - 5.0*x[ijk+4*aStride] + x[ijk+5*aStride] +
             5.0*x[ijk+bStride] - 10.0*x[ijk+2*bStride] + 10.0*x[ijk+3*bStride] - 5.0*x[ijk+4*bStride] + x[ijk+5*bStride] +
             5.0*x[ijk+cStride] - 10.0*x[ijk+2*cStride] + 10.0*x[ijk+3*cStride] - 5.0*x[ijk+4*cStride] + x[ijk+5*cStride] ) / 3.0;
}
//------------------------------------------------------------------------------------------------------------------------------
void apply_BCs_4thOrder(level_type * level, int x_id){
  if(level->domain_boundary_condition == BC_PERIODIC)return; // no BC's to apply !
  if(level->box_dim<4){apply_BCs_2ndOrder(level,x_id);return;}

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

  uint64_t _timeStart = CycleTime();
  int omp_across_boxes = 1;
  int omp_within_a_box = 0;
  int box;
  double OneThirtyFifth = 1.0/35.0;

  #pragma omp parallel for private(box) num_threads(level->concurrent_boxes)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k,s;
    int jStride = level->my_boxes[box].jStride;int jStride2 = jStride*2;int jStride3 = jStride*3;int jStride4 = jStride*4;
    int kStride = level->my_boxes[box].kStride;int kStride2 = kStride*2;int kStride3 = kStride*3;int kStride4 = kStride*4;
    int  ghosts = level->my_boxes[box].ghosts;
    int     dim = level->my_boxes[box].dim;
    double h2inv = 1.0/(level->h*level->h);
    double * __restrict__ x      = level->my_boxes[box].components[     x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double * __restrict__  valid = level->my_boxes[box].components[  STENCIL_VALID] + ghosts*(1+jStride+kStride);

    if(level->domain_boundary_condition == BC_DIRICHLET){
                  k= -1;if((level->my_boxes[box].low.k     ==            0))for(j=0;j<dim;j++)for(i=0;i<dim;i++){BC_4thOrder_face(x,i+j*jStride+k*kStride, kStride);} // face
            j= -1;      if((level->my_boxes[box].low.j     ==            0))for(k=0;k<dim;k++)for(i=0;i<dim;i++){BC_4thOrder_face(x,i+j*jStride+k*kStride, jStride);} // face
      i= -1;            if((level->my_boxes[box].low.i     ==            0))for(k=0;k<dim;k++)for(j=0;j<dim;j++){BC_4thOrder_face(x,i+j*jStride+k*kStride, 1      );} // face
      i=dim;            if((level->my_boxes[box].low.i+dim == level->dim.i))for(k=0;k<dim;k++)for(j=0;j<dim;j++){BC_4thOrder_face(x,i+j*jStride+k*kStride,-1      );} // face
            j=dim;      if((level->my_boxes[box].low.j+dim == level->dim.j))for(k=0;k<dim;k++)for(i=0;i<dim;i++){BC_4thOrder_face(x,i+j*jStride+k*kStride,-jStride);} // face
                  k=dim;if((level->my_boxes[box].low.k+dim == level->dim.k))for(j=0;j<dim;j++)for(i=0;i<dim;i++){BC_4thOrder_face(x,i+j*jStride+k*kStride,-kStride);} // face

      #if 0
            j= -1;k= -1;if((level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0))for(i=0;i<dim;i++){BC_4thOrder_edge(x,i+j*jStride+k*kStride, jStride, kStride);} // edge
      i= -1;      k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.k     ==            0))for(j=0;j<dim;j++){BC_4thOrder_edge(x,i+j*jStride+k*kStride, 1      , kStride);} // edge
      i=dim;      k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.k     ==            0))for(j=0;j<dim;j++){BC_4thOrder_edge(x,i+j*jStride+k*kStride,-1      , kStride);} // edge
            j=dim;k= -1;if((level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0))for(i=0;i<dim;i++){BC_4thOrder_edge(x,i+j*jStride+k*kStride,-jStride, kStride);} // edge
      i= -1;j= -1;      if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0))for(k=0;k<dim;k++){BC_4thOrder_edge(x,i+j*jStride+k*kStride, 1      , jStride);} // edge
      i=dim;j= -1;      if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0))for(k=0;k<dim;k++){BC_4thOrder_edge(x,i+j*jStride+k*kStride,-1      , jStride);} // edge
      i= -1;j=dim;      if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j))for(k=0;k<dim;k++){BC_4thOrder_edge(x,i+j*jStride+k*kStride, 1      ,-jStride);} // edge
      i=dim;j=dim;      if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j))for(k=0;k<dim;k++){BC_4thOrder_edge(x,i+j*jStride+k*kStride,-1      ,-jStride);} // edge
            j= -1;k=dim;if((level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(i=0;i<dim;i++){BC_4thOrder_edge(x,i+j*jStride+k*kStride, jStride,-kStride);} // edge
      i= -1;      k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(j=0;j<dim;j++){BC_4thOrder_edge(x,i+j*jStride+k*kStride, 1      ,-kStride);} // edge
      i=dim;      k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(j=0;j<dim;j++){BC_4thOrder_edge(x,i+j*jStride+k*kStride,-1      ,-kStride);} // edge
            j=dim;k=dim;if((level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(i=0;i<dim;i++){BC_4thOrder_edge(x,i+j*jStride+k*kStride,-jStride,-kStride);} // edge

      i= -1;j= -1;k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0)){BC_4thOrder_corner(x,i+j*jStride+k*kStride, 1, jStride, kStride);} // corner
      i=dim;j= -1;k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0)){BC_4thOrder_corner(x,i+j*jStride+k*kStride,-1, jStride, kStride);} // corner
      i= -1;j=dim;k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0)){BC_4thOrder_corner(x,i+j*jStride+k*kStride, 1,-jStride, kStride);} // corner
      i=dim;j=dim;k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0)){BC_4thOrder_corner(x,i+j*jStride+k*kStride,-1,-jStride, kStride);} // corner
      i= -1;j= -1;k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_4thOrder_corner(x,i+j*jStride+k*kStride, 1, jStride,-kStride);} // corner
      i=dim;j= -1;k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_4thOrder_corner(x,i+j*jStride+k*kStride,-1, jStride,-kStride);} // corner
      i= -1;j=dim;k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_4thOrder_corner(x,i+j*jStride+k*kStride, 1,-jStride,-kStride);} // corner
      i=dim;j=dim;k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_4thOrder_corner(x,i+j*jStride+k*kStride,-1,-jStride,-kStride);} // corner
      #else
            j= -1;k= -1;if((level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0))for(i=0;i<dim;i++){BC_4thOrder_face(x,i+j*jStride+k*kStride, jStride+kStride);} // edge
      i= -1;      k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.k     ==            0))for(j=0;j<dim;j++){BC_4thOrder_face(x,i+j*jStride+k*kStride, 1      +kStride);} // edge
      i=dim;      k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.k     ==            0))for(j=0;j<dim;j++){BC_4thOrder_face(x,i+j*jStride+k*kStride,-1      +kStride);} // edge
            j=dim;k= -1;if((level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0))for(i=0;i<dim;i++){BC_4thOrder_face(x,i+j*jStride+k*kStride,-jStride+kStride);} // edge
      i= -1;j= -1;      if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0))for(k=0;k<dim;k++){BC_4thOrder_face(x,i+j*jStride+k*kStride, 1      +jStride);} // edge
      i=dim;j= -1;      if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0))for(k=0;k<dim;k++){BC_4thOrder_face(x,i+j*jStride+k*kStride,-1      +jStride);} // edge
      i= -1;j=dim;      if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j))for(k=0;k<dim;k++){BC_4thOrder_face(x,i+j*jStride+k*kStride, 1      -jStride);} // edge
      i=dim;j=dim;      if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j))for(k=0;k<dim;k++){BC_4thOrder_face(x,i+j*jStride+k*kStride,-1      -jStride);} // edge
            j= -1;k=dim;if((level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(i=0;i<dim;i++){BC_4thOrder_face(x,i+j*jStride+k*kStride, jStride-kStride);} // edge
      i= -1;      k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(j=0;j<dim;j++){BC_4thOrder_face(x,i+j*jStride+k*kStride, 1      -kStride);} // edge
      i=dim;      k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(j=0;j<dim;j++){BC_4thOrder_face(x,i+j*jStride+k*kStride,-1      -kStride);} // edge
            j=dim;k=dim;if((level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(i=0;i<dim;i++){BC_4thOrder_face(x,i+j*jStride+k*kStride,-jStride-kStride);} // edge

      i= -1;j= -1;k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0)){BC_4thOrder_face(x,i+j*jStride+k*kStride, 1+jStride+kStride);} // corner
      i=dim;j= -1;k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0)){BC_4thOrder_face(x,i+j*jStride+k*kStride,-1+jStride+kStride);} // corner
      i= -1;j=dim;k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0)){BC_4thOrder_face(x,i+j*jStride+k*kStride, 1-jStride+kStride);} // corner
      i=dim;j=dim;k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0)){BC_4thOrder_face(x,i+j*jStride+k*kStride,-1-jStride+kStride);} // corner
      i= -1;j= -1;k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_4thOrder_face(x,i+j*jStride+k*kStride, 1+jStride-kStride);} // corner
      i=dim;j= -1;k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_4thOrder_face(x,i+j*jStride+k*kStride,-1+jStride-kStride);} // corner
      i= -1;j=dim;k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_4thOrder_face(x,i+j*jStride+k*kStride, 1-jStride-kStride);} // corner
      i=dim;j=dim;k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_4thOrder_face(x,i+j*jStride+k*kStride,-1-jStride-kStride);} // corner
      #endif
    }
  }
  level->cycles.boundary_conditions += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
static inline void BC_2ndOrder(double * __restrict__ x, int ijk, int aStride){
  x[ijk] = ( -6.0*x[ijk+aStride] + x[ijk+2*aStride] ) / 3.0;
}

//------------------------------------------------------------------------------------------------------------------------------
void apply_BCs_2ndOrder(level_type * level, int x_id){
  if(level->domain_boundary_condition == BC_PERIODIC)return; // no BC's to apply !
  if(level->box_dim<2){apply_BCs_linear(level,x_id);return;}

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

  uint64_t _timeStart = CycleTime();
  int omp_across_boxes = 1;
  int omp_within_a_box = 0;
  int box;

  #pragma omp parallel for private(box) num_threads(level->concurrent_boxes)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k,s;
    int jStride = level->my_boxes[box].jStride;int jStride2 = jStride*2;int jStride3 = jStride*3;int jStride4 = jStride*4;
    int kStride = level->my_boxes[box].kStride;int kStride2 = kStride*2;int kStride3 = kStride*3;int kStride4 = kStride*4;
    int  ghosts = level->my_boxes[box].ghosts;
    int     dim = level->my_boxes[box].dim;
    double h2inv = 1.0/(level->h*level->h);
    double * __restrict__ x      = level->my_boxes[box].components[           x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double * __restrict__  valid = level->my_boxes[box].components[  STENCIL_VALID] + ghosts*(1+jStride+kStride);

    if(level->domain_boundary_condition == BC_DIRICHLET){
                  k= -1;if((level->my_boxes[box].low.k     ==            0))for(j=0;j<dim;j++)for(i=0;i<dim;i++){BC_2ndOrder(x,i+j*jStride+k*kStride, kStride);} // face
            j= -1;      if((level->my_boxes[box].low.j     ==            0))for(k=0;k<dim;k++)for(i=0;i<dim;i++){BC_2ndOrder(x,i+j*jStride+k*kStride, jStride);} // face
      i= -1;            if((level->my_boxes[box].low.i     ==            0))for(k=0;k<dim;k++)for(j=0;j<dim;j++){BC_2ndOrder(x,i+j*jStride+k*kStride, 1      );} // face
      i=dim;            if((level->my_boxes[box].low.i+dim == level->dim.i))for(k=0;k<dim;k++)for(j=0;j<dim;j++){BC_2ndOrder(x,i+j*jStride+k*kStride,-1      );} // face
            j=dim;      if((level->my_boxes[box].low.j+dim == level->dim.j))for(k=0;k<dim;k++)for(i=0;i<dim;i++){BC_2ndOrder(x,i+j*jStride+k*kStride,-jStride);} // face
                  k=dim;if((level->my_boxes[box].low.k+dim == level->dim.k))for(j=0;j<dim;j++)for(i=0;i<dim;i++){BC_2ndOrder(x,i+j*jStride+k*kStride,-kStride);} // face

            j= -1;k= -1;if((level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0))for(i=0;i<dim;i++){BC_2ndOrder(x,i+j*jStride+k*kStride, jStride+kStride);} // edge
      i= -1;      k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.k     ==            0))for(j=0;j<dim;j++){BC_2ndOrder(x,i+j*jStride+k*kStride, 1      +kStride);} // edge
      i=dim;      k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.k     ==            0))for(j=0;j<dim;j++){BC_2ndOrder(x,i+j*jStride+k*kStride,-1      +kStride);} // edge
            j=dim;k= -1;if((level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0))for(i=0;i<dim;i++){BC_2ndOrder(x,i+j*jStride+k*kStride,-jStride+kStride);} // edge
      i= -1;j= -1;      if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0))for(k=0;k<dim;k++){BC_2ndOrder(x,i+j*jStride+k*kStride, 1      +jStride);} // edge
      i=dim;j= -1;      if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0))for(k=0;k<dim;k++){BC_2ndOrder(x,i+j*jStride+k*kStride,-1      +jStride);} // edge
      i= -1;j=dim;      if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j))for(k=0;k<dim;k++){BC_2ndOrder(x,i+j*jStride+k*kStride, 1      -jStride);} // edge
      i=dim;j=dim;      if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j))for(k=0;k<dim;k++){BC_2ndOrder(x,i+j*jStride+k*kStride,-1      -jStride);} // edge
            j= -1;k=dim;if((level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(i=0;i<dim;i++){BC_2ndOrder(x,i+j*jStride+k*kStride, jStride-kStride);} // edge
      i= -1;      k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(j=0;j<dim;j++){BC_2ndOrder(x,i+j*jStride+k*kStride, 1      -kStride);} // edge
      i=dim;      k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(j=0;j<dim;j++){BC_2ndOrder(x,i+j*jStride+k*kStride,-1      -kStride);} // edge
            j=dim;k=dim;if((level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k))for(i=0;i<dim;i++){BC_2ndOrder(x,i+j*jStride+k*kStride,-jStride-kStride);} // edge

      i= -1;j= -1;k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0)){BC_2ndOrder(x,i+j*jStride+k*kStride, 1+jStride+kStride);} // corner
      i=dim;j= -1;k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k     ==            0)){BC_2ndOrder(x,i+j*jStride+k*kStride,-1+jStride+kStride);} // corner
      i= -1;j=dim;k= -1;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0)){BC_2ndOrder(x,i+j*jStride+k*kStride, 1-jStride+kStride);} // corner
      i=dim;j=dim;k= -1;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k     ==            0)){BC_2ndOrder(x,i+j*jStride+k*kStride,-1-jStride+kStride);} // corner
      i= -1;j= -1;k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_2ndOrder(x,i+j*jStride+k*kStride, 1+jStride-kStride);} // corner
      i=dim;j= -1;k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j     ==            0)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_2ndOrder(x,i+j*jStride+k*kStride,-1+jStride-kStride);} // corner
      i= -1;j=dim;k=dim;if((level->my_boxes[box].low.i     ==            0)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_2ndOrder(x,i+j*jStride+k*kStride, 1-jStride-kStride);} // corner
      i=dim;j=dim;k=dim;if((level->my_boxes[box].low.i+dim == level->dim.i)&&(level->my_boxes[box].low.j+dim == level->dim.j)&&(level->my_boxes[box].low.k+dim == level->dim.k)){BC_2ndOrder(x,i+j*jStride+k*kStride,-1-jStride-kStride);} // corner
    }
  }
  level->cycles.boundary_conditions += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
