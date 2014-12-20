//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void zero_vector(level_type * level, int component_id){
  // zero's the entire grid INCLUDING ghost zones...
  uint64_t _timeStart = CycleTime();
  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid = level->my_boxes[box].vectors[component_id] + ghosts*(1+jStride+kStride);
    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=-ghosts;k<dim+ghosts;k++){
    for(j=-ghosts;j<dim+ghosts;j++){
    for(i=-ghosts;i<dim+ghosts;i++){
      int ijk = i + j*jStride + k*kStride;
      grid[ijk] = 0.0;
    }}}
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void initialize_valid_region(level_type * level){
  uint64_t _timeStart = CycleTime();
  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ valid = level->my_boxes[box].vectors[VECTOR_VALID] + ghosts*(1+jStride+kStride);
    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=-ghosts;k<dim+ghosts;k++){
    for(j=-ghosts;j<dim+ghosts;j++){
    for(i=-ghosts;i<dim+ghosts;i++){
      int ijk = i + j*jStride + k*kStride;
      valid[ijk] = 1.0; // i.e. all cells including ghosts are valid for periodic BC's
      if(level->boundary_condition.type == BC_DIRICHLET){ // cells outside the domain boundaries are not valid
        if(i + level->my_boxes[box].low.i <             0)valid[ijk] = 0.0;
        if(j + level->my_boxes[box].low.j <             0)valid[ijk] = 0.0;
        if(k + level->my_boxes[box].low.k <             0)valid[ijk] = 0.0;
        if(i + level->my_boxes[box].low.i >= level->dim.i)valid[ijk] = 0.0;
        if(j + level->my_boxes[box].low.j >= level->dim.j)valid[ijk] = 0.0;
        if(k + level->my_boxes[box].low.k >= level->dim.k)valid[ijk] = 0.0;
      }
    }}}
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void initialize_grid_to_scalar(level_type * level, int component_id, double scalar){
  // initializes the grid to a scalar while zero'ing the ghost zones...
  uint64_t _timeStart = CycleTime();
  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid = level->my_boxes[box].vectors[component_id] + ghosts*(1+jStride+kStride);
    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=-ghosts;k<dim+ghosts;k++){
    for(j=-ghosts;j<dim+ghosts;j++){
    for(i=-ghosts;i<dim+ghosts;i++){
        int ijk = i + j*jStride + k*kStride;
        int ghostZone = (i<0) || (j<0) || (k<0) || (i>=dim) || (j>=dim) || (k>=dim);
        grid[ijk] = ghostZone ? 0.0 : scalar;
    }}}
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void add_vectors(level_type * level, int id_c, double scale_a, int id_a, double scale_b, int id_b){ // c=scale_a*id_a + scale_b*id_b
  uint64_t _timeStart = CycleTime();

  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_b = level->my_boxes[box].vectors[id_b] + ghosts*(1+jStride+kStride);
    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
        int ijk = i + j*jStride + k*kStride;
        grid_c[ijk] = scale_a*grid_a[ijk] + scale_b*grid_b[ijk];
    }}}
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void mul_vectors(level_type * level, int id_c, double scale, int id_a, int id_b){ // id_c=scale*id_a*id_b
  uint64_t _timeStart = CycleTime();

  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_b = level->my_boxes[box].vectors[id_b] + ghosts*(1+jStride+kStride);
    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
        int ijk = i + j*jStride + k*kStride;
        grid_c[ijk] = scale*grid_a[ijk]*grid_b[ijk];
    }}}
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void invert_vector(level_type * level, int id_c, double scale_a, int id_a){ // c[]=scale_a/a[]
  uint64_t _timeStart = CycleTime();

  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);
    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
        int ijk = i + j*jStride + k*kStride;
        grid_c[ijk] = scale_a/grid_a[ijk];
    }}}
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
void scale_vector(level_type * level, int id_c, double scale_a, int id_a){ // c[]=scale_a*a[]
  uint64_t _timeStart = CycleTime();

  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride);
    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
        int ijk = i + j*jStride + k*kStride;
        grid_c[ijk] = scale_a*grid_a[ijk];
    }}}
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
double dot(level_type * level, int id_a, int id_b){
  uint64_t _timeStart = CycleTime();


  int box;
  double a_dot_b_level =  0.0;
  // FIX, schedule(static) is a stand in to guarantee reproducibility...
  PRAGMA_THREAD_ACROSS_BOXES_SUM(level,box,a_dot_b_level)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double * __restrict__ grid_b = level->my_boxes[box].vectors[id_b] + ghosts*(1+jStride+kStride);
    double a_dot_b_box = 0.0;
    PRAGMA_THREAD_WITHIN_A_BOX_SUM(level,i,j,k,a_dot_b_box)
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
      int ijk = i + j*jStride + k*kStride;
      a_dot_b_box += grid_a[ijk]*grid_b[ijk];
    }}}
    a_dot_b_level+=a_dot_b_box;
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);

  #ifdef USE_MPI
  uint64_t _timeStartAllReduce = CycleTime();
  double send = a_dot_b_level;
  MPI_Allreduce(&send,&a_dot_b_level,1,MPI_DOUBLE,MPI_SUM,level->MPI_COMM_ALLREDUCE);
  uint64_t _timeEndAllReduce = CycleTime();
  level->cycles.collectives   += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce);
  #endif

  return(a_dot_b_level);
}

//------------------------------------------------------------------------------------------------------------------------------
double norm(level_type * level, int component_id){ // implements the max norm
  uint64_t _timeStart = CycleTime();

  int box;
  double max_norm =  0.0;
  // FIX, schedule(static) is a stand in to guarantee reproducibility...
  PRAGMA_THREAD_ACROSS_BOXES_MAX(level,box,max_norm)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid   = level->my_boxes[box].vectors[component_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double box_norm = 0.0;
    PRAGMA_THREAD_WITHIN_A_BOX_MAX(level,i,j,k,box_norm)
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
      int ijk = i + j*jStride + k*kStride;
      double fabs_grid_ijk = fabs(grid[ijk]);
      if(fabs_grid_ijk>box_norm){box_norm=fabs_grid_ijk;} // max norm
    }}}
    if(box_norm>max_norm){max_norm = box_norm;}
  } // box list
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);

  #ifdef USE_MPI
  uint64_t _timeStartAllReduce = CycleTime();
  double send = max_norm;
  MPI_Allreduce(&send,&max_norm,1,MPI_DOUBLE,MPI_MAX,level->MPI_COMM_ALLREDUCE);
  uint64_t _timeEndAllReduce = CycleTime();
  level->cycles.collectives   += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce);
  #endif
  return(max_norm);
}


//------------------------------------------------------------------------------------------------------------------------------
double mean(level_type * level, int id_a){
  uint64_t _timeStart = CycleTime();


  int box;
  double sum_level =  0.0;
  PRAGMA_THREAD_ACROSS_BOXES_SUM(level,box,sum_level)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double sum_box = 0.0;
    PRAGMA_THREAD_WITHIN_A_BOX_SUM(level,i,j,k,sum_box)
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
      int ijk = i + j*jStride + k*kStride;
      sum_box += grid_a[ijk];
    }}}
    sum_level+=sum_box;
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
  double ncells_level = (double)level->dim.i*(double)level->dim.j*(double)level->dim.k;

  #ifdef USE_MPI
  uint64_t _timeStartAllReduce = CycleTime();
  double send = sum_level;
  MPI_Allreduce(&send,&sum_level,1,MPI_DOUBLE,MPI_SUM,level->MPI_COMM_ALLREDUCE);
  uint64_t _timeEndAllReduce = CycleTime();
  level->cycles.collectives   += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce);
  #endif

  double mean_level = sum_level / ncells_level;
  return(mean_level);
}


void shift_vector(level_type * level, int id_c, int id_a, double shift_a){
  uint64_t _timeStart = CycleTime();


  int box;
  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid_c = level->my_boxes[box].vectors[id_c] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
    double * __restrict__ grid_a = level->my_boxes[box].vectors[id_a] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point

    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
      int ijk = i + j*jStride + k*kStride;
      grid_c[ijk] = grid_a[ijk] + shift_a;
    }}}
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
}

//------------------------------------------------------------------------------------------------------------------------------
void project_cell_to_face(level_type * level, int id_cell, int id_face, int dir){
  uint64_t _timeStart = CycleTime();
  int box;

  PRAGMA_THREAD_ACROSS_BOXES(level,box)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    const int     dim = level->my_boxes[box].dim;
    double * __restrict__ grid_cell = level->my_boxes[box].vectors[id_cell] + ghosts*(1+jStride+kStride);
    double * __restrict__ grid_face = level->my_boxes[box].vectors[id_face] + ghosts*(1+jStride+kStride);
    int stride;
    switch(dir){
      case 0: stride =       1;break;//i-direction
      case 1: stride = jStride;break;//j-direction
      case 2: stride = kStride;break;//k-direction
    }
    PRAGMA_THREAD_WITHIN_A_BOX(level,i,j,k)
    for(k=0;k<=dim;k++){ // <= to ensure you do low and high faces
    for(j=0;j<=dim;j++){
    for(i=0;i<=dim;i++){
      int ijk = i + j*jStride + k*kStride;
      grid_face[ijk] = 0.5*(grid_cell[ijk-stride] + grid_cell[ijk]); // simple linear interpolation
    }}}
  }

  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
double error(level_type * level, int id_a, int id_b){
  double h3 = level->h * level->h * level->h;
               add_vectors(level,VECTOR_TEMP,1.0,id_a,-1.0,id_b);            // VECTOR_TEMP = id_a - id_b
  double   max =      norm(level,VECTOR_TEMP);                return(max);   // max norm of error function
  double    L2 = sqrt( dot(level,VECTOR_TEMP,VECTOR_TEMP)*h3);return( L2);   // normalized L2 error ?
}
