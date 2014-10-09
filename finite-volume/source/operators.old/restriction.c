//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
static inline void RestrictBlock(level_type *level_c, int id_c, level_type *level_f, int id_f, blockCopy_type *block, int restrictionType){
  // restrict 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int   dim_i       = block->dim.i; // calculate the dimensions of the resultant coarse block
  int   dim_j       = block->dim.j;
  int   dim_k       = block->dim.k;

  int  read_i       = block->read.i;
  int  read_j       = block->read.j;
  int  read_k       = block->read.k;
  int  read_jStride = block->read.jStride;
  int  read_kStride = block->read.kStride;

  int write_i       = block->write.i;
  int write_j       = block->write.j;
  int write_k       = block->write.k;
  int write_jStride = block->write.jStride;
  int write_kStride = block->write.kStride;

  double * __restrict__  read = block->read.ptr;
  double * __restrict__ write = block->write.ptr;
  if(block->read.box >=0){
     read = level_f->my_boxes[ block->read.box].vectors[id_f] + level_f->my_boxes[ block->read.box].ghosts*(1+level_f->my_boxes[ block->read.box].jStride+level_f->my_boxes[ block->read.box].kStride);
     read_jStride = level_f->my_boxes[block->read.box ].jStride;
     read_kStride = level_f->my_boxes[block->read.box ].kStride;
  }
  if(block->write.box>=0){
    write = level_c->my_boxes[block->write.box].vectors[id_c] + level_c->my_boxes[block->write.box].ghosts*(1+level_c->my_boxes[block->write.box].jStride+level_c->my_boxes[block->write.box].kStride);
    write_jStride = level_c->my_boxes[block->write.box].jStride;
    write_kStride = level_c->my_boxes[block->write.box].kStride;
  }



  int i,j,k;
  switch(restrictionType){
    case RESTRICT_CELL:
         for(k=0;k<dim_k;k++){
         for(j=0;j<dim_j;j++){
         for(i=0;i<dim_i;i++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
           write[write_ijk] = ( read[read_ijk                            ]+read[read_ijk+1                          ] +
                                read[read_ijk  +read_jStride             ]+read[read_ijk+1+read_jStride             ] +
                                read[read_ijk               +read_kStride]+read[read_ijk+1             +read_kStride] +
                                read[read_ijk  +read_jStride+read_kStride]+read[read_ijk+1+read_jStride+read_kStride] ) * 0.125;
         }}}break;
    case RESTRICT_FACE_I:
         for(k=0;k<dim_k;k++){
         for(j=0;j<dim_j;j++){
         for(i=0;i<dim_i;i++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
           write[write_ijk] = ( read[read_ijk                          ] +
                                read[read_ijk+read_jStride             ] +
                                read[read_ijk             +read_kStride] +
                                read[read_ijk+read_jStride+read_kStride] ) * 0.25;
         }}}break;
    case RESTRICT_FACE_J:
         for(k=0;k<dim_k;k++){
         for(j=0;j<dim_j;j++){
         for(i=0;i<dim_i;i++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
           write[write_ijk] = ( read[read_ijk               ] +
                                read[read_ijk+1             ] +
                                read[read_ijk  +read_kStride] +
                                read[read_ijk+1+read_kStride] ) * 0.25;
         }}}break;
    case RESTRICT_FACE_K:
         for(k=0;k<dim_k;k++){
         for(j=0;j<dim_j;j++){
         for(i=0;i<dim_i;i++){
           int write_ijk = ((i   )+write_i) + ((j   )+write_j)*write_jStride + ((k   )+write_k)*write_kStride;
           int  read_ijk = ((i<<1)+ read_i) + ((j<<1)+ read_j)* read_jStride + ((k<<1)+ read_k)* read_kStride;
           write[write_ijk] = ( read[read_ijk               ] +
                                read[read_ijk+1             ] +
                                read[read_ijk  +read_jStride] +
                                read[read_ijk+1+read_jStride] ) * 0.25;
         }}}break;
  }

}


//------------------------------------------------------------------------------------------------------------------------------
// perform a (inter-level) restriction
void restriction(level_type * level_c, int id_c, level_type *level_f, int id_f, int restrictionType){
  uint64_t _timeCommunicationStart = CycleTime();
  uint64_t _timeStart,_timeEnd;
  int buffer=0;
  int n;




  #ifdef USE_MPI
  // by convention, level_f allocates a combined array of requests for both level_f sends and level_c recvs...
  int nMessages = level_c->restriction[restrictionType].num_recvs + level_f->restriction[restrictionType].num_sends;
  MPI_Request *recv_requests = level_f->restriction[restrictionType].requests;
  MPI_Request *send_requests = level_f->restriction[restrictionType].requests + level_c->restriction[restrictionType].num_recvs;


  // loop through packed list of MPI receives and prepost Irecv's...
  _timeStart = CycleTime();
  #ifdef USE_MPI_THREAD_MULTIPLE
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for(n=0;n<level_c->restriction[restrictionType].num_recvs;n++){
    MPI_Irecv(level_c->restriction[restrictionType].recv_buffers[n],
              level_c->restriction[restrictionType].recv_sizes[n],
              MPI_DOUBLE,
              level_c->restriction[restrictionType].recv_ranks[n],
              5, // by convention, restriction uses tag=5
              MPI_COMM_WORLD,
              &recv_requests[n]
    );
  }
  _timeEnd = CycleTime();
  level_f->cycles.restriction_recv += (_timeEnd-_timeStart);


  // pack MPI send buffers...
  _timeStart = CycleTime();
  #pragma omp parallel for private(buffer) if(level_f->restriction[restrictionType].num_blocks[0]>1) schedule(static,1)
  for(buffer=0;buffer<level_f->restriction[restrictionType].num_blocks[0];buffer++){RestrictBlock(level_c,id_c,level_f,id_f,&level_f->restriction[restrictionType].blocks[0][buffer],restrictionType);}
  _timeEnd = CycleTime();
  level_f->cycles.restriction_pack += (_timeEnd-_timeStart);

 
  // loop through MPI send buffers and post Isend's...
  _timeStart = CycleTime();
  #ifdef USE_MPI_THREAD_MULTIPLE
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for(n=0;n<level_f->restriction[restrictionType].num_sends;n++){
    MPI_Isend(level_f->restriction[restrictionType].send_buffers[n],
              level_f->restriction[restrictionType].send_sizes[n],
              MPI_DOUBLE,
              level_f->restriction[restrictionType].send_ranks[n],
              5, // by convention, restriction uses tag=5
              MPI_COMM_WORLD,
              &send_requests[n]
    );
  }
  _timeEnd = CycleTime();
  level_f->cycles.restriction_send += (_timeEnd-_timeStart);
  #endif


  // perform local restriction[restrictionType]... try and hide within Isend latency... 
  _timeStart = CycleTime();
  #pragma omp parallel for private(buffer) if(level_f->restriction[restrictionType].num_blocks[1]>1) schedule(static,1)
  for(buffer=0;buffer<level_f->restriction[restrictionType].num_blocks[1];buffer++){RestrictBlock(level_c,id_c,level_f,id_f,&level_f->restriction[restrictionType].blocks[1][buffer],restrictionType);}
  _timeEnd = CycleTime();
  level_f->cycles.restriction_local += (_timeEnd-_timeStart);


  // wait for MPI to finish...
  #ifdef USE_MPI 
  _timeStart = CycleTime();
  if(nMessages)MPI_Waitall(nMessages,level_f->restriction[restrictionType].requests,level_f->restriction[restrictionType].status);
  _timeEnd = CycleTime();
  level_f->cycles.restriction_wait += (_timeEnd-_timeStart);


  // unpack MPI receive buffers 
  _timeStart = CycleTime();
  #pragma omp parallel for private(buffer) if(level_c->restriction[restrictionType].num_blocks[2]>1) schedule(static,1)
  for(buffer=0;buffer<level_c->restriction[restrictionType].num_blocks[2];buffer++){CopyBlock(level_c,id_c,&level_c->restriction[restrictionType].blocks[2][buffer]);}
  _timeEnd = CycleTime();
  level_f->cycles.restriction_unpack += (_timeEnd-_timeStart);


  #endif
 
 
  level_f->cycles.restriction_total += (uint64_t)(CycleTime()-_timeCommunicationStart);
}
