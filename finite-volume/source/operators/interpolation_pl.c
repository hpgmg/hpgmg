//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
static inline void InterpolateBlock_PL(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, blockCopy_type *block, int threads_per_block){
  // interpolate 3D array from read_i,j,k of read[] to write_i,j,k in write[]
  int write_dim_i   = block->dim.i<<1; // calculate the dimensions of the resultant fine block
  int write_dim_j   = block->dim.j<<1;
  int write_dim_k   = block->dim.k<<1;

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
  double * __restrict__ valid;
  if(block->read.box >=0){
     read = level_c->my_boxes[ block->read.box].components[         id_c] + level_c->my_boxes[ block->read.box].ghosts*(1+level_c->my_boxes[ block->read.box].jStride+level_c->my_boxes[ block->read.box].kStride);
     valid= level_c->my_boxes[ block->read.box].components[STENCIL_VALID] + level_c->my_boxes[ block->read.box].ghosts*(1+level_c->my_boxes[ block->read.box].jStride+level_c->my_boxes[ block->read.box].kStride);
     read_jStride = level_c->my_boxes[block->read.box ].jStride;
     read_kStride = level_c->my_boxes[block->read.box ].kStride;
  }
  if(block->write.box>=0){
    write = level_f->my_boxes[block->write.box].components[id_f] + level_f->my_boxes[block->write.box].ghosts*(1+level_f->my_boxes[block->write.box].jStride+level_f->my_boxes[block->write.box].kStride);
    write_jStride = level_f->my_boxes[block->write.box].jStride;
    write_kStride = level_f->my_boxes[block->write.box].kStride;
  }
 
 
  int i,j,k;
  double cim1,cip1;
  double cjm1,cjp1;
  double ckm1,ckp1;
  #pragma omp parallel for num_threads(threads_per_block) OMP_COLLAPSE
  for(k=0;k<write_dim_k;k++){
  for(j=0;j<write_dim_j;j++){
  for(i=0;i<write_dim_i;i++){
    int write_ijk = ((i   )+write_i) + (((j   )+write_j)*write_jStride) + (((k   )+write_k)*write_kStride);
    int  read_ijk = ((i>>1)+ read_i) + (((j>>1)+ read_j)* read_jStride) + (((k>>1)+ read_k)* read_kStride);
    //
    // |   o   |   o   |
    // +---+---+---+---+
    // |   | x | x |   |
    //
    // CAREFUL !!!  you must guarantee you zero'd the MPI buffers(write[]) and destination boxes at some point to avoid 0.0*NaN or 0.0*inf
    #if 1
    // piecewise linear interpolation... NOTE, BC's must have been previously applied
    int delta_i=           -1;if(i&0x1)delta_i=           1; // i.e. even points look backwards while odd points look forward
    int delta_j=-read_jStride;if(j&0x1)delta_j=read_jStride;
    int delta_k=-read_kStride;if(k&0x1)delta_k=read_kStride;
    write[write_ijk] = prescale_f*write[write_ijk] + 
        0.421875*read[read_ijk                        ] +
        0.140625*read[read_ijk                +delta_k] +
        0.140625*read[read_ijk        +delta_j        ] +
        0.046875*read[read_ijk        +delta_j+delta_k] +
        0.140625*read[read_ijk+delta_i                ] +
        0.046875*read[read_ijk+delta_i        +delta_k] +
        0.046875*read[read_ijk+delta_i+delta_j        ] +
        0.015625*read[read_ijk+delta_i+delta_j+delta_k];
    #endif
    #if 0
    // piecewise linear interpolation... NOTE, linear BC's have been fused into the stencil ala Mohr paper
    int delta_i=           -1;if(i&0x1)delta_i=           1; // i.e. even points look backwards while odd points look forward
    int delta_j=-read_jStride;if(j&0x1)delta_j=read_jStride;
    int delta_k=-read_kStride;if(k&0x1)delta_k=read_kStride;
    write[write_ijk] = prescale_f*write[write_ijk] + 
      0.015625*(     valid[read_ijk+delta_i] )*(     valid[read_ijk+delta_j] )*(     valid[read_ijk+delta_k] )*read[read_ijk+delta_i+delta_j+delta_k]+
      0.015625*(     valid[read_ijk+delta_i] )*(     valid[read_ijk+delta_j] )*( 2.0+valid[read_ijk+delta_k] )*read[read_ijk+delta_i+delta_j        ]+
      0.015625*(     valid[read_ijk+delta_i] )*( 2.0+valid[read_ijk+delta_j] )*(     valid[read_ijk+delta_k] )*read[read_ijk+delta_i        +delta_k]+
      0.015625*( 2.0+valid[read_ijk+delta_i] )*(     valid[read_ijk+delta_j] )*(     valid[read_ijk+delta_k] )*read[read_ijk        +delta_j+delta_k]+
      0.015625*(     valid[read_ijk+delta_i] )*( 2.0+valid[read_ijk+delta_j] )*( 2.0+valid[read_ijk+delta_k] )*read[read_ijk+delta_i                ]+
      0.015625*( 2.0+valid[read_ijk+delta_i] )*(     valid[read_ijk+delta_j] )*( 2.0+valid[read_ijk+delta_k] )*read[read_ijk        +delta_j        ]+
      0.015625*( 2.0+valid[read_ijk+delta_i] )*( 2.0+valid[read_ijk+delta_j] )*(     valid[read_ijk+delta_k] )*read[read_ijk                +delta_k]+
      0.015625*( 2.0+valid[read_ijk+delta_i] )*( 2.0+valid[read_ijk+delta_j] )*( 2.0+valid[read_ijk+delta_k] )*read[read_ijk                        ];
    #endif
    #if 0
    // gradient across cell... NOTE, BC's must have been previously applied
    double coefi = -0.125;if(i&0x1)coefi = 0.125;
    double coefj = -0.125;if(j&0x1)coefj = 0.125;
    double coefk = -0.125;if(k&0x1)coefk = 0.125;
    write[write_ijk] = prescale_f*write[write_ijk] +
                               read[read_ijk             ] + 
                        coefi*(read[read_ijk+           1]-read[read_ijk-           1]) +
                        coefj*(read[read_ijk+read_jStride]-read[read_ijk-read_jStride]) +
                        coefk*(read[read_ijk+read_kStride]-read[read_ijk-read_kStride]);
    #endif
    #if 0
    // Kwak interpolation... NOTE, BC's must have been previously applied
    int delta_i=           -1;if(i&0x1)delta_i=           1; // i.e. even points look backwards while odd points look forward
    int delta_j=-read_jStride;if(j&0x1)delta_j=read_jStride;
    int delta_k=-read_kStride;if(k&0x1)delta_k=read_kStride;
    write[write_ijk] = prescale_f*write[write_ijk] + 
        0.250*( read[read_ijk        ] ) +
        0.250*( read[read_ijk+delta_i] ) +
        0.250*( read[read_ijk+delta_j] ) +
        0.250*( read[read_ijk+delta_k] ) ;
    #endif
  }}}

}


//------------------------------------------------------------------------------------------------------------------------------
// perform a (inter-level) piecewise linear interpolation
void interpolation_pl(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){
  exchange_boundary(level_c,id_c,0);
  apply_BCs_linear(level_c,id_c);
  //apply_BCs_2ndOrder(level_c,id_c);
  //apply_BCs_4thOrder(level_c,id_c);

  uint64_t _timeCommunicationStart = CycleTime();
  uint64_t _timeStart,_timeEnd;
  int buffer=0;
  int sendBox,recvBox,n;

  #ifdef USE_MPI

  // loop through packed list of MPI receives and prepost Irecv's...
  _timeStart = CycleTime();
  #ifdef USE_MPI_THREAD_MULTIPLE
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for(n=0;n<level_f->interpolation.num_recvs;n++){
    MPI_Irecv(level_f->interpolation.recv_buffers[n],
              level_f->interpolation.recv_sizes[n],
              MPI_DOUBLE,
              level_f->interpolation.recv_ranks[n],
              0, // only one message should be received from each neighboring process
              MPI_COMM_WORLD,
              &level_f->interpolation.requests[n]
    );
  }
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_recv += (_timeEnd-_timeStart);


  // pack MPI send buffers...
  _timeStart = CycleTime();
  #pragma omp parallel for num_threads(level_f->concurrent_boxes) schedule(static,1)
  for(buffer=0;buffer<level_c->interpolation.num_blocks[0];buffer++){InterpolateBlock_PL(level_f,id_f,0.0,level_c,id_c,&level_c->interpolation.blocks[0][buffer],level_f->threads_per_box);} // !!! prescale==0 because you don't want to increment the MPI buffer
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_pack += (_timeEnd-_timeStart);

 
  // loop through MPI send buffers and post Isend's...
  _timeStart = CycleTime();
  #ifdef USE_MPI_THREAD_MULTIPLE
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for(n=0;n<level_c->interpolation.num_sends;n++){
    MPI_Isend(level_c->interpolation.send_buffers[n],
              level_c->interpolation.send_sizes[n],
              MPI_DOUBLE,
              level_c->interpolation.send_ranks[n],
              0, // only one message should be sent to each neighboring process
              MPI_COMM_WORLD,
              &level_c->interpolation.requests[n]
    );
  }
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_send += (_timeEnd-_timeStart);
  #endif


  // perform local interpolation... try and hide within Isend latency... 
  _timeStart = CycleTime();
  #pragma omp parallel for num_threads(level_f->concurrent_boxes) schedule(static,1)
  for(buffer=0;buffer<level_c->interpolation.num_blocks[1];buffer++){InterpolateBlock_PL(level_f,id_f,prescale_f,level_c,id_c,&level_c->interpolation.blocks[1][buffer],level_f->threads_per_box);}
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_local += (_timeEnd-_timeStart);


  // wait for MPI to finish...
  #ifdef USE_MPI 
  _timeStart = CycleTime();
  if(level_c->interpolation.num_sends)MPI_Waitall(level_c->interpolation.num_sends,level_c->interpolation.requests,level_c->interpolation.status);
  if(level_f->interpolation.num_recvs)MPI_Waitall(level_f->interpolation.num_recvs,level_f->interpolation.requests,level_f->interpolation.status);
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_wait += (_timeEnd-_timeStart);

  // unpack MPI receive buffers 
  _timeStart = CycleTime();
  #pragma omp parallel for num_threads(level_f->concurrent_boxes) schedule(static,1)
  for(buffer=0;buffer<level_f->interpolation.num_blocks[2];buffer++){IncrementBlock(level_f,id_f,prescale_f,&level_f->interpolation.blocks[2][buffer],level_f->threads_per_box);}
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_unpack += (_timeEnd-_timeStart);
  #endif 
 
 
  level_f->cycles.interpolation_total += (uint64_t)(CycleTime()-_timeCommunicationStart);
}
