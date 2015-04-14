//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
static inline void interpolation_v4_block(level_type *level_f, int id_f, double prescale_f, level_type *level_c, int id_c, blockCopy_type *block){
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
  if(block->read.box >=0){
     read_jStride = level_c->my_boxes[block->read.box ].jStride;
     read_kStride = level_c->my_boxes[block->read.box ].kStride;
     read = level_c->my_boxes[ block->read.box].vectors[id_c] + level_c->my_boxes[ block->read.box].ghosts*(1+ read_jStride+ read_kStride);
  }
  if(block->write.box>=0){
    write_jStride = level_f->my_boxes[block->write.box].jStride;
    write_kStride = level_f->my_boxes[block->write.box].kStride;
    write = level_f->my_boxes[block->write.box].vectors[id_f] + level_f->my_boxes[block->write.box].ghosts*(1+write_jStride+write_kStride);
  }
 

  int i,j,k;
  double c2 = -3.0/128.0;
  double c1 = 22.0/128.0;
  int dj  =   read_jStride;
  int dk  =   read_kStride;
  int dj2 = 2*read_jStride;
  int dk2 = 2*read_kStride;
  for(k=0;k<write_dim_k;k++){double sk1=c1,sk2=c2;if(k&0x1){sk1=-c1;sk2=-c2;}
  for(j=0;j<write_dim_j;j++){double sj1=c1,sj2=c2;if(j&0x1){sj1=-c1;sj2=-c2;}
  for(i=0;i<write_dim_i;i++){double si1=c1,si2=c2;if(i&0x1){si1=-c1;si2=-c2;}
    int write_ijk = ((i   )+write_i) + (((j   )+write_j)*write_jStride) + (((k   )+write_k)*write_kStride);
    int  read_ijk = ((i>>1)+ read_i) + (((j>>1)+ read_j)* read_jStride) + (((k>>1)+ read_k)* read_kStride);
    //
    // |   -3/128  |  +22/128  |    1.0    |  -22/128  |   +3/128  | coarse grid
    // |-----+-----|-----+-----|-----+-----|-----+-----|-----+-----|
    // |     |     |     |     |?????|     |     |     |     |     | fine grid
    //
    write[write_ijk] = prescale_f*write[write_ijk] +
                       + sk2*( + sj2*( si2*read[read_ijk-2-dj2-dk2] + si1*read[read_ijk-1-dj2-dk2] + read[read_ijk-dj2-dk2] - si1*read[read_ijk+1-dj2-dk2] - si2*read[read_ijk+2-dj2-dk2] )
                               + sj1*( si2*read[read_ijk-2-dj -dk2] + si1*read[read_ijk-1-dj -dk2] + read[read_ijk-dj -dk2] - si1*read[read_ijk+1-dj -dk2] - si2*read[read_ijk+2-dj -dk2] )
                               +     ( si2*read[read_ijk-2    -dk2] + si1*read[read_ijk-1    -dk2] + read[read_ijk    -dk2] - si1*read[read_ijk+1    -dk2] - si2*read[read_ijk+2    -dk2] )
                               - sj1*( si2*read[read_ijk-2+dj -dk2] + si1*read[read_ijk-1+dj -dk2] + read[read_ijk+dj -dk2] - si1*read[read_ijk+1+dj -dk2] - si2*read[read_ijk+2+dj -dk2] )
                               - sj2*( si2*read[read_ijk-2+dj2-dk2] + si1*read[read_ijk-1+dj2-dk2] + read[read_ijk+dj2-dk2] - si1*read[read_ijk+1+dj2-dk2] - si2*read[read_ijk+2+dj2-dk2] ) )
                       + sk1*( + sj2*( si2*read[read_ijk-2-dj2-dk ] + si1*read[read_ijk-1-dj2-dk ] + read[read_ijk-dj2-dk ] - si1*read[read_ijk+1-dj2-dk ] - si2*read[read_ijk+2-dj2-dk ] )
                               + sj1*( si2*read[read_ijk-2-dj -dk ] + si1*read[read_ijk-1-dj -dk ] + read[read_ijk-dj -dk ] - si1*read[read_ijk+1-dj -dk ] - si2*read[read_ijk+2-dj -dk ] )
                               +     ( si2*read[read_ijk-2    -dk ] + si1*read[read_ijk-1    -dk ] + read[read_ijk    -dk ] - si1*read[read_ijk+1    -dk ] - si2*read[read_ijk+2    -dk ] )
                               - sj1*( si2*read[read_ijk-2+dj -dk ] + si1*read[read_ijk-1+dj -dk ] + read[read_ijk+dj -dk ] - si1*read[read_ijk+1+dj -dk ] - si2*read[read_ijk+2+dj -dk ] )
                               - sj2*( si2*read[read_ijk-2+dj2-dk ] + si1*read[read_ijk-1+dj2-dk ] + read[read_ijk+dj2-dk ] - si1*read[read_ijk+1+dj2-dk ] - si2*read[read_ijk+2+dj2-dk ] ) )
                       +     ( + sj2*( si2*read[read_ijk-2-dj2    ] + si1*read[read_ijk-1-dj2    ] + read[read_ijk-dj2    ] - si1*read[read_ijk+1-dj2    ] - si2*read[read_ijk+2-dj2    ] )
                               + sj1*( si2*read[read_ijk-2-dj     ] + si1*read[read_ijk-1-dj     ] + read[read_ijk-dj     ] - si1*read[read_ijk+1-dj     ] - si2*read[read_ijk+2-dj     ] )
                               +     ( si2*read[read_ijk-2        ] + si1*read[read_ijk-1        ] + read[read_ijk        ] - si1*read[read_ijk+1        ] - si2*read[read_ijk+2        ] )
                               - sj1*( si2*read[read_ijk-2+dj     ] + si1*read[read_ijk-1+dj     ] + read[read_ijk+dj     ] - si1*read[read_ijk+1+dj     ] - si2*read[read_ijk+2+dj     ] )
                               - sj2*( si2*read[read_ijk-2+dj2    ] + si1*read[read_ijk-1+dj2    ] + read[read_ijk+dj2    ] - si1*read[read_ijk+1+dj2    ] - si2*read[read_ijk+2+dj2    ] ) )
                       - sk1*( + sj2*( si2*read[read_ijk-2-dj2+dk ] + si1*read[read_ijk-1-dj2+dk ] + read[read_ijk-dj2+dk ] - si1*read[read_ijk+1-dj2+dk ] - si2*read[read_ijk+2-dj2+dk ] )
                               + sj1*( si2*read[read_ijk-2-dj +dk ] + si1*read[read_ijk-1-dj +dk ] + read[read_ijk-dj +dk ] - si1*read[read_ijk+1-dj +dk ] - si2*read[read_ijk+2-dj +dk ] )
                               +     ( si2*read[read_ijk-2    +dk ] + si1*read[read_ijk-1    +dk ] + read[read_ijk    +dk ] - si1*read[read_ijk+1    +dk ] - si2*read[read_ijk+2    +dk ] )
                               - sj1*( si2*read[read_ijk-2+dj +dk ] + si1*read[read_ijk-1+dj +dk ] + read[read_ijk+dj +dk ] - si1*read[read_ijk+1+dj +dk ] - si2*read[read_ijk+2+dj +dk ] )
                               - sj2*( si2*read[read_ijk-2+dj2+dk ] + si1*read[read_ijk-1+dj2+dk ] + read[read_ijk+dj2+dk ] - si1*read[read_ijk+1+dj2+dk ] - si2*read[read_ijk+2+dj2+dk ] ) )
                       - sk2*( + sj2*( si2*read[read_ijk-2-dj2+dk2] + si1*read[read_ijk-1-dj2+dk2] + read[read_ijk-dj2+dk2] - si1*read[read_ijk+1-dj2+dk2] - si2*read[read_ijk+2-dj2+dk2] )
                               + sj1*( si2*read[read_ijk-2-dj +dk2] + si1*read[read_ijk-1-dj +dk2] + read[read_ijk-dj +dk2] - si1*read[read_ijk+1-dj +dk2] - si2*read[read_ijk+2-dj +dk2] )
                               +     ( si2*read[read_ijk-2    +dk2] + si1*read[read_ijk-1    +dk2] + read[read_ijk    +dk2] - si1*read[read_ijk+1    +dk2] - si2*read[read_ijk+2    +dk2] )
                               - sj1*( si2*read[read_ijk-2+dj +dk2] + si1*read[read_ijk-1+dj +dk2] + read[read_ijk+dj +dk2] - si1*read[read_ijk+1+dj +dk2] - si2*read[read_ijk+2+dj +dk2] )
                               - sj2*( si2*read[read_ijk-2+dj2+dk2] + si1*read[read_ijk-1+dj2+dk2] + read[read_ijk+dj2+dk2] - si1*read[read_ijk+1+dj2+dk2] - si2*read[read_ijk+2+dj2+dk2] ) );
  }}}

}


//------------------------------------------------------------------------------------------------------------------------------
// perform a (inter-level) volumetric quartic interpolation
void interpolation_v4(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){
    exchange_boundary(level_c,id_c,0);
         apply_BCs_v4(level_c,id_c,0);

  uint64_t _timeCommunicationStart = CycleTime();
  uint64_t _timeStart,_timeEnd;
  int buffer=0;
  int n;
  int my_tag = (level_f->tag<<4) | 0x7;


  #ifdef USE_MPI
  // by convention, level_f allocates a combined array of requests for both level_f recvs and level_c sends...
  int nMessages = level_c->interpolation.num_sends + level_f->interpolation.num_recvs;
  MPI_Request *recv_requests = level_f->interpolation.requests;
  MPI_Request *send_requests = level_f->interpolation.requests + level_f->interpolation.num_recvs;


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
              my_tag,
              MPI_COMM_WORLD,
              &recv_requests[n]
    );
  }
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_recv += (_timeEnd-_timeStart);


  // pack MPI send buffers...
  _timeStart = CycleTime();
  PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_c->interpolation.num_blocks[0])
  for(buffer=0;buffer<level_c->interpolation.num_blocks[0];buffer++){
    // !!! prescale==0 because you don't want to increment the MPI buffer
    interpolation_v4_block(level_f,id_f,0.0,level_c,id_c,&level_c->interpolation.blocks[0][buffer]);
  }
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
              my_tag,
              MPI_COMM_WORLD,
              &send_requests[n]
    );
  }
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_send += (_timeEnd-_timeStart);
  #endif


  // perform local interpolation... try and hide within Isend latency... 
  _timeStart = CycleTime();
  PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_c->interpolation.num_blocks[1])
  for(buffer=0;buffer<level_c->interpolation.num_blocks[1];buffer++){
    interpolation_v4_block(level_f,id_f,prescale_f,level_c,id_c,&level_c->interpolation.blocks[1][buffer]);
  }
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_local += (_timeEnd-_timeStart);


  // wait for MPI to finish...
  #ifdef USE_MPI 
  _timeStart = CycleTime();
  if(nMessages)MPI_Waitall(nMessages,level_f->interpolation.requests,level_f->interpolation.status);
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_wait += (_timeEnd-_timeStart);


  // unpack MPI receive buffers 
  _timeStart = CycleTime();
  PRAGMA_THREAD_ACROSS_BLOCKS(level_f,buffer,level_f->interpolation.num_blocks[2])
  for(buffer=0;buffer<level_f->interpolation.num_blocks[2];buffer++){
    IncrementBlock(level_f,id_f,prescale_f,&level_f->interpolation.blocks[2][buffer]);
  }
  _timeEnd = CycleTime();
  level_f->cycles.interpolation_unpack += (_timeEnd-_timeStart);


  #endif 
 
 
  level_f->cycles.interpolation_total += (uint64_t)(CycleTime()-_timeCommunicationStart);
}
