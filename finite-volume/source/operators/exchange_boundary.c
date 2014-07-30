//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// perform a (intra-level) ghost zone exchange
//  NOTE exchange_boundary() only exchanges the boundary.  
//  It will not enforce any boundary conditions
//  BC's are either the responsibility of a separate function or should be fused into the stencil
void exchange_boundary(level_type * level, int id, int justFaces){
  uint64_t _timeCommunicationStart = CycleTime();
  uint64_t _timeStart,_timeEnd;
  int buffer=0;
  int n;

  if(justFaces)justFaces=1;else justFaces=0;  // must be 0 or 1 in order to index into exchange_ghosts[]

  #ifdef USE_MPI
  int nMessages = level->exchange_ghosts[justFaces].num_recvs + level->exchange_ghosts[justFaces].num_sends;
  MPI_Request *recv_requests = level->exchange_ghosts[justFaces].requests;
  MPI_Request *send_requests = level->exchange_ghosts[justFaces].requests + level->exchange_ghosts[justFaces].num_recvs;

  // loop through packed list of MPI receives and prepost Irecv's...
  _timeStart = CycleTime();
  #ifdef USE_MPI_THREAD_MULTIPLE
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for(n=0;n<level->exchange_ghosts[justFaces].num_recvs;n++){
    MPI_Irecv(level->exchange_ghosts[justFaces].recv_buffers[n],
              level->exchange_ghosts[justFaces].recv_sizes[n],
              MPI_DOUBLE,
              level->exchange_ghosts[justFaces].recv_ranks[n],
              0, // by convention, ghost zone exchanges use tag=0
              MPI_COMM_WORLD,
              //&level->exchange_ghosts[justFaces].requests[n]
              &recv_requests[n]
    );
  }
  _timeEnd = CycleTime();
  level->cycles.ghostZone_recv += (_timeEnd-_timeStart);


  // pack MPI send buffers...
  _timeStart = CycleTime();
  #pragma omp parallel for if(level->exchange_ghosts[justFaces].num_blocks[0]>1) schedule(static,1)
  for(buffer=0;buffer<level->exchange_ghosts[justFaces].num_blocks[0];buffer++){CopyBlock(level,id,&level->exchange_ghosts[justFaces].blocks[0][buffer]);}
  _timeEnd = CycleTime();
  level->cycles.ghostZone_pack += (_timeEnd-_timeStart);

 
  // loop through MPI send buffers and post Isend's...
  _timeStart = CycleTime();
  #ifdef USE_MPI_THREAD_MULTIPLE
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for(n=0;n<level->exchange_ghosts[justFaces].num_sends;n++){
    MPI_Isend(level->exchange_ghosts[justFaces].send_buffers[n],
              level->exchange_ghosts[justFaces].send_sizes[n],
              MPI_DOUBLE,
              level->exchange_ghosts[justFaces].send_ranks[n],
              0, // by convention, ghost zone exchanges use tag=0
              MPI_COMM_WORLD,
              &send_requests[n]
              //&level->exchange_ghosts[justFaces].requests[n+level->exchange_ghosts[justFaces].num_recvs]
                                              // requests[0..num_recvs-1] were used by recvs.  So sends start at num_recvs
    ); 
  }
  _timeEnd = CycleTime();
  level->cycles.ghostZone_send += (_timeEnd-_timeStart);
  #endif


  // exchange locally... try and hide within Isend latency... 
  _timeStart = CycleTime();
  #pragma omp parallel for if(level->exchange_ghosts[justFaces].num_blocks[1]>1) schedule(static,1)
  for(buffer=0;buffer<level->exchange_ghosts[justFaces].num_blocks[1];buffer++){CopyBlock(level,id,&level->exchange_ghosts[justFaces].blocks[1][buffer]);}
  _timeEnd = CycleTime();
  level->cycles.ghostZone_local += (_timeEnd-_timeStart);


  // wait for MPI to finish...
  #ifdef USE_MPI 
  _timeStart = CycleTime();
  if(nMessages)MPI_Waitall(nMessages,level->exchange_ghosts[justFaces].requests,level->exchange_ghosts[justFaces].status);
  _timeEnd = CycleTime();
  level->cycles.ghostZone_wait += (_timeEnd-_timeStart);


  // unpack MPI receive buffers 
  _timeStart = CycleTime();
  #pragma omp parallel for if(level->exchange_ghosts[justFaces].num_blocks[2]>1) schedule(static,1)
  for(buffer=0;buffer<level->exchange_ghosts[justFaces].num_blocks[2];buffer++){CopyBlock(level,id,&level->exchange_ghosts[justFaces].blocks[2][buffer]);}
  _timeEnd = CycleTime();
  level->cycles.ghostZone_unpack += (_timeEnd-_timeStart);
  #endif

 
  level->cycles.ghostZone_total += (uint64_t)(CycleTime()-_timeCommunicationStart);
}
