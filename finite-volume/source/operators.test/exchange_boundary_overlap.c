//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// perform a (intra-level) ghost zone exchange
//  NOTE exchange_boundary() only exchanges the boundary.  
//  It will not enforce any boundary conditions
//  BC's are either the responsibility of a separate function or should be fused into the stencil
void exchange_boundary(level_type * level, int id, int shape){
  double _timeCommunicationStart = getTime();

  if(shape>=STENCIL_MAX_SHAPES)shape=STENCIL_SHAPE_BOX;  // shape must be < STENCIL_MAX_SHAPES in order to safely index into exchange_ghosts[]
  #pragma omp parallel
  {
  double _timeStart;
  int tid = omp_get_thread_num();
  int nth = omp_get_num_threads();

  int my_tag = (level->tag<<4) | shape;
  int buffer=0;
  int n;

  #ifdef USE_MPI
  int nMessages = level->exchange_ghosts[shape].num_recvs + level->exchange_ghosts[shape].num_sends;
  MPI_Request *recv_requests = level->exchange_ghosts[shape].requests;
  MPI_Request *send_requests = level->exchange_ghosts[shape].requests + level->exchange_ghosts[shape].num_recvs;

  // loop through packed list of MPI receives and prepost Irecv's...
  if(tid==0){
  if(level->exchange_ghosts[shape].num_recvs>0){
    _timeStart = getTime();
    for(n=0;n<level->exchange_ghosts[shape].num_recvs;n++){
      MPI_Irecv(level->exchange_ghosts[shape].recv_buffers[n],
                level->exchange_ghosts[shape].recv_sizes[n],
                MPI_DOUBLE,
                level->exchange_ghosts[shape].recv_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &recv_requests[n]
      );
    }
    level->timers.ghostZone_recv += (getTime()-_timeStart);
  }
  }


  // pack MPI send buffers...
  if(tid>0){
  if(level->exchange_ghosts[shape].num_blocks[0]){
    if(tid==1)_timeStart = getTime();
    int lo = ((tid-1)*level->exchange_ghosts[shape].num_blocks[0])/(nth-1);
    int hi = ((tid  )*level->exchange_ghosts[shape].num_blocks[0])/(nth-1);
    for(buffer=lo;buffer<hi;buffer++){
      CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[0][buffer]);
    }
    if(tid==1)level->timers.ghostZone_pack += (getTime()-_timeStart);
  }}
  #pragma omp barrier // wait for packing


  // loop through MPI send buffers and post Isend's...
  if(tid==0){
  if(level->exchange_ghosts[shape].num_sends>0){
    _timeStart = getTime();
    for(n=0;n<level->exchange_ghosts[shape].num_sends;n++){
      MPI_Isend(level->exchange_ghosts[shape].send_buffers[n],
                level->exchange_ghosts[shape].send_sizes[n],
                MPI_DOUBLE,
                level->exchange_ghosts[shape].send_ranks[n],
                my_tag,
                MPI_COMM_WORLD,
                &send_requests[n]
      ); 
    }
    level->timers.ghostZone_send += (getTime()-_timeStart);
  }
  }
  #endif


  // exchange locally... try and hide within Isend latency... 
  if(tid>0){
  if(level->exchange_ghosts[shape].num_blocks[1]){
    if(tid==1)_timeStart = getTime();
    int lo = ((tid-1)*level->exchange_ghosts[shape].num_blocks[1])/(nth-1);
    int hi = ((tid  )*level->exchange_ghosts[shape].num_blocks[1])/(nth-1);
    if(hi>level->exchange_ghosts[shape].num_blocks[1])hi=level->exchange_ghosts[shape].num_blocks[1];
    for(buffer=lo;buffer<hi;buffer++){
      CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[1][buffer]);
    }
    if(tid==1)level->timers.ghostZone_local += (getTime()-_timeStart);
  }}


  // wait for MPI to finish...
  #ifdef USE_MPI 
  if(tid==0){
  if(nMessages){
    _timeStart = getTime();
    MPI_Waitall(nMessages,level->exchange_ghosts[shape].requests,level->exchange_ghosts[shape].status);
    level->timers.ghostZone_wait += (getTime()-_timeStart);
  }
  }

  #pragma omp barrier // wait for master to WaitAll

  // unpack MPI receive buffers 
  if(level->exchange_ghosts[shape].num_blocks[2]){
    if(tid==0)_timeStart = getTime();
    #pragma omp for private(buffer) schedule(static,1)
    for(buffer=0;buffer<level->exchange_ghosts[shape].num_blocks[2];buffer++){
      CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[2][buffer]);
    }
    if(tid==0)level->timers.ghostZone_unpack += (getTime()-_timeStart);
  }
  #endif

  } 
  level->timers.ghostZone_total += (double)(getTime()-_timeCommunicationStart);
}
