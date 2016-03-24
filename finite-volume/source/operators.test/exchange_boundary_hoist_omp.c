//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// perform a (intra-level) ghost zone exchange on vector id
//  NOTE exchange_boundary() only exchanges the boundary.  
//  It will not enforce any boundary conditions
//  BC's are either the responsibility of a separate function or should be fused into the stencil
// The argument shape indicates which of faces, edges, and corners on each box must be exchanged
//  If the specified shape exceeds the range of defined shapes, the code will default to STENCIL_SHAPE_BOX (i.e. exchange faces, edges, and corners)
void exchange_boundary(level_type * level, int id, int shape){
  double _timeCommunicationStart = getTime();

  if( level->exchange_ghosts[shape].num_blocks[0] ||
      level->exchange_ghosts[shape].num_blocks[1] ||
      level->exchange_ghosts[shape].num_blocks[2] )
  {
  #pragma omp parallel
  {
  int threadID = omp_get_thread_num();

  double _timeStart;

  if(shape>=STENCIL_MAX_SHAPES)shape=STENCIL_SHAPE_BOX;  // shape must be < STENCIL_MAX_SHAPES in order to safely index into exchange_ghosts[]
  int my_tag = (level->tag<<4) | shape;
  int buffer=0;
  int n;

  #ifdef USE_MPI
  int nMessages = level->exchange_ghosts[shape].num_recvs + level->exchange_ghosts[shape].num_sends;
  MPI_Request *recv_requests = level->exchange_ghosts[shape].requests;
  MPI_Request *send_requests = level->exchange_ghosts[shape].requests + level->exchange_ghosts[shape].num_recvs;

  // loop through packed list of MPI receives and prepost Irecv's...
  if(level->exchange_ghosts[shape].num_recvs>0){
  if(threadID==0){
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
  }}


  // pack MPI send buffers...
  if(level->exchange_ghosts[shape].num_blocks[0]){
    if(threadID==1)_timeStart = getTime();
    if(threadID>0) // let thread 0 keep churning on Irecv's
    for(buffer=threadID-1;buffer<level->exchange_ghosts[shape].num_blocks[0];buffer+=(level->num_threads-1)){
    //#pragma omp for schedule(static,1)
    //for(buffer=0;buffer<level->exchange_ghosts[shape].num_blocks[0];buffer++){
      CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[0][buffer]);
    }
    if(threadID==1)level->timers.ghostZone_pack += (getTime()-_timeStart);
  }

 
  // loop through MPI send buffers and post Isend's...
  if(level->exchange_ghosts[shape].num_sends>0){
    #pragma omp barrier // wait for threads to finish packing...
    if(threadID==0){
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
  if(level->exchange_ghosts[shape].num_blocks[1]){
    if(threadID==1)_timeStart = getTime();
    if(threadID>0) 
    for(buffer=threadID-1;buffer<level->exchange_ghosts[shape].num_blocks[1];buffer+=(level->num_threads-1)){
    //#pragma omp for schedule(static,1)
    //for(buffer=0;buffer<level->exchange_ghosts[shape].num_blocks[1];buffer++){
      CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[1][buffer]);
    }
    if(threadID==1)level->timers.ghostZone_local += (getTime()-_timeStart);
  }


  // wait for MPI to finish...
  #ifdef USE_MPI 
  if(nMessages){
  if(threadID==0){
    _timeStart = getTime();
    MPI_Waitall(nMessages,level->exchange_ghosts[shape].requests,level->exchange_ghosts[shape].status);
    level->timers.ghostZone_wait += (getTime()-_timeStart);
  }}


  // unpack MPI receive buffers 
  if(level->exchange_ghosts[shape].num_blocks[2]){
    #pragma omp barrier // wait for thread 0 to finish WaitAll
    if(threadID==0)_timeStart = getTime();
    #pragma omp for schedule(static,1)
    for(buffer=0;buffer<level->exchange_ghosts[shape].num_blocks[2];buffer++){
      CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[2][buffer]);
    }
    if(threadID==0)level->timers.ghostZone_unpack += (getTime()-_timeStart);
  }
  #endif

  }} 
  level->timers.ghostZone_total += (double)(getTime()-_timeCommunicationStart);
}
