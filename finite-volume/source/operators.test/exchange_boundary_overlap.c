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

  if(shape>=STENCIL_MAX_SHAPES)shape=STENCIL_SHAPE_BOX;  // shape must be < STENCIL_MAX_SHAPES in order to safely index into exchange_ghosts[]
  int my_tag = (level->tag<<4) | shape;

  // short circuit if no MPI to do...
  if( (level->exchange_ghosts[shape].num_blocks[0] + level->exchange_ghosts[shape].num_blocks[1] + level->exchange_ghosts[shape].num_blocks[2]) == 0)return;

  #ifdef USE_MPI
  int nMessages = level->exchange_ghosts[shape].num_recvs + level->exchange_ghosts[shape].num_sends;
  MPI_Request *recv_requests = level->exchange_ghosts[shape].requests;
  MPI_Request *send_requests = level->exchange_ghosts[shape].requests + level->exchange_ghosts[shape].num_recvs;
  #endif

  #ifdef _OPENMP
  #warning exchange_boundary_overlap.c must be run with at least 2 threads per process
  #pragma omp parallel
  {
    double _timeStart;
    int buffer=0;
    int n;
    int threadID   = omp_get_thread_num();
    //int numThreads = omp_get_num_threads();
    int numThreads = level->num_threads;

    #ifdef USE_MPI
    if(threadID==0){
    // loop through packed list of MPI receives and prepost Irecv's...
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
    } // master thread
    #endif


    // n-1 threads race thru to pack MPI send buffers... overlap packing with irecv
    #ifdef USE_MPI
    if(threadID>0){
    if(level->exchange_ghosts[shape].num_blocks[0]){
      if(threadID==1)_timeStart = getTime();
      for(buffer=threadID-1;buffer<level->exchange_ghosts[shape].num_blocks[0];buffer+=(numThreads-1)){ // like schedule(static,1) with n-1 threads
        CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[0][buffer]);
      }
      if(threadID==1)level->timers.ghostZone_pack += (getTime()-_timeStart);
    }
    }
    #endif


    #pragma omp barrier // master must wait for n-1 threads to complete MPI buffer packing


    #ifdef USE_MPI
    if(threadID==0){
    // master does MPI while other threads race thru
    // loop through MPI send buffers and post Isend's...
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
    // master waits for MPI to finish...
    if(nMessages){
      _timeStart = getTime();
      MPI_Waitall(nMessages,level->exchange_ghosts[shape].requests,level->exchange_ghosts[shape].status);
      level->timers.ghostZone_wait += (getTime()-_timeStart);
    }
    } // master
    #endif


    // exchange locally using n-1 threads... overlap local exchange with isend/waitall
    if(threadID>0){
    if(level->exchange_ghosts[shape].num_blocks[1]){
      if(threadID==1)_timeStart = getTime();
      for(buffer=threadID-1;buffer<level->exchange_ghosts[shape].num_blocks[1];buffer+=(numThreads-1)){ // like schedule(static,1) with n-1 threads
        CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[1][buffer]);
      }
      if(threadID==1)level->timers.ghostZone_local += (getTime()-_timeStart);
    }
    }


    #pragma omp barrier // all threads must wait for MPI to complete before unpacking


    // unpack MPI receive buffers
    #ifdef USE_MPI
    if(level->exchange_ghosts[shape].num_blocks[2]){
      if(threadID==0)_timeStart = getTime();
      //can't use PRAGMA_THREAD_ACROSS_BLOCKS as it would create a nested parallel region
      //#pragma omp for schedule(static,1)
      //for(buffer=0;buffer<level->exchange_ghosts[shape].num_blocks[2];buffer++){
      for(buffer=threadID;buffer<level->exchange_ghosts[shape].num_blocks[2];buffer+=numThreads){ // like schedule(static,1) with n threads
        CopyBlock(level,id,&level->exchange_ghosts[shape].blocks[2][buffer]);
      }
      if(threadID==0)level->timers.ghostZone_unpack += (getTime()-_timeStart);
    }
    #endif
  } // omp parallel
  #else
  #error exchange_boundary_overlap.c must be compiled with OpenMP and run with at least 2 threads per process
  #endif
  level->timers.ghostZone_total += (double)(getTime()-_timeCommunicationStart);
}
