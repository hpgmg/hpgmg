//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// perform a (intra-level) ghost zone exchange
//  NOTE exchange_boundary() only exchanges the boundary.  
//  It will not enforce any boundary conditions
//  BC's are either the responsibility of a separate function or should be fused into the stencil
//------------------------------------------------------------------------------------------------------------------------------
//#define MAX_MESSAGE_SIZE 262144  // 256^3
//#define MAX_MESSAGE_SIZE 131072  // 
//#define MAX_MESSAGE_SIZE  65536  // 256^3
//#define MAX_MESSAGE_SIZE  32768  //
//#define MAX_MESSAGE_SIZE  16384  // 128^3
//#define MAX_MESSAGE_SIZE   8192  //
#define MAX_MESSAGE_SIZE   4096  //  64^3
//#define MAX_MESSAGE_SIZE   2048  //
//#define MAX_MESSAGE_SIZE   1024  //  32^3
//------------------------------------------------------------------------------------------------------------------------------
void exchange_boundary(level_type * level, int id, int justFaces){
  uint64_t _timeCommunicationStart = CycleTime();
  uint64_t _timeStart,_timeEnd;
  int buffer=0;
  int n,nn,ofs;

  if(justFaces)justFaces=1;else justFaces=0;  // must be 0 or 1 in order to index into exchange_ghosts[]

  #ifdef USE_MPI
  // count the actual number of messages...
  int nMessages = 0;
  for(n=0;n<level->exchange_ghosts[justFaces].num_recvs;n++){nMessages+=(level->exchange_ghosts[justFaces].recv_sizes[n]+MAX_MESSAGE_SIZE-1)/MAX_MESSAGE_SIZE;}

  // allocate temporary arrays...
  MPI_Request *recv_requests = NULL;MPI_Request *send_requests = NULL;
  MPI_Status  *recv_status   = NULL;MPI_Status  *send_status   = NULL;
       double**recv_buffers  = NULL;     double**send_buffers  = NULL;
          int *recv_ranks    = NULL;        int *send_ranks    = NULL;
          int *recv_sizes    = NULL;        int *send_sizes    = NULL;
          int *recv_tags     = NULL;        int *send_tags     = NULL;
  if(nMessages){
    recv_requests = (MPI_Request*)malloc(2*nMessages*sizeof(MPI_Request));send_requests = recv_requests + nMessages;
    recv_status   = (MPI_Status* )malloc(2*nMessages*sizeof(MPI_Status ));send_status   = recv_status   + nMessages;
    recv_buffers  = (double**    )malloc(2*nMessages*sizeof(double*    ));send_buffers  = recv_buffers  + nMessages;
    recv_ranks    = (int*        )malloc(6*nMessages*sizeof(int        ));send_ranks    = recv_ranks  + 1*nMessages;
                                                                          recv_sizes    = recv_ranks  + 2*nMessages;
                                                                          send_sizes    = recv_ranks  + 3*nMessages;
                                                                          recv_tags     = recv_ranks  + 4*nMessages;
                                                                          send_tags     = recv_ranks  + 5*nMessages;
  }

  // populate temoprary arrays...
  for(nn=0,n=0;n<level->exchange_ghosts[justFaces].num_recvs;n++){
  for(ofs=0;ofs<level->exchange_ghosts[justFaces].recv_sizes[n];ofs+=MAX_MESSAGE_SIZE){
    recv_buffers[nn] = level->exchange_ghosts[justFaces].recv_buffers[n]+ofs;
      recv_sizes[nn] = level->exchange_ghosts[justFaces].recv_sizes[n]-ofs;if(recv_sizes[nn]>MAX_MESSAGE_SIZE)recv_sizes[nn]=MAX_MESSAGE_SIZE;
      recv_ranks[nn] = level->exchange_ghosts[justFaces].recv_ranks[n];
       recv_tags[nn] = ofs;
                 nn++;
  }}
  for(nn=0,n=0;n<level->exchange_ghosts[justFaces].num_sends;n++){
  for(ofs=0;ofs<level->exchange_ghosts[justFaces].send_sizes[n];ofs+=MAX_MESSAGE_SIZE){
    send_buffers[nn] = level->exchange_ghosts[justFaces].send_buffers[n]+ofs;
      send_sizes[nn] = level->exchange_ghosts[justFaces].send_sizes[n]-ofs;if(send_sizes[nn]>MAX_MESSAGE_SIZE)send_sizes[nn]=MAX_MESSAGE_SIZE;
      send_ranks[nn] = level->exchange_ghosts[justFaces].send_ranks[n];
       send_tags[nn] = ofs;
                 nn++;
  }}

  // loop through packed list of MPI receives and prepost Irecv's...
  _timeStart = CycleTime();
  for(n=0;n<nMessages;n++){MPI_Irecv(recv_buffers[n],recv_sizes[n],MPI_DOUBLE,recv_ranks[n],recv_tags[n],MPI_COMM_WORLD,&recv_requests[n]);}
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
  for(n=0;n<nMessages;n++){MPI_Isend(send_buffers[n],send_sizes[n],MPI_DOUBLE,send_ranks[n],send_tags[n],MPI_COMM_WORLD,&send_requests[n]);}
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
  if(2*nMessages)MPI_Waitall(2*nMessages,recv_requests,recv_status); // send_requests/status are simply appended to recv_requests/status
  _timeEnd = CycleTime();
  level->cycles.ghostZone_wait += (_timeEnd-_timeStart);


  // unpack MPI receive buffers 
  _timeStart = CycleTime();
  #pragma omp parallel for if(level->exchange_ghosts[justFaces].num_blocks[2]>1) schedule(static,1)
  for(buffer=0;buffer<level->exchange_ghosts[justFaces].num_blocks[2];buffer++){CopyBlock(level,id,&level->exchange_ghosts[justFaces].blocks[2][buffer]);}
  _timeEnd = CycleTime();
  level->cycles.ghostZone_unpack += (_timeEnd-_timeStart);


  // free buffers
  if(nMessages){
  free(recv_ranks);
  free(recv_buffers);
  free(recv_status);
  free(recv_requests);
  }

  #endif


  level->cycles.ghostZone_total += (uint64_t)(CycleTime()-_timeCommunicationStart);
}
