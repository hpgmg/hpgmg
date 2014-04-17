//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <omp.h>
//------------------------------------------------------------------------------------------------------------------------------
#include "level.h"
//#include "operators.h"
//------------------------------------------------------------------------------------------------------------------------------
//
//    / 24 25 26 /
//   / 21 22 23 /
//  / 18 19 20 /
//
//    / 15 16 17 /
//   / 12 13 14 /
//  /  9 10 11 /
//
//    /  6  7  8 /
//   /  3  4  5 /
//  /  0  1  2 /
//
int    faces[27] = {0,0,0,0,1,0,0,0,0,  0,1,0,1,0,1,0,1,0,  0,0,0,0,1,0,0,0,0};
int    edges[27] = {0,1,0,1,0,1,0,1,0,  1,0,1,0,0,0,1,0,1,  0,1,0,1,0,1,0,1,0};
int  corners[27] = {1,0,1,0,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,  1,0,1,0,0,0,1,0,1};


//------------------------------------------------------------------------------------------------------------------------------
print_communicator(int printSendRecv, int rank, int level, communicator_type *comm){
  int i;
  printf("rank=%2d level=%d ",rank,level);
  if(printSendRecv & 0x1){
  printf("num_sends=%2d ",comm->num_sends);
  printf("send_ranks=[ ");for(i=0;i<comm->num_sends;i++)printf("%2d ",comm->send_ranks[i]);printf("] ");
  printf("send_sizes=[ ");for(i=0;i<comm->num_sends;i++)printf("%2d ",comm->send_sizes[i]);printf("] ");
  printf("send_buffers=[ ");for(i=0;i<comm->num_sends;i++)printf("%08lx ",(uint64_t)comm->send_buffers[i]);printf("] ");
  for(i=0;i<comm->num_blocks[0];i++)printf("[ %dx%dx%d from %d %d %d %d %d to %d %d %d %d %d ] ",comm->blocks[0][i].dim.i,comm->blocks[0][i].dim.j,comm->blocks[0][i].dim.k,comm->blocks[0][i].read.i,comm->blocks[0][i].read.j,comm->blocks[0][i].read.k,comm->blocks[0][i].read.jStride,comm->blocks[0][i].read.kStride,comm->blocks[0][i].write.i,comm->blocks[0][i].write.j,comm->blocks[0][i].write.k,comm->blocks[0][i].write.jStride,comm->blocks[0][i].write.kStride);
  printf("\n");
  }
  if(printSendRecv & 0x2){
  for(i=0;i<comm->num_blocks[1];i++)printf("[ %dx%dx%d from %d %d %d %d %d to %d %d %d %d %d ] ",comm->blocks[1][i].dim.i,comm->blocks[1][i].dim.j,comm->blocks[1][i].dim.k,comm->blocks[1][i].read.i,comm->blocks[1][i].read.j,comm->blocks[1][i].read.k,comm->blocks[1][i].read.jStride,comm->blocks[1][i].read.kStride,comm->blocks[1][i].write.i,comm->blocks[1][i].write.j,comm->blocks[1][i].write.k,comm->blocks[1][i].write.jStride,comm->blocks[1][i].write.kStride);
  printf("\n");
  }
  if(printSendRecv & 0x4){
  printf("num_recvs=%2d ",comm->num_recvs);
  printf("recv_ranks=[ ");for(i=0;i<comm->num_recvs;i++)printf("%2d ",comm->recv_ranks[i]);printf("] ");
  printf("recv_sizes=[ ");for(i=0;i<comm->num_recvs;i++)printf("%2d ",comm->recv_sizes[i]);printf("] ");
  printf("recv_buffers=[ ");for(i=0;i<comm->num_recvs;i++)printf("%08lx ",(uint64_t)comm->recv_buffers[i]);printf("] ");
  for(i=0;i<comm->num_blocks[2];i++)printf("[ %dx%dx%d from %d %d %d %d %d to %d %d %d %d %d ] ",comm->blocks[2][i].dim.i,comm->blocks[2][i].dim.j,comm->blocks[2][i].dim.k,comm->blocks[2][i].read.i,comm->blocks[2][i].read.j,comm->blocks[2][i].read.k,comm->blocks[2][i].read.jStride,comm->blocks[2][i].read.kStride,comm->blocks[2][i].write.i,comm->blocks[2][i].write.j,comm->blocks[2][i].write.k,comm->blocks[2][i].write.jStride,comm->blocks[2][i].write.kStride);
  printf("\n");
  }
  fflush(stdout);
}
//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  int sendRank;
  int sendBoxID;
  int sendBox;
  int sendDir;
  int recvRank;
  int recvBoxID;
  int recvBox;
} GZ_type;


int qsortGZ(const void *a, const void*b){
  GZ_type *gza = (GZ_type*)a;
  GZ_type *gzb = (GZ_type*)b;
  // by convention, MPI buffers are first sorted by sendRank
  if(gza->sendRank < gzb->sendRank)return(-1);
  if(gza->sendRank > gzb->sendRank)return( 1);
  // then by sendBoxID
  if(gza->sendBoxID < gzb->sendBoxID)return(-1);
  if(gza->sendBoxID > gzb->sendBoxID)return( 1);
  // and finally by the direction sent
  if(gza->sendDir < gzb->sendDir)return(-1);
  if(gza->sendDir > gzb->sendDir)return( 1);
  return(0);
}


int qsortInt(const void *a, const void *b){
  int *ia = (int*)a;
  int *ib = (int*)b;
  if(*ia < *ib)return(-1);
  if(*ia > *ib)return( 1);
               return( 0);
}

//------------------------------------------------------------------------------------------------------------------------------
int create_box(box_type *box, int numVectors, int dim, int ghosts){
  uint64_t memory_allocated = 0;
  box->numVectors = numVectors;
  box->dim = dim;
  box->ghosts = ghosts;
  #ifndef BOX_SIMD_ALIGNMENT
  #define BOX_SIMD_ALIGNMENT 1 // allignment requirement for j+/-1
  #endif
  box->jStride = (dim+2*ghosts);while(box->jStride % BOX_SIMD_ALIGNMENT)box->jStride++; // pencil
  #ifndef BOX_PLANE_PADDING
  #define BOX_PLANE_PADDING  8 // scratch space to avoid unrolled loop cleanup
  #endif
  box->kStride = box->jStride*(dim+2*ghosts); // plane
  if(box->jStride<BOX_PLANE_PADDING)box->kStride += (BOX_PLANE_PADDING-box->jStride); // allow the ghost zone to be clobbered...
  while(box->kStride % BOX_SIMD_ALIGNMENT)box->kStride++;

  #if 0
  // pad each plane such that 
  //   1. it is greater than a multiple of the maximum unrolling (or vector length)
  //   2. it carries the same alignement as the plane above/below.   i.e.  if ijk is aligned, ijkk+/-plane is aligned
//int MaxUnrolling =  8; // 2-way SIMD x unroll by 4 =  8/thread
  int MaxUnrolling = 16; // 4-way SIMD x unroll by 4 = 16/thread
//int MaxUnrolling = 32; // 8-way SIMD x unroll by 4 = 32/thread
  int paddingToAvoidStencilCleanup = 0;
  if(box->jStride < (MaxUnrolling-1)){paddingToAvoidStencilCleanup = (MaxUnrolling-1)-(box->jStride);} // allows partial clobbering of ghost zone
//box->kStride  =( ((dim+2*ghosts)*box->jStride)+paddingToAvoidStencilCleanup+0xF) & ~0xF; // multiple of  128 bytes
  box->kStride  =( ((dim+2*ghosts)*box->jStride)+paddingToAvoidStencilCleanup+0x7) & ~0x7; // multiple of   64 bytes (required for MIC)
//box->kStride  =( ((dim+2*ghosts)*box->jStride)+paddingToAvoidStencilCleanup+0x3) & ~0x3; // multiple of   32 bytes (required for AVX/QPX)
//box->kStride  =( ((dim+2*ghosts)*box->jStride)+paddingToAvoidStencilCleanup+0x1) & ~0x1; // multiple of   16 bytes (required for SSE)
  #endif
  box->volume = (dim+2*ghosts)*box->kStride;


  // allocate an array of pointers to vectors...
  box->vectors = (double **)malloc(box->numVectors*sizeof(double*));
               memory_allocated += box->numVectors*sizeof(double*);
  // allocate one aligned, double-precision array and divide it among vectors...
  uint64_t malloc_size = box->volume*box->numVectors*sizeof(double) + 4096; // shift pointer by up to 1 TLB page...
  box->vectors_base = (double*)malloc(malloc_size);
                  memory_allocated += malloc_size;
  double * tmpbuf = box->vectors_base;
  while( (uint64_t)(tmpbuf+box->ghosts*(1+box->jStride+box->kStride)) & (4096-1) ){tmpbuf++;} // allign first *non-ghost* zone element of first component to page boundary...
  memset(tmpbuf,0,box->volume*box->numVectors*sizeof(double)); // zero to avoid 0.0*NaN or 0.0*Inf
  int c;for(c=0;c<box->numVectors;c++){box->vectors[c] = tmpbuf + c*box->volume;}


  // done... return the total amount of memory allocated...
  return(memory_allocated);
}


//------------------------------------------------------------------------------------------------------------------------------
void add_vectors_to_box(box_type *box, int numAdditionalVectors){
  if(numAdditionalVectors<=0)return;									// nothing to do
  double * old_bp = box->vectors_base;									// save a pointer to the base pointer for subsequent free...
  double * old_v0 = box->vectors[0];									// save a pointer to the old FP data...
  double ** old_v = box->vectors;									// save a pointer to the old array of pointers...
  box->numVectors+=numAdditionalVectors;								//
  box->vectors = (double **)malloc(box->numVectors*sizeof(double*));				// new array of pointers vectors
  // NOTE !!! realloc() cannot guarantee the same alignment... malloc, allign, copy...
  uint64_t malloc_size = box->volume*box->numVectors*sizeof(double) + 4096; // shift pointer by up to 1 TLB page...
  box->vectors_base = (double*)malloc(malloc_size);
  double * tmpbuf = box->vectors_base;
  while( (uint64_t)(tmpbuf+box->ghosts*(1+box->jStride+box->kStride)) & (4096-1) ){tmpbuf++;} // allign first *non-ghost* zone element of first component to page boundary...
  memset(tmpbuf,0,box->volume*box->numVectors*sizeof(double));	// zero to avoid 0.0*NaN or 0.0*Inf
  memcpy(tmpbuf,old_v0,box->volume*(box->numVectors-numAdditionalVectors)*sizeof(double));	// copy any existant data over...
  int c;for(c=0;c<box->numVectors;c++){box->vectors[c] = tmpbuf + c*box->volume;}			// pointer arithmetic...
  free(old_bp);												// all vectors were created from one malloc + pointer arithmetic...
  free(old_v );												// free the list of pointers...
}


//------------------------------------------------------------------------------------------------------------------------------
void destroy_box(box_type *box){
  free(box->vectors_base); // single allocate with pointer arithmetic...
  free(box->vectors);
}


//------------------------------------------------------------------------------------------------------------------------------
// should implement a 3D hilbert curve on non pow2 (but cubical) domain sizes
//void decompose_level_hilbert(int *rank_of_box, int jStride, int kStride, int ilo, int jlo, int klo, int idim, int jdim, int kdim, int rank_lo, int ranks){
//}
//---------------------------------------------------------------------------------------------------------------------------------------------------
void decompose_level_kd_tree(int *rank_of_box, int jStride, int kStride, int ilo, int jlo, int klo, int idim, int jdim, int kdim, int rank_lo, int ranks){
  // recursive bisection (or prime-section) of the domain
  // can lead to imbalance unless the number of processes and number of boxes per process are chosen well

  #define numPrimes 13
  int primes[numPrimes] = {41,37,31,29,23,19,17,13,11,7,5,3,2};
  int i,j,k,p,f,ff;


  // base case, no further recursion...
  if( (ranks==1)|| ((idim==1)&&(jdim==1)&&(kdim==1)) ){
    for(i=ilo;i<ilo+idim;i++){
    for(j=jlo;j<jlo+jdim;j++){
    for(k=klo;k<klo+kdim;k++){
      int b = i + j*jStride + k*kStride;
      rank_of_box[b] = rank_lo;
    }}}
    return;
  }


  for(p=0;p<numPrimes;p++){f=primes[p];
    // don't partition if the aspect ratio would be extreme
    if( (1.5*kdim>=idim)&&(1.5*kdim>=jdim) )if( (kdim%f==0) && (ranks%f==0) ){for(ff=0;ff<f;ff++)decompose_level_kd_tree(rank_of_box,jStride,kStride,ilo,jlo,klo+ff*kdim/f,idim,jdim,kdim/f,rank_lo+ff*ranks/f,ranks/f);return;}
    if( (1.5*jdim>=idim)&&(1.5*jdim>=kdim) )if( (jdim%f==0) && (ranks%f==0) ){for(ff=0;ff<f;ff++)decompose_level_kd_tree(rank_of_box,jStride,kStride,ilo,jlo+ff*jdim/f,klo,idim,jdim/f,kdim,rank_lo+ff*ranks/f,ranks/f);return;}
    if( (1.5*idim>=jdim)&&(1.5*idim>=kdim) )if( (idim%f==0) && (ranks%f==0) ){for(ff=0;ff<f;ff++)decompose_level_kd_tree(rank_of_box,jStride,kStride,ilo+ff*idim/f,jlo,klo,idim/f,jdim,kdim,rank_lo+ff*ranks/f,ranks/f);return;}
  }


  // try and bisect the domain in the i-dimension
  if( (idim>=jdim)&&(idim>=kdim) ){
    int dim0 = (int)(0.5*(double)idim + 0.50);
    int dim1 = idim-dim0;
    int r0 = (int)( 0.5 + (double)ranks*(double)dim0/(double)idim );
    int r1 = ranks-r0;
    decompose_level_kd_tree(rank_of_box,jStride,kStride,ilo     ,jlo,klo,dim0,jdim,kdim,rank_lo   ,r0); // lo
    decompose_level_kd_tree(rank_of_box,jStride,kStride,ilo+dim0,jlo,klo,dim1,jdim,kdim,rank_lo+r0,r1); // hi
    return;
  }
  // try and bisect the domain in the j-dimension
  if( (jdim>=idim)&&(jdim>=kdim) ){
    int dim0 = (int)(0.5*(double)jdim + 0.50);
    int dim1 = jdim-dim0;
    int r0 = (int)( 0.5 + (double)ranks*(double)dim0/(double)jdim );
    int r1 = ranks-r0;
    decompose_level_kd_tree(rank_of_box,jStride,kStride,ilo,jlo     ,klo,idim,dim0,kdim,rank_lo   ,r0); // lo
    decompose_level_kd_tree(rank_of_box,jStride,kStride,ilo,jlo+dim0,klo,idim,dim1,kdim,rank_lo+r0,r1); // hi
    return;
  }
  // try and bisect the domain in the k-dimension
  if( (kdim>=idim)&&(kdim>=jdim) ){
    int dim0 = (int)(0.5*(double)kdim + 0.50);
    int dim1 = kdim-dim0;
    int r0 = (int)( 0.5 + (double)ranks*(double)dim0/(double)kdim );
    int r1 = ranks-r0;
    decompose_level_kd_tree(rank_of_box,jStride,kStride,ilo,jlo,klo     ,idim,jdim,dim0,rank_lo   ,r0); // lo
    decompose_level_kd_tree(rank_of_box,jStride,kStride,ilo,jlo,klo+dim0,idim,jdim,dim1,rank_lo+r0,r1); // hi
    return;
  }
    printf("decompose_level_kd_tree failed !!!\n");
    exit(0);
}


//---------------------------------------------------------------------------------------------------------------------------------------------------
void print_decomposition(level_type *level){
  if(level->my_rank!=0)return;
  printf("\n");
  int i,j,k;
  int jStride = level->boxes_in.i;
  int kStride = level->boxes_in.i*level->boxes_in.j;
  for(k=level->boxes_in.k-1;k>=0;k--){ // (i,j,k)=(0,0,0) is bottom left corner
  for(j=level->boxes_in.j-1;j>=0;j--){ // (i,j)=(0,0) is bottom left corner
  for(i=0;i<j;i++)printf(" ");
  for(i=0;i<level->boxes_in.i;i++){
    int b = i + j*jStride + k*kStride;
    printf("%4d ",level->rank_of_box[b]);
  }printf("\n");
  }printf("\n\n");
  }
  fflush(stdout);
}


//------------------------------------------------------------------------------------------------------------------------------
#ifndef BLOCKCOPY_TILE_J
#define BLOCKCOPY_TILE_J 8
#endif
#ifndef BLOCKCOPY_TILE_K
#define BLOCKCOPY_TILE_K 8
#endif
void append_block_to_list(blockCopy_type * blocks, int *tail, int doWrite,
                          int dim_i, int dim_j, int dim_k,
                          int  read_box, double*  read_ptr, int  read_i, int  read_j, int  read_k, int  read_jStride, int  read_kStride,
                          int write_box, double* write_ptr, int write_i, int write_j, int write_k, int write_jStride, int write_kStride
                         ){
  int jj,kk;
  // Take a dim_j x dim_k iteration space and tile it into smaller faces of size BLOCKCOPY_TILE_J x BLOCKCOPY_TILE_K
  // This increases the number of blockCopies in the ghost zone exchange and thereby increases the thread-level parallelism
  for(kk=0;kk<dim_k;kk+=BLOCKCOPY_TILE_K){
  for(jj=0;jj<dim_j;jj+=BLOCKCOPY_TILE_J){
    int dim_k_mod = dim_k-kk;if(dim_k_mod>BLOCKCOPY_TILE_K)dim_k_mod=BLOCKCOPY_TILE_K;
    int dim_j_mod = dim_j-jj;if(dim_j_mod>BLOCKCOPY_TILE_J)dim_j_mod=BLOCKCOPY_TILE_J;
    if(doWrite){
      blocks[*tail].dim.i         = dim_i;
      blocks[*tail].dim.j         = dim_j_mod;
      blocks[*tail].dim.k         = dim_k_mod;
      blocks[*tail].read.box      = read_box;
      blocks[*tail].read.ptr      = read_ptr;
      blocks[*tail].read.i        = read_i;
      blocks[*tail].read.j        = read_j + jj;
      blocks[*tail].read.k        = read_k + kk;
      blocks[*tail].read.jStride  = read_jStride;
      blocks[*tail].read.kStride  = read_kStride;
      blocks[*tail].write.box     = write_box;
      blocks[*tail].write.ptr     = write_ptr;
      blocks[*tail].write.i       = write_i;
      blocks[*tail].write.j       = write_j + jj;
      blocks[*tail].write.k       = write_k + kk;
      blocks[*tail].write.jStride = write_jStride;
      blocks[*tail].write.kStride = write_kStride;
    }       (*tail)++;
  }}
}

//----------------------------------------------------------------------------------------------------------------------------------------------------
// create a mini program that packs data into MPI recv buffers, exchanges local data, and unpacks the MPI send buffers
//   broadly speaking... 
//   1. traverse my list of Boxes and create a list of ghosts that must be sent
//   2. create a list of neighbors to send to
//   3. allocate and populate the pack list and allocate the send buffers
//   4. allocate and populate the local exchange list
//   5. traverse my list of Boxes and create a list of ghosts that must be received
//   6. create a list of neighbors to receive from
//   7. allocate and populate the unpack list and allocate the recv buffers
//
//   thus a ghost zone exchange is
//   1. prepost a Irecv for each MPI recv buffer (1 per neighbor)
//   2. traverse the pack list
//   3. post the Isends for each MPI send buffer (1 per neighbor)
//   4. traverse the local copy list
//   5. waitall
//   6. traverse the unpack list
void build_exchange_ghosts(level_type *level, int justFaces){
  level->exchange_ghosts[justFaces].num_recvs     = 0;
  level->exchange_ghosts[justFaces].num_sends     = 0;
  level->exchange_ghosts[justFaces].num_blocks[0] = 0;
  level->exchange_ghosts[justFaces].num_blocks[1] = 0;
  level->exchange_ghosts[justFaces].num_blocks[2] = 0;

  int CommunicateThisDir[27];
         int n;for(n=0;n<27;n++)CommunicateThisDir[n]=1;CommunicateThisDir[13]=0;
  if(justFaces)for(n=0;n<27;n++)CommunicateThisDir[n]=faces[n];

  int sendBox,recvBox;
  int stage;
  int _rank;
  int ghost,numGhosts,numGhostsRemote;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // traverse my list of boxes and create a lists of neighboring boxes and neighboring ranks
  GZ_type *ghostsToSend = (GZ_type*)malloc(26*level->num_my_boxes*sizeof(GZ_type)); // There are at most 26 neighbors per box.
         int *sendRanks = (    int*)malloc(26*level->num_my_boxes*sizeof(    int)); // There are at most 26 neighbors per box.
  numGhosts       = 0;
  numGhostsRemote = 0;
  for(sendBox=0;sendBox<level->num_my_boxes;sendBox++){
    int di,dj,dk;
    for(dk=-1;dk<=1;dk++){
    for(dj=-1;dj<=1;dj++){
    for(di=-1;di<=1;di++){
      int dir = 13+di+3*dj+9*dk;if(CommunicateThisDir[dir]){
      int       myBoxID = level->my_boxes[sendBox].global_box_id;
      int       myBox_i = level->my_boxes[sendBox].low.i / level->box_dim;
      int       myBox_j = level->my_boxes[sendBox].low.j / level->box_dim;
      int       myBox_k = level->my_boxes[sendBox].low.k / level->box_dim;
      int neighborBox_i = (  myBox_i + di + level->boxes_in.i) % level->boxes_in.i;
      int neighborBox_j = (  myBox_j + dj + level->boxes_in.j) % level->boxes_in.j;
      int neighborBox_k = (  myBox_k + dk + level->boxes_in.k) % level->boxes_in.k;
      int neighborBoxID =  neighborBox_i + neighborBox_j*level->boxes_in.i + neighborBox_k*level->boxes_in.i*level->boxes_in.j;
      if( level->rank_of_box[neighborBoxID] != -1 ){
        ghostsToSend[numGhosts].sendRank  = level->my_rank;
        ghostsToSend[numGhosts].sendBoxID = myBoxID;
        ghostsToSend[numGhosts].sendBox   = sendBox;
        ghostsToSend[numGhosts].sendDir   = dir;
        ghostsToSend[numGhosts].recvRank  = level->rank_of_box[neighborBoxID];
        ghostsToSend[numGhosts].recvBoxID = neighborBoxID;
        ghostsToSend[numGhosts].recvBox   = -1;
        if( level->rank_of_box[neighborBoxID] != level->my_rank )sendRanks[numGhostsRemote++] = level->rank_of_box[neighborBoxID];else{
          int recvBox=0;while(level->my_boxes[recvBox].global_box_id!=neighborBoxID)recvBox++; // search my list of boxes for the appropriate recvBox index
          ghostsToSend[numGhosts].recvBox   = recvBox;
        }
        numGhosts++;
      }
    }}}}
  }
  // sort boxes by sendRank(==my rank) then by sendBoxID... ensures the sends and receive buffers are always sorted by sendBoxID...
  qsort(ghostsToSend,numGhosts      ,sizeof(GZ_type),qsortGZ );
  // sort the lists of neighboring ranks and remove duplicates...
  qsort(sendRanks   ,numGhostsRemote,sizeof(    int),qsortInt);
  int numSendRanks=0;_rank=-1;for(ghost=0;ghost<numGhostsRemote;ghost++)if(sendRanks[ghost] != _rank){_rank=sendRanks[ghost];sendRanks[numSendRanks++]=sendRanks[ghost];}


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in a two-stage process, traverse the list of ghosts and allocate the pack/local lists as well as the MPI buffers, and then populate the pack/local lists
  level->exchange_ghosts[justFaces].num_sends     =                  numSendRanks;
  level->exchange_ghosts[justFaces].send_ranks    =     (int*)malloc(numSendRanks*sizeof(int));
  level->exchange_ghosts[justFaces].send_sizes    =     (int*)malloc(numSendRanks*sizeof(int));
  level->exchange_ghosts[justFaces].send_buffers  = (double**)malloc(numSendRanks*sizeof(double*));
  for(stage=0;stage<=1;stage++){
    // stage=0... traverse the list and calculate the buffer sizes and number of pack/local blockCopies
    // stage=1... allocate structures, traverse the list, and populate the unpack/local lists...
    int neighbor;
    if(stage==1){
      for(neighbor=0;neighbor<numSendRanks;neighbor++){
             level->exchange_ghosts[justFaces].send_buffers[neighbor] = (double*)malloc(level->exchange_ghosts[justFaces].send_sizes[neighbor]*sizeof(double));
      memset(level->exchange_ghosts[justFaces].send_buffers[neighbor],                0,level->exchange_ghosts[justFaces].send_sizes[neighbor]*sizeof(double));
      }
      level->exchange_ghosts[justFaces].blocks[0]      = (blockCopy_type*)malloc(level->exchange_ghosts[justFaces].num_blocks[0]*sizeof(blockCopy_type));;
      level->exchange_ghosts[justFaces].blocks[1]      = (blockCopy_type*)malloc(level->exchange_ghosts[justFaces].num_blocks[1]*sizeof(blockCopy_type));;
    }
    for(neighbor=0;neighbor<numSendRanks;neighbor++){
      level->exchange_ghosts[justFaces].send_ranks[neighbor]=sendRanks[neighbor];
      level->exchange_ghosts[justFaces].send_sizes[neighbor]=0;
    }
    level->exchange_ghosts[justFaces].num_blocks[0] = 0;
    level->exchange_ghosts[justFaces].num_blocks[1] = 0;
    for(ghost=0;ghost<numGhosts;ghost++){
      int dim_i,dim_j,dim_k;
      int send_i,send_j,send_k;
      int recv_i,recv_j,recv_k;
  
      // decode ghostsToSend[ghost].sendDir (direction sent) into di/dj/dk 
      int di = ((ghostsToSend[ghost].sendDir % 3)  )-1;
      int dj = ((ghostsToSend[ghost].sendDir % 9)/3)-1;
      int dk = ((ghostsToSend[ghost].sendDir / 9)  )-1;
      switch(di){ // direction relative to sender
        case -1:send_i=0;                               dim_i=level->box_ghosts;recv_i=  level->box_dim;   break;
        case  0:send_i=0;                               dim_i=level->box_dim;   recv_i=0;                  break;
        case  1:send_i=level->box_dim-level->box_ghosts;dim_i=level->box_ghosts;recv_i=0-level->box_ghosts;break;
      }
      switch(dj){ // direction relative to sender
        case -1:send_j=0;                               dim_j=level->box_ghosts;recv_j=  level->box_dim;   break;
        case  0:send_j=0;                               dim_j=level->box_dim;   recv_j=0;                  break;
        case  1:send_j=level->box_dim-level->box_ghosts;dim_j=level->box_ghosts;recv_j=0-level->box_ghosts;break;
      }
      switch(dk){ // direction relative to sender
        case -1:send_k=0;                               dim_k=level->box_ghosts;recv_k=  level->box_dim;   break;
        case  0:send_k=0;                               dim_k=level->box_dim;   recv_k=0;                  break;
        case  1:send_k=level->box_dim-level->box_ghosts;dim_k=level->box_ghosts;recv_k=0-level->box_ghosts;break;
      }
 
      // determine if this ghost requires a pack or local exchange 
      int LocalExchange; // 0 = pack list, 1 = local exchange list
      if(ghostsToSend[ghost].recvRank != level->my_rank){
        LocalExchange=0; // pack
        neighbor=0;while(level->exchange_ghosts[justFaces].send_ranks[neighbor] != ghostsToSend[ghost].recvRank)neighbor++;
      }else{
        LocalExchange=1; // local
        neighbor=-1;
      }
    
      #if 0 
      if(stage==1){ // lists have been allocated... populate them...
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].dim.i         = dim_i;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].dim.j         = dim_j;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].dim.k         = dim_k;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].read.box      = ghostsToSend[ghost].sendBox;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].read.ptr      = NULL;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].read.i        = send_i;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].read.j        = send_j;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].read.k        = send_k;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].read.jStride  = level->my_boxes[ghostsToSend[ghost].sendBox].jStride;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].read.kStride  = level->my_boxes[ghostsToSend[ghost].sendBox].kStride;
        if(LocalExchange==0){
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.box     = -1;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.ptr     = level->exchange_ghosts[justFaces].send_buffers[neighbor];
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.i       = level->exchange_ghosts[justFaces].send_sizes[neighbor]; // current offset in the MPI send buffer
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.j       = 0;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.k       = 0;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.jStride = dim_i;       // contiguous block
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.kStride = dim_i*dim_j; // contiguous block
        }else{
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.box     = ghostsToSend[ghost].recvBox;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.ptr     = NULL;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.i       = recv_i;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.j       = recv_j;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.k       = recv_k;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.jStride = level->my_boxes[ghostsToSend[ghost].recvBox].jStride;
        level->exchange_ghosts[justFaces].blocks[LocalExchange][level->exchange_ghosts[justFaces].num_blocks[LocalExchange]].write.kStride = level->my_boxes[ghostsToSend[ghost].recvBox].kStride;
        }
      } 
                     level->exchange_ghosts[justFaces].num_blocks[LocalExchange]++;
      if(neighbor>=0)level->exchange_ghosts[justFaces].send_sizes[neighbor]+=dim_i*dim_j*dim_k;
      #else
      if(LocalExchange) // append to the local exchange list...
      append_block_to_list(&(level->exchange_ghosts[justFaces].blocks[1][0]),&(level->exchange_ghosts[justFaces].num_blocks[1]),(stage==1),
        /* dim.i         = */ dim_i,
        /* dim.j         = */ dim_j,
        /* dim.k         = */ dim_k,
        /* read.box      = */ ghostsToSend[ghost].sendBox,
        /* read.ptr      = */ NULL,
        /* read.i        = */ send_i,
        /* read.j        = */ send_j,
        /* read.k        = */ send_k,
        /* read.jStride  = */ level->my_boxes[ghostsToSend[ghost].sendBox].jStride,
        /* read.kStride  = */ level->my_boxes[ghostsToSend[ghost].sendBox].kStride,
        /* write.box     = */ ghostsToSend[ghost].recvBox,
        /* write.ptr     = */ NULL,
        /* write.i       = */ recv_i,
        /* write.j       = */ recv_j,
        /* write.k       = */ recv_k,
        /* write.jStride = */ level->my_boxes[ghostsToSend[ghost].recvBox].jStride,
        /* write.kStride = */ level->my_boxes[ghostsToSend[ghost].recvBox].kStride
      );else // append to the MPI pack list...
      append_block_to_list(&(level->exchange_ghosts[justFaces].blocks[0][0]),&(level->exchange_ghosts[justFaces].num_blocks[0]),(stage==1),
        /* dim.i         = */ dim_i,
        /* dim.j         = */ dim_j,
        /* dim.k         = */ dim_k,
        /* read.box      = */ ghostsToSend[ghost].sendBox,
        /* read.ptr      = */ NULL,
        /* read.i        = */ send_i,
        /* read.j        = */ send_j,
        /* read.k        = */ send_k,
        /* read.jStride  = */ level->my_boxes[ghostsToSend[ghost].sendBox].jStride,
        /* read.kStride  = */ level->my_boxes[ghostsToSend[ghost].sendBox].kStride,
        /* write.box     = */ -1,
        /* write.ptr     = */ level->exchange_ghosts[justFaces].send_buffers[neighbor],
        /* write.i       = */ level->exchange_ghosts[justFaces].send_sizes[neighbor], // current offset in the MPI send buffer
        /* write.j       = */ 0,
        /* write.k       = */ 0,
        /* write.jStride = */ dim_i,       // contiguous block
        /* write.kStride = */ dim_i*dim_j  // contiguous block
      );
      if(neighbor>=0)level->exchange_ghosts[justFaces].send_sizes[neighbor]+=dim_i*dim_j*dim_k;
      #endif
    } // ghost for-loop
  } // stage for-loop


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // free temporary storage...
  free(ghostsToSend);
  free(sendRanks);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // traverse my list of boxes and create a lists of neighboring boxes and neighboring ranks
  GZ_type *ghostsToRecv = (GZ_type*)malloc(26*level->num_my_boxes*sizeof(GZ_type)); // There are at most 26 neighbors per box.
         int *recvRanks = (    int*)malloc(26*level->num_my_boxes*sizeof(    int)); // There are at most 26 neighbors per box.
  numGhosts       = 0;
  numGhostsRemote = 0;
  for(recvBox=0;recvBox<level->num_my_boxes;recvBox++){
    int di,dj,dk;
    for(dk=-1;dk<=1;dk++){
    for(dj=-1;dj<=1;dj++){
    for(di=-1;di<=1;di++){
      int dir = 13+di+3*dj+9*dk;if(CommunicateThisDir[dir]){
      int       myBoxID = level->my_boxes[recvBox].global_box_id;
      int       myBox_i = level->my_boxes[recvBox].low.i / level->box_dim;
      int       myBox_j = level->my_boxes[recvBox].low.j / level->box_dim;
      int       myBox_k = level->my_boxes[recvBox].low.k / level->box_dim;
      int neighborBox_i = (  myBox_i + di + level->boxes_in.i) % level->boxes_in.i;
      int neighborBox_j = (  myBox_j + dj + level->boxes_in.j) % level->boxes_in.j;
      int neighborBox_k = (  myBox_k + dk + level->boxes_in.k) % level->boxes_in.k;
      int neighborBoxID =  neighborBox_i + neighborBox_j*level->boxes_in.i + neighborBox_k*level->boxes_in.i*level->boxes_in.j;
      if( (level->rank_of_box[neighborBoxID] != -1) && (level->rank_of_box[neighborBoxID] != level->my_rank)  ){
        ghostsToRecv[numGhosts].sendRank  = level->rank_of_box[neighborBoxID];
        ghostsToRecv[numGhosts].sendBoxID = neighborBoxID;
        ghostsToRecv[numGhosts].sendBox   = -1;
        ghostsToRecv[numGhosts].sendDir   = 26-dir;
        ghostsToRecv[numGhosts].recvRank  = level->my_rank;
        ghostsToRecv[numGhosts].recvBoxID = myBoxID;
        ghostsToRecv[numGhosts].recvBox   = recvBox;
                     numGhosts++;
        recvRanks[numGhostsRemote++] = level->rank_of_box[neighborBoxID];
      }
    }}}}
  }
  // sort boxes by sendRank then by sendBoxID... ensures the recvs and receive buffers are always sorted by sendBoxID...
  qsort(ghostsToRecv,numGhosts      ,sizeof(GZ_type),qsortGZ );
  // sort the lists of neighboring ranks and remove duplicates...
  qsort(recvRanks   ,numGhostsRemote,sizeof(    int),qsortInt);
  int numRecvRanks=0;_rank=-1;for(ghost=0;ghost<numGhostsRemote;ghost++)if(recvRanks[ghost] != _rank){_rank=recvRanks[ghost];recvRanks[numRecvRanks++]=recvRanks[ghost];}


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // in a two-stage process, traverse the list of ghosts and allocate the unpack lists as well as the MPI buffers, and then populate the unpack list
  level->exchange_ghosts[justFaces].num_recvs     =                  numRecvRanks;
  level->exchange_ghosts[justFaces].recv_ranks    =     (int*)malloc(numRecvRanks*sizeof(int));
  level->exchange_ghosts[justFaces].recv_sizes    =     (int*)malloc(numRecvRanks*sizeof(int));
  level->exchange_ghosts[justFaces].recv_buffers  = (double**)malloc(numRecvRanks*sizeof(double*));
  for(stage=0;stage<=1;stage++){
    // stage=0... traverse the list and calculate the buffer sizes and number of pack/local blockCopies
    // stage=1... allocate structures, traverse the list, and populate the unpack/local lists...
    int neighbor;
    if(stage==1){
      for(neighbor=0;neighbor<numRecvRanks;neighbor++){
             level->exchange_ghosts[justFaces].recv_buffers[neighbor] = (double*)malloc(level->exchange_ghosts[justFaces].recv_sizes[neighbor]*sizeof(double));
      memset(level->exchange_ghosts[justFaces].recv_buffers[neighbor],                0,level->exchange_ghosts[justFaces].recv_sizes[neighbor]*sizeof(double));
      }
      level->exchange_ghosts[justFaces].blocks[2]      = (blockCopy_type*)malloc(level->exchange_ghosts[justFaces].num_blocks[2]*sizeof(blockCopy_type));;
    }
    for(neighbor=0;neighbor<numRecvRanks;neighbor++){
      level->exchange_ghosts[justFaces].recv_ranks[neighbor]=recvRanks[neighbor];
      level->exchange_ghosts[justFaces].recv_sizes[neighbor]=0;
    }
    level->exchange_ghosts[justFaces].num_blocks[2] = 0;
    for(ghost=0;ghost<numGhosts;ghost++){
      int dim_i,dim_j,dim_k;
      int send_i,send_j,send_k;
      int recv_i,recv_j,recv_k;
  
      // decode ghostsToRecv[ghost].sendDir (direction sent) into di/dj/dk 
      int di = ((ghostsToRecv[ghost].sendDir % 3)  )-1;
      int dj = ((ghostsToRecv[ghost].sendDir % 9)/3)-1;
      int dk = ((ghostsToRecv[ghost].sendDir / 9)  )-1;
      switch(di){ // direction relative to sender
        case -1:dim_i=level->box_ghosts;recv_i=  level->box_dim;   break;
        case  0:dim_i=level->box_dim;   recv_i=0;                  break;
        case  1:dim_i=level->box_ghosts;recv_i=0-level->box_ghosts;break;
      }
      switch(dj){ // direction relative to sender
        case -1:dim_j=level->box_ghosts;recv_j=  level->box_dim;   break;
        case  0:dim_j=level->box_dim;   recv_j=0;                  break;
        case  1:dim_j=level->box_ghosts;recv_j=0-level->box_ghosts;break;
      }
      switch(dk){ // direction relative to sender
        case -1:dim_k=level->box_ghosts;recv_k=  level->box_dim;   break;
        case  0:dim_k=level->box_dim;   recv_k=0;                  break;
        case  1:dim_k=level->box_ghosts;recv_k=0-level->box_ghosts;break;
      }
 
      // determine if this ghost requires a pack or local exchange 
      neighbor=0;while(level->exchange_ghosts[justFaces].recv_ranks[neighbor] != ghostsToRecv[ghost].sendRank)neighbor++;
    
      #if 0 
      if(stage==1){ // lists have been allocated... populate them...
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].dim.i         = dim_i;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].dim.j         = dim_j;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].dim.k         = dim_k;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].read.box      = -1;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].read.ptr      = level->exchange_ghosts[justFaces].recv_buffers[neighbor];
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].read.i        = level->exchange_ghosts[justFaces].recv_sizes[neighbor]; // current offset in the MPI recv buffer
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].read.j        = 0;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].read.k        = 0;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].read.jStride  = dim_i;       // contiguous block
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].read.kStride  = dim_i*dim_j; // contiguous block
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].write.box     = ghostsToRecv[ghost].recvBox;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].write.ptr     = NULL;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].write.i       = recv_i;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].write.j       = recv_j;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].write.k       = recv_k;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].write.jStride = level->my_boxes[ghostsToRecv[ghost].recvBox].jStride;
        level->exchange_ghosts[justFaces].blocks[2][level->exchange_ghosts[justFaces].num_blocks[2]].write.kStride = level->my_boxes[ghostsToRecv[ghost].recvBox].kStride;

      } 
                     level->exchange_ghosts[justFaces].num_blocks[2]++;
      if(neighbor>=0)level->exchange_ghosts[justFaces].recv_sizes[neighbor]+=dim_i*dim_j*dim_k;
      #else
      append_block_to_list(&(level->exchange_ghosts[justFaces].blocks[2][0]),&(level->exchange_ghosts[justFaces].num_blocks[2]),(stage==1),
      /*dim.i         = */ dim_i,
      /*dim.j         = */ dim_j,
      /*dim.k         = */ dim_k,
      /*read.box      = */ -1,
      /*read.ptr      = */ level->exchange_ghosts[justFaces].recv_buffers[neighbor],
      /*read.i        = */ level->exchange_ghosts[justFaces].recv_sizes[neighbor], // current offset in the MPI recv buffer
      /*read.j        = */ 0,
      /*read.k        = */ 0,
      /*read.jStride  = */ dim_i,       // contiguous block
      /*read.kStride  = */ dim_i*dim_j, // contiguous block
      /*write.box     = */ ghostsToRecv[ghost].recvBox,
      /*write.ptr     = */ NULL,
      /*write.i       = */ recv_i,
      /*write.j       = */ recv_j,
      /*write.k       = */ recv_k,
      /*write.jStride = */ level->my_boxes[ghostsToRecv[ghost].recvBox].jStride,
      /*write.kStride = */ level->my_boxes[ghostsToRecv[ghost].recvBox].kStride
      );
      if(neighbor>=0)level->exchange_ghosts[justFaces].recv_sizes[neighbor]+=dim_i*dim_j*dim_k;
      #endif
    } // ghost for-loop
  } // stage for-loop


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // free temporary storage...
  free(ghostsToRecv);
  free(recvRanks);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // malloc MPI requests/status arrays
  #ifdef USE_MPI
  level->exchange_ghosts[justFaces].requests = (MPI_Request*)malloc((level->exchange_ghosts[justFaces].num_sends+level->exchange_ghosts[justFaces].num_recvs)*sizeof(MPI_Request));
  level->exchange_ghosts[justFaces].status   = (MPI_Status *)malloc((level->exchange_ghosts[justFaces].num_sends+level->exchange_ghosts[justFaces].num_recvs)*sizeof(MPI_Status ));
  #endif


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //print_communicator(4,level->my_rank,0,&level->exchange_ghosts[justFaces]);
}



//---------------------------------------------------------------------------------------------------------------------------------------------------
void create_level(level_type *level, int boxes_in_i, int box_dim, int box_ghosts, int box_vectors, int domain_boundary_condition, int MPI_Rank, int MPI_Tasks){
  int box;
  int TotalBoxes = boxes_in_i*boxes_in_i*boxes_in_i;

  //if(MPI_Rank==0){printf("attempting to create a %5d^3 level using a %3d^3 grid of %3d^3 boxes with a target of %6.3f boxes per process...\n",box_dim*boxes_in_i,boxes_in_i,box_dim,(double)TotalBoxes/(double)MPI_Tasks);fflush(stdout);}
  if(MPI_Rank==0){printf("\nattempting to create a %5d^3 level using a %3d^3 grid of %3d^3 boxes...\n",box_dim*boxes_in_i,boxes_in_i,box_dim);fflush(stdout);}

  int omp_threads = 1;
  int omp_nested  = 0;
  #pragma omp parallel 
  {
    #pragma omp master
    {
      omp_threads = omp_get_num_threads();
      omp_nested  = omp_get_nested();
    }
  }

  level->memory_allocated = 0;
  level->box_dim        = box_dim;
  level->box_ghosts     = box_ghosts;
  level->box_vectors    = box_vectors;
  level->boxes_in.i     = boxes_in_i;
  level->boxes_in.j     = boxes_in_i;
  level->boxes_in.k     = boxes_in_i;
  level->dim.i          = box_dim*level->boxes_in.i;
  level->dim.j          = box_dim*level->boxes_in.j;
  level->dim.k          = box_dim*level->boxes_in.k;
  level->active         = 1;
  level->my_rank        = MPI_Rank;
  level->num_ranks      = MPI_Tasks;
  level->domain_boundary_condition = domain_boundary_condition;
  level->alpha_is_zero  = -1;
  level->num_threads      = omp_threads;
  // intra-box threading...
  level->threads_per_box  = omp_threads;
  level->concurrent_boxes = 1;
  // inter-box threading...
  //level->threads_per_box  = 1;
  //level->concurrent_boxes = omp_threads;

  // allocate 3D array of integers to hold the MPI rank of the corresponding box
  level->rank_of_box = (int*)malloc(level->boxes_in.i*level->boxes_in.j*level->boxes_in.k*sizeof(int));
  level->memory_allocated +=       (level->boxes_in.i*level->boxes_in.j*level->boxes_in.k*sizeof(int));
  for(box=0;box<level->boxes_in.i*level->boxes_in.j*level->boxes_in.k;box++){level->rank_of_box[box]=-1;}  // -1 denotes that there is no actual box assigned to this region

  // parallelize the grid...
  decompose_level_kd_tree(level->rank_of_box,level->boxes_in.i,level->boxes_in.i*level->boxes_in.j,0,0,0,level->boxes_in.i,level->boxes_in.j,level->boxes_in.k,0,MPI_Tasks);
  //print_decomposition(level);// for debug purposes only


  // build my list of boxes...
  level->num_my_boxes=0;
  for(box=0;box<level->boxes_in.i*level->boxes_in.j*level->boxes_in.k;box++){if(level->rank_of_box[box]==level->my_rank)level->num_my_boxes++;} 
  level->my_boxes = (box_type*)malloc(level->num_my_boxes*sizeof(box_type));
  box=0;
  int i,j,k;
  for(k=0;k<level->boxes_in.k;k++){
  for(j=0;j<level->boxes_in.j;j++){
  for(i=0;i<level->boxes_in.i;i++){
    int jStride = level->boxes_in.i;
    int kStride = level->boxes_in.i*level->boxes_in.j;
    int b=i + j*jStride + k*kStride;
    if(level->rank_of_box[b]==level->my_rank){
      level->memory_allocated += create_box(&level->my_boxes[box],level->box_vectors,level->box_dim,level->box_ghosts);
      level->my_boxes[box].low.i = i*level->box_dim;
      level->my_boxes[box].low.j = j*level->box_dim;
      level->my_boxes[box].low.k = k*level->box_dim;
      level->my_boxes[box].global_box_id = b;
      box++;
  }}}}
  initialize_valid_region(level); // define which cells are within the domain


  // Tune the OpenMP style of parallelism...
  if(omp_nested){
  #ifndef OMP_STENCILS_PER_THREAD
  #define OMP_STENCILS_PER_THREAD 64
  #endif
                                           level->concurrent_boxes = level->num_my_boxes;
  if(level->concurrent_boxes > omp_threads)level->concurrent_boxes = omp_threads;
  if(level->concurrent_boxes <           1)level->concurrent_boxes = 1;
  level->threads_per_box = omp_threads / level->concurrent_boxes;
  if(level->threads_per_box > level->box_dim*level->box_dim)
     level->threads_per_box = level->box_dim*level->box_dim; // JK collapse
  if(level->threads_per_box > level->box_dim*level->box_dim*level->box_dim/OMP_STENCILS_PER_THREAD )
     level->threads_per_box = level->box_dim*level->box_dim*level->box_dim/OMP_STENCILS_PER_THREAD;
  if(level->threads_per_box<1)level->threads_per_box = 1;
  }
  if(MPI_Rank==0){
    if(omp_nested)printf("  OMP_NESTED=TRUE  OMP_NUM_THREADS=%d ... %d teams of %d threads\n",omp_threads,level->concurrent_boxes,level->threads_per_box);
             else printf("  OMP_NESTED=FALSE OMP_NUM_THREADS=%d ... %d teams of %d threads\n",omp_threads,level->concurrent_boxes,level->threads_per_box);
    fflush(stdout);
  }


  // build an assist structure for Gauss Seidel Red Black that would facilitate unrolling and SIMDization...
  if(level->num_my_boxes){
    int i,j;
    int kStride = level->my_boxes[0].kStride;
    int jStride = level->my_boxes[0].jStride;

    posix_memalign((void**)&(level->RedBlack_FP[0]  ),64,kStride*sizeof(double  )); // even planes
    posix_memalign((void**)&(level->RedBlack_FP[1]  ),64,kStride*sizeof(double  ));
                      memset(level->RedBlack_FP[0]  ,  0,kStride*sizeof(double  ));
                      memset(level->RedBlack_FP[1]  ,  0,kStride*sizeof(double  ));
                              level->memory_allocated += kStride*sizeof(double  );
                              level->memory_allocated += kStride*sizeof(double  );
  
    for(j=0-level->box_ghosts;j<level->box_dim+level->box_ghosts;j++){
    for(i=0-level->box_ghosts;i<level->box_dim+level->box_ghosts;i++){
      int ij = (i+level->box_ghosts) + (j+level->box_ghosts)*jStride;
  //  if((i^j)&0x1)level->RedBlack_64bMask[ij]= ~0;else level->RedBlack_64bMask[ij]=  0; // useful for blend instructions
      if((i^j^1)&0x1)level->RedBlack_FP[  0][ij]=1.0;else level->RedBlack_FP[  0][ij]=0.0;
      if((i^j^1)&0x1)level->RedBlack_FP[  1][ij]=0.0;else level->RedBlack_FP[  1][ij]=1.0;
    }}
  }




  // create mini programs that affect ghost zone exchanges
  for(i=0;i<2;i++){
    level->exchange_ghosts[i].num_recvs    =0;
    level->exchange_ghosts[i].num_sends    =0;
    level->exchange_ghosts[i].num_blocks[0]=0;
    level->exchange_ghosts[i].num_blocks[1]=0;
    level->exchange_ghosts[i].num_blocks[2]=0;
  }
  build_exchange_ghosts(level,0); // faces, edges, corners
  build_exchange_ghosts(level,1); // justFaces


  // duplicate MPI_COMM_WORLD to be the communicator for each level
  #ifdef USE_MPI
  MPI_Comm_dup(MPI_COMM_WORLD,&level->MPI_COMM_LEVEL);
  #endif
    
  // report on potential load imbalance
  uint64_t BoxesPerProcess = level->num_my_boxes;
  #ifdef USE_MPI
  uint64_t send = level->num_my_boxes;
  MPI_Allreduce(&send,&BoxesPerProcess,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  #endif
  if(MPI_Rank==0){printf("  calculating boxes per process... target=%0.3f, max=%ld\n",(double)TotalBoxes/(double)MPI_Tasks,BoxesPerProcess);fflush(stdout);}
}



//---------------------------------------------------------------------------------------------------------------------------------------------------
void reset_level_timers(level_type *level){
  // cycle counters information...
  level->cycles.smooth                  = 0;
  level->cycles.apply_op                = 0;
  level->cycles.residual                = 0;
  level->cycles.blas1                   = 0;
  level->cycles.blas3                   = 0;
  level->cycles.boundary_conditions     = 0;
  level->cycles.restriction_total       = 0;
  level->cycles.restriction_pack        = 0;
  level->cycles.restriction_local       = 0;
  level->cycles.restriction_unpack      = 0;
  level->cycles.restriction_recv        = 0;
  level->cycles.restriction_send        = 0;
  level->cycles.restriction_wait        = 0;
  level->cycles.interpolation_total     = 0;
  level->cycles.interpolation_pack      = 0;
  level->cycles.interpolation_local     = 0;
  level->cycles.interpolation_unpack    = 0;
  level->cycles.interpolation_recv      = 0;
  level->cycles.interpolation_send      = 0;
  level->cycles.interpolation_wait      = 0;
  level->cycles.ghostZone_total         = 0;
  level->cycles.ghostZone_pack          = 0;
  level->cycles.ghostZone_local         = 0;
  level->cycles.ghostZone_unpack        = 0;
  level->cycles.ghostZone_recv          = 0;
  level->cycles.ghostZone_send          = 0;
  level->cycles.ghostZone_wait          = 0;
  level->cycles.collectives             = 0;
  level->cycles.Total                   = 0;
  // solver events information...
  level->Krylov_iterations              = 0;
  level->CAKrylov_formations_of_G       = 0;
  level->vcycles_from_this_level        = 0;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------
void destroy_level(level_type *level){
  int box;
  for(box=0;box>level->num_my_boxes;box++)destroy_box(&level->my_boxes[box]);
  free(level->rank_of_box);
}
