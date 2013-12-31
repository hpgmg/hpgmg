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
#include <unistd.h>
//------------------------------------------------------------------------------------------------------------------------------
#include "cuda_runtime.h"
#include "cuda.h"
//------------------------------------------------------------------------------------------------------------------------------
#include <omp.h>
#ifdef __MPI
#include <mpi.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include "defines.h"
#include "box.h"
#include "mg.h"
#include "kernels.h"
//------------------------------------------------------------------------------------------------------------------------------
uint64_t CycleTime(){
  uint64_t lo, hi;
  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
  return( (((uint64_t)hi) << 32) | ((uint64_t)lo) );
}

//------------------------------------------------------------------------------------------------------------------------------
#define enqueueEvent(id) domain->cudaEvents[domain->num_cudaEvents].level=level;domain->cudaEvents[domain->num_cudaEvents].type=id;cudaEventRecord(domain->cudaEvents[domain->num_cudaEvents].event,0);domain->num_cudaEvents++;
//------------------------------------------------------------------------------------------------------------------------------
int create_subdomain(subdomain_type * box, int subdomain_low_i, int subdomain_low_j, int subdomain_low_k,  
                                       int subdomain_dim_i, int subdomain_dim_j, int subdomain_dim_k, 
                                       int numGrids, int ghosts, int numLevels){
  int level;
  uint64_t memory_allocated=0;
  box->numLevels=numLevels;
  box->ghosts=ghosts;
  box->low.i = subdomain_low_i;
  box->low.j = subdomain_low_j;
  box->low.k = subdomain_low_k;
  box->dim.i = subdomain_dim_i;
  box->dim.j = subdomain_dim_j;
  box->dim.k = subdomain_dim_k;
  posix_memalign((void**)&(box->levels),64,numLevels*sizeof(box_type));
  memory_allocated += numLevels*sizeof(box_type);
  for(level=0;level<numLevels;level++){
    int numGridsActual = numGrids;
    #ifdef __USE_BICGSTAB
    if(level == (numLevels-1))numGridsActual+=6; // BiCGStab requires additional grids r0,r,p,s,Ap,As
    #endif
    memory_allocated += create_box(&box->levels[level],numGridsActual,subdomain_low_i>>level,subdomain_low_j>>level,subdomain_low_k>>level,
                                                                      subdomain_dim_i>>level,subdomain_dim_j>>level,subdomain_dim_k>>level,ghosts);
  }
  return(memory_allocated);
}

void destroy_subdomain(subdomain_type * box){
  int level;
  for(level=0;level<box->numLevels;level++){
    destroy_box(&box->levels[level]);
  }
  free(box->levels);
}



//------------------------------------------------------------------------------------------------------------------------------
int calculate_neighboring_subdomain_index(int bi, int bj, int bk, int di, int dj, int dk, int subdomains_per_rank_in_i, int subdomains_per_rank_in_j, int subdomains_per_rank_in_k){
  int ni,nj,nk;
  ni=(bi+subdomains_per_rank_in_i+di)%subdomains_per_rank_in_i;
  nj=(bj+subdomains_per_rank_in_j+dj)%subdomains_per_rank_in_j;
  nk=(bk+subdomains_per_rank_in_k+dk)%subdomains_per_rank_in_k;
  int index = ni + nj*subdomains_per_rank_in_i + nk*subdomains_per_rank_in_i*subdomains_per_rank_in_j;
  return(index);
}

int calculate_neighboring_subdomain_rank(int bi, int bj, int bk, int di, int dj, int dk, int subdomains_per_rank_in_i, int subdomains_per_rank_in_j, int subdomains_per_rank_in_k,
                                   int ri, int rj, int rk, int ranks_in_i, int ranks_in_j, int ranks_in_k){
  if((bi+di)<0)ri--;if((bi+di)>=subdomains_per_rank_in_i)ri++;ri=(ri+ranks_in_i)%ranks_in_i;
  if((bj+dj)<0)rj--;if((bj+dj)>=subdomains_per_rank_in_j)rj++;rj=(rj+ranks_in_j)%ranks_in_j;
  if((bk+dk)<0)rk--;if((bk+dk)>=subdomains_per_rank_in_k)rk++;rk=(rk+ranks_in_k)%ranks_in_k;
  int rank = ri + rj*ranks_in_i + rk*ranks_in_i*ranks_in_j;
  return(rank);
}

int on_same_face(int n1, int n2){
  int i1 = ((n1  )%3)-1;
  int j1 = ((n1/3)%3)-1;
  int k1 = ((n1/9)%3)-1;
  int i2 = ((n2  )%3)-1;
  int j2 = ((n2/3)%3)-1;
  int k2 = ((n2/9)%3)-1;
  if( (i1!=0)&&(i1==i2) )return(1);
  if( (j1!=0)&&(j1==j2) )return(1);
  if( (k1!=0)&&(k1==k2) )return(1);
  return(0);
}

int on_same_edge(int n1, int n2){
  int i1 = ((n1  )%3)-1;
  int j1 = ((n1/3)%3)-1;
  int k1 = ((n1/9)%3)-1;
  int i2 = ((n2  )%3)-1;
  int j2 = ((n2/3)%3)-1;
  int k2 = ((n2/9)%3)-1;
  if( (i1!=0) && (j1!=0) && (i1==i2) && (j1==j2) )return(1);
  if( (i1!=0) && (k1!=0) && (i1==i2) && (k1==k2) )return(1);
  if( (j1!=0) && (k1!=0) && (j1==j2) && (k1==k2) )return(1);
  return(0);
}
//------------------------------------------------------------------------------------------------------------------------------
int create_domain(domain_type * domain, 
              int subdomain_dim_i,  int subdomain_dim_j,  int subdomain_dim_k,  
              int subdomains_per_rank_in_i, int subdomains_per_rank_in_j, int subdomains_per_rank_in_k, 
              int ranks_in_i,      int ranks_in_j,      int ranks_in_k, 
              int rank,
              int numGrids, int ghosts, int numLevels
             ){
  int  i, j, k;
  int di,dj,dk;

  domain->rank = rank;
  if(domain->rank==0){printf("creating domain...       ");fflush(stdout);}
  if(ghosts>1){if(domain->rank==0)printf("#ghosts(%d)>1\n",ghosts);exit(0);}
  if(ghosts>subdomain_dim_i>>(numLevels-1)){if(domain->rank==0)printf("#ghosts(%d)>bottom grid size(%d)\n",ghosts,subdomain_dim_i>>(numLevels-1));exit(0);}

  if( (subdomain_dim_i!=subdomain_dim_j)||(subdomain_dim_j!=subdomain_dim_k)||(subdomain_dim_i!=subdomain_dim_k) ){
  if(domain->rank==0)printf("subdomain_dim's must be equal\n");exit(0);
  }
  uint64_t memory_allocated =0;
  // processes are laid out in x, then y, then z
  int my_k_rank = (rank                                                              ) / (ranks_in_i*ranks_in_j);
  int my_j_rank = (rank - (ranks_in_i*ranks_in_j*my_k_rank)                          ) / (ranks_in_i           );
  int my_i_rank = (rank - (ranks_in_i*ranks_in_j*my_k_rank) - (ranks_in_i*my_j_rank) );
  //printf("%2d: (%2d,%2d,%2d)\n",domain->rank,my_k_rank,my_j_rank,my_i_rank);

  for(dk=-1;dk<=1;dk++){
  for(dj=-1;dj<=1;dj++){
  for(di=-1;di<=1;di++){
    int n = 13+di+3*dj+9*dk;
    int neighbor_rank_in_i = (my_i_rank+di+ranks_in_i)%ranks_in_i;
    int neighbor_rank_in_j = (my_j_rank+dj+ranks_in_j)%ranks_in_j;
    int neighbor_rank_in_k = (my_k_rank+dk+ranks_in_k)%ranks_in_k;
    domain->rank_of_neighbor[n] = neighbor_rank_in_i + ranks_in_i*neighbor_rank_in_j + ranks_in_i*ranks_in_j*neighbor_rank_in_k;
  }}}

  domain->subdomains_per_rank_in.i = subdomains_per_rank_in_i;
  domain->subdomains_per_rank_in.j = subdomains_per_rank_in_j;
  domain->subdomains_per_rank_in.k = subdomains_per_rank_in_k;
  domain->numsubdomains = subdomains_per_rank_in_i*subdomains_per_rank_in_j*subdomains_per_rank_in_k;
  posix_memalign((void**)&(domain->subdomains),64,domain->numsubdomains*sizeof(subdomain_type));
  memory_allocated+=domain->numsubdomains*sizeof(subdomain_type);

  domain->subdomains_in.i = domain->subdomains_per_rank_in.i * ranks_in_i;
  domain->subdomains_in.j = domain->subdomains_per_rank_in.j * ranks_in_j;
  domain->subdomains_in.k = domain->subdomains_per_rank_in.k * ranks_in_k;
  domain->dim.i = domain->subdomains_in.i * subdomain_dim_i;
  domain->dim.j = domain->subdomains_in.j * subdomain_dim_j;
  domain->dim.k = domain->subdomains_in.k * subdomain_dim_k;
  domain->numLevels = numLevels;
  domain->numGrids  = numGrids;
  domain->ghosts = ghosts;

  // subdomains within a process are laid out in i, then j, then k
  for(k=0;k<subdomains_per_rank_in_k;k++){
  for(j=0;j<subdomains_per_rank_in_j;j++){
  for(i=0;i<subdomains_per_rank_in_i;i++){
    int box = i + j*subdomains_per_rank_in_i + k*subdomains_per_rank_in_i*subdomains_per_rank_in_j;
    int low_i = subdomain_dim_i * (i + subdomains_per_rank_in_i*my_i_rank);
    int low_j = subdomain_dim_j * (j + subdomains_per_rank_in_j*my_j_rank);
    int low_k = subdomain_dim_k * (k + subdomains_per_rank_in_k*my_k_rank);
    memory_allocated += create_subdomain(&domain->subdomains[box],low_i,low_j,low_k,
                                         //(i+subdomains_per_rank_in_i)*subdomain_dim_i,(j+subdomains_per_rank_in_j)*subdomain_dim_j,(k+subdomains_per_rank_in_k)*subdomain_dim_k,
                                            subdomain_dim_i,                             subdomain_dim_j,                             subdomain_dim_k,
                                            numGrids,ghosts,numLevels);
    for(dk=-1;dk<=1;dk++){
    for(dj=-1;dj<=1;dj++){
    for(di=-1;di<=1;di++){
      int n = 13+di+3*dj+9*dk;
    domain->subdomains[box].neighbors[n].rank        = calculate_neighboring_subdomain_rank( i,j,k,di,dj,dk, subdomains_per_rank_in_i,subdomains_per_rank_in_j,subdomains_per_rank_in_k,
                                                                                                               my_i_rank,my_j_rank,my_k_rank,ranks_in_i,ranks_in_j,ranks_in_k);
    domain->subdomains[box].neighbors[n].local_index = calculate_neighboring_subdomain_index(i,j,k,di,dj,dk, subdomains_per_rank_in_i,subdomains_per_rank_in_j,subdomains_per_rank_in_k);
      //printf("rank=%2d, box[%2d].neighbors[%3d]=(%3d,%3d)\n",domain->rank,box,n,domain->subdomains[box].neighbors[n].rank,domain->subdomains[box].neighbors[n].local_index);
    }}}
  }}}


  #ifdef __MPI
  int   FaceSizeAtLevel0 = subdomain_dim_i*subdomain_dim_i*ghosts;
  int   EdgeSizeAtLevel0 = subdomain_dim_i*ghosts*ghosts;
  int CornerSizeAtLevel0 = ghosts*ghosts*ghosts;

  if(ghosts==1){
        EdgeSizeAtLevel0 = 0;
      CornerSizeAtLevel0 = 0;
  }

  int   FaceBufSizePerSubdomain = FaceSizeAtLevel0 + 4*EdgeSizeAtLevel0 + 4*CornerSizeAtLevel0;
  int   EdgeBufSizePerSubdomain =                      EdgeSizeAtLevel0 + 2*CornerSizeAtLevel0;
  int CornerBufSizePerSubdomain =                                           CornerSizeAtLevel0;

  int    faces[27] = {0,0,0,0,1,0,0,0,0,  0,1,0,1,0,1,0,1,0,  0,0,0,0,1,0,0,0,0};
  int    edges[27] = {0,1,0,1,0,1,0,1,0,  1,0,1,0,0,0,1,0,1,  0,1,0,1,0,1,0,1,0};
  int  corners[27] = {1,0,1,0,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,  1,0,1,0,0,0,1,0,1};

  #define __USE_CUDAMALLOCHOST
  #ifndef __USE_CUDAMALLOCHOST
  #warning You should #define __USE_CUDAMALLOCHOST to accelerate cudaMemcpy over PCIe
  #endif
  // allocate MPI send/recv buffers for the 26 neighbors....
  for(dk=-1;dk<=1;dk++){
  for(dj=-1;dj<=1;dj++){
  for(di=-1;di<=1;di++){
    int n = 13+di+3*dj+9*dk;
    if(faces[n]){
      int SubdomainsWritingToThisBuffer;
      if(di!=0)SubdomainsWritingToThisBuffer=domain->subdomains_per_rank_in.j*domain->subdomains_per_rank_in.k;
      if(dj!=0)SubdomainsWritingToThisBuffer=domain->subdomains_per_rank_in.i*domain->subdomains_per_rank_in.k;
      if(dk!=0)SubdomainsWritingToThisBuffer=domain->subdomains_per_rank_in.i*domain->subdomains_per_rank_in.j;
      #ifdef __USE_CUDAMALLOCHOST
      cudaMallocHost((void**)&(    domain->send_buffer[n]),   FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      cudaMallocHost((void**)&(    domain->recv_buffer[n]),   FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      #else
      posix_memalign((void**)&(    domain->send_buffer[n]),64,FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      posix_memalign((void**)&(    domain->recv_buffer[n]),64,FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      #endif
          cudaMalloc((void**)&(domain->gpu_send_buffer[n]),   FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
          cudaMalloc((void**)&(domain->gpu_recv_buffer[n]),   FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
                            memset(domain->send_buffer[n],0,  FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
                            memset(domain->recv_buffer[n],0,  FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
                                          memory_allocated+=  FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double);
                                          memory_allocated+=  FaceBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double);
      domain->buffer_size[n].faces   = 1*SubdomainsWritingToThisBuffer;
      domain->buffer_size[n].edges   = 4*SubdomainsWritingToThisBuffer;
      domain->buffer_size[n].corners = 4*SubdomainsWritingToThisBuffer;
    }
    if(edges[n]){
      int SubdomainsWritingToThisBuffer;
      if(di==0)SubdomainsWritingToThisBuffer=domain->subdomains_per_rank_in.i;
      if(dj==0)SubdomainsWritingToThisBuffer=domain->subdomains_per_rank_in.j;
      if(dk==0)SubdomainsWritingToThisBuffer=domain->subdomains_per_rank_in.k;
      #ifdef __USE_CUDAMALLOCHOST
      cudaMallocHost((void**)&(    domain->send_buffer[n]),   EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      cudaMallocHost((void**)&(    domain->recv_buffer[n]),   EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      #else
      posix_memalign((void**)&(    domain->send_buffer[n]),64,EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      posix_memalign((void**)&(    domain->recv_buffer[n]),64,EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      #endif
          cudaMalloc((void**)&(domain->gpu_send_buffer[n]),   EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
          cudaMalloc((void**)&(domain->gpu_recv_buffer[n]),   EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
                            memset(domain->send_buffer[n],0,  EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
                            memset(domain->recv_buffer[n],0,  EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
                                          memory_allocated+=  EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double);
                                          memory_allocated+=  EdgeBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double);
      domain->buffer_size[n].faces   = 0;
      domain->buffer_size[n].edges   = 1*SubdomainsWritingToThisBuffer;
      domain->buffer_size[n].corners = 2*SubdomainsWritingToThisBuffer;
    }
    if(corners[n]){
      int SubdomainsWritingToThisBuffer=1;
      #ifdef __USE_CUDAMALLOCHOST
      cudaMallocHost((void**)&(    domain->send_buffer[n]),   CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      cudaMallocHost((void**)&(    domain->recv_buffer[n]),   CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      #else
      posix_memalign((void**)&(    domain->send_buffer[n]),64,CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      posix_memalign((void**)&(    domain->recv_buffer[n]),64,CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
      #endif
          cudaMalloc((void**)&(domain->gpu_send_buffer[n]),   CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
          cudaMalloc((void**)&(domain->gpu_recv_buffer[n]),   CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
                            memset(domain->send_buffer[n],0,  CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
                            memset(domain->recv_buffer[n],0,  CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double));
                                          memory_allocated+=  CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double);
                                          memory_allocated+=  CornerBufSizePerSubdomain*SubdomainsWritingToThisBuffer*sizeof(double);
      domain->buffer_size[n].faces   = 0;
      domain->buffer_size[n].edges   = 0;
      domain->buffer_size[n].corners = 1;
    }
  }}}
  cudaMalloc((void**)&(domain->gpu_pointers_to_gpu_send_buffer),27*sizeof(double*));
  cudaMalloc((void**)&(domain->gpu_pointers_to_gpu_recv_buffer),27*sizeof(double*));
  cudaMemcpy( domain->gpu_pointers_to_gpu_send_buffer, domain->gpu_send_buffer, 27*sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy( domain->gpu_pointers_to_gpu_recv_buffer, domain->gpu_recv_buffer, 27*sizeof(double*), cudaMemcpyHostToDevice);

  //if(domain->rank==0)printf("\n");
  struct{int faces,edges,corners;}send_offset[27][28]; // offset[buffer][neighbor]... ie location of neighbor
  struct{int faces,edges,corners;}recv_offset[27][28]; // offset[buffer][neighbor]... ie location of neighbor

  int n1,n2,n3;
  for(n1=0;n1<27;n1++){
      send_offset[n1][27].faces  =0;
      send_offset[n1][27].edges  =0;
      send_offset[n1][27].corners=0;
      recv_offset[n1][27].faces  =0;
      recv_offset[n1][27].edges  =0;
      recv_offset[n1][27].corners=0;
    for(n2=0;n2<27;n2++){
      n3 = 26-n2;
      send_offset[n1][n2].faces  =0;
      send_offset[n1][n2].edges  =0;
      send_offset[n1][n2].corners=0;
      recv_offset[n1][n3].faces  =0;
      recv_offset[n1][n3].edges  =0;
      recv_offset[n1][n3].corners=0;
      if( (faces[n1] && on_same_face(n1,n2)) || (edges[n1] && on_same_edge(n1,n2)) || (corners[n1] && (n1==n2)) ){
        //if(domain->rank==0)printf("%d (%d,%d) (%d,%d) (%d,%d)\n",n2,faces[n1],on_same_face(n1,n2),edges[n1],on_same_edge(n1,n2),corners[n1],(n1==n2));
        send_offset[n1][n2].faces  =send_offset[n1][27].faces;  send_offset[n1][27].faces  +=  faces[n2];
        send_offset[n1][n2].edges  =send_offset[n1][27].edges;  send_offset[n1][27].edges  +=  edges[n2];
        send_offset[n1][n2].corners=send_offset[n1][27].corners;send_offset[n1][27].corners+=corners[n2];
      }
      if( (faces[n1] && on_same_face(n1,n3)) || (edges[n1] && on_same_edge(n1,n3)) || (corners[n1] && (n1==n3)) ){
        //if(domain->rank==0)printf("%d (%d,%d) (%d,%d) (%d,%d)\n",n3,faces[n1],on_same_face(n1,n3),edges[n1],on_same_edge(n1,n3),corners[n1],(n1==n3));
        recv_offset[n1][n3].faces  =recv_offset[n1][27].faces;  recv_offset[n1][27].faces  +=  faces[n3];
        recv_offset[n1][n3].edges  =recv_offset[n1][27].edges;  recv_offset[n1][27].edges  +=  edges[n3];
        recv_offset[n1][n3].corners=recv_offset[n1][27].corners;recv_offset[n1][27].corners+=corners[n3];
      }
    }
    //if(domain->rank==0)printf("%2d: ",n1);
    //for(n2=0;n2<27;n2++){
    //  if(domain->rank==0)printf("(%d,%d,%d) ",recv_offset[n1][n2].faces,recv_offset[n1][n2].edges,recv_offset[n1][n2].corners);
    //}
    //if(domain->rank==0)printf("\n");
  }

  // for all ijk
  // for all dijk
  for(k=0;k<subdomains_per_rank_in_k;k++){
  for(j=0;j<subdomains_per_rank_in_j;j++){
  for(i=0;i<subdomains_per_rank_in_i;i++){
    int box = i + j*subdomains_per_rank_in_i + k*subdomains_per_rank_in_i*subdomains_per_rank_in_j;
    for(dk=-1;dk<=1;dk++){
    for(dj=-1;dj<=1;dj++){
    for(di=-1;di<=1;di++){
      int n = 13+di+3*dj+9*dk;
         domain->subdomains[box].neighbors[n].send.buf = -1;
         domain->subdomains[box].neighbors[n].recv.buf = -1;
      if(domain->subdomains[box].neighbors[n].rank != domain->rank){
        int ni = i+di;
        int nj = j+dj;
        int nk = k+dk;
        int buf_di=0;if(ni<0)buf_di=-1;if(ni>=subdomains_per_rank_in_i)buf_di=1;
        int buf_dj=0;if(nj<0)buf_dj=-1;if(nj>=subdomains_per_rank_in_j)buf_dj=1;
        int buf_dk=0;if(nk<0)buf_dk=-1;if(nk>=subdomains_per_rank_in_k)buf_dk=1;
        // now decide where this box writes surface data to and where it reads ghost data from
        int buf;
        int base=-1000;
        buf = 13+buf_di+3*buf_dj+9*buf_dk;
        if(faces[buf]){ 
          // send
          if(buf_di!=0){base = (j+k*subdomains_per_rank_in_j);}// jk plane
          if(buf_dj!=0){base = (i+k*subdomains_per_rank_in_i);}// ik plane
          if(buf_dk!=0){base = (i+j*subdomains_per_rank_in_i);}// ij plane
          domain->subdomains[box].neighbors[n].send.buf = buf;
          domain->subdomains[box].neighbors[n].send.offset.faces   = 1*base+send_offset[buf][n].faces;
          domain->subdomains[box].neighbors[n].send.offset.edges   = 4*base+send_offset[buf][n].edges;
          domain->subdomains[box].neighbors[n].send.offset.corners = 4*base+send_offset[buf][n].corners;
          // recv
          if(buf_di!=0){base = (nj+nk*subdomains_per_rank_in_j);}// jk plane
          if(buf_dj!=0){base = (ni+nk*subdomains_per_rank_in_i);}// ik plane
          if(buf_dk!=0){base = (ni+nj*subdomains_per_rank_in_i);}// ij plane
          domain->subdomains[box].neighbors[n].recv.buf = buf;
          domain->subdomains[box].neighbors[n].recv.offset.faces   = 1*base+recv_offset[buf][n].faces;
          domain->subdomains[box].neighbors[n].recv.offset.edges   = 4*base+recv_offset[buf][n].edges;
          domain->subdomains[box].neighbors[n].recv.offset.corners = 4*base+recv_offset[buf][n].corners;
        }
        if(edges[buf]){ 
          // send
          if(buf_di==0){base = i;} // edge along i
          if(buf_dj==0){base = j;} // edge along j
          if(buf_dk==0){base = k;} // edge along k
          domain->subdomains[box].neighbors[n].send.buf = buf;
          domain->subdomains[box].neighbors[n].send.offset.faces   = 0*base+send_offset[buf][n].faces;
          domain->subdomains[box].neighbors[n].send.offset.edges   = 1*base+send_offset[buf][n].edges;
          domain->subdomains[box].neighbors[n].send.offset.corners = 2*base+send_offset[buf][n].corners;
          // recv
          if(buf_di==0){base = ni;} // edge along i
          if(buf_dj==0){base = nj;} // edge along j
          if(buf_dk==0){base = nk;} // edge along k
          domain->subdomains[box].neighbors[n].recv.buf = buf;
          domain->subdomains[box].neighbors[n].recv.offset.faces   = 0*base+recv_offset[buf][n].faces;
          domain->subdomains[box].neighbors[n].recv.offset.edges   = 1*base+recv_offset[buf][n].edges;
          domain->subdomains[box].neighbors[n].recv.offset.corners = 2*base+recv_offset[buf][n].corners;
        }
        if(corners[buf]){
          // send
          domain->subdomains[box].neighbors[n].send.buf = buf;
          domain->subdomains[box].neighbors[n].send.offset.faces   = send_offset[buf][n].faces;
          domain->subdomains[box].neighbors[n].send.offset.edges   = send_offset[buf][n].edges;
          domain->subdomains[box].neighbors[n].send.offset.corners = send_offset[buf][n].corners;
          // recv
          domain->subdomains[box].neighbors[n].recv.buf = buf;
          domain->subdomains[box].neighbors[n].recv.offset.faces   = recv_offset[buf][n].faces;
          domain->subdomains[box].neighbors[n].recv.offset.edges   = recv_offset[buf][n].edges;
          domain->subdomains[box].neighbors[n].recv.offset.corners = recv_offset[buf][n].corners;
        }
  
        //if(domain->rank==0)printf("box=%2d, n=%2d: send.buf=%2d, send.offset={%2d,%2d,%2d}\n",box,n,
        //  domain->subdomains[box].neighbors[n].send.buf,
        //  domain->subdomains[box].neighbors[n].send.offset.faces,   
        //  domain->subdomains[box].neighbors[n].send.offset.edges,  
        //  domain->subdomains[box].neighbors[n].send.offset.corners);
        //if(domain->rank==0)printf("box=%2d, n=%2d: recv.buf=%2d, recv.offset={%2d,%2d,%2d} ijk=(%2d,%2d,%2d), nijk=(%2d,%2d,%2d)\n",box,n,
        //  domain->subdomains[box].neighbors[n].recv.buf,
        //  domain->subdomains[box].neighbors[n].recv.offset.faces,   
        //  domain->subdomains[box].neighbors[n].recv.offset.edges,  
        //  domain->subdomains[box].neighbors[n].recv.offset.corners,i,j,k,ni,nj,nk);

      }
    }}}
  }}}
  #endif

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // create a list of events for timing purposes...
  int ev;
  domain->num_cudaEvents = 0; 
  domain->cudaEvents = (cudaMGEvent_type*)malloc(max_cudaEvents*sizeof(cudaMGEvent_type));
  for(ev=0;ev<max_cudaEvents;ev++){
    domain->cudaEvents[ev].type  = -1;
    domain->cudaEvents[ev].level = -1;
    cudaEventCreate(&domain->cudaEvents[ev].event);
  }
  cudaEventCreate(&domain->cycles.communicationStart);
  cudaEventCreate(&domain->cycles.communicationEnd);


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // create a buffer on the device for reductions (norm, dot, etc...)
  cudaMalloc((void**)&domain->gpu_reduction_buffer,100*sizeof(double));


  int box,level,grid,n;
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // create domain->gpu_subdomains (gpu shadow copy of all domain meta data...)
  subdomain_type * temp_subdomains = (subdomain_type *)malloc(domain->numsubdomains*sizeof(subdomain_type));
                    memcpy(temp_subdomains,domain->subdomains,domain->numsubdomains*sizeof(subdomain_type));
  for(box=0;box<domain->numsubdomains;box++){
    box_type *temp_levels =        (box_type *)malloc(domain->numLevels*sizeof(box_type));
    memcpy(temp_levels,domain->subdomains[box].levels,domain->numLevels*sizeof(box_type));
    //temp_subdomains[box].levels = temp_levels;
    for(level=0;level<domain->numLevels;level++){
      //printf("%d\n",domain->subdomains[box].levels[level].numGrids);
      cudaMalloc((void**)&temp_levels[level].grids,                                            domain->subdomains[box].levels[level].numGrids*sizeof(double*));
      cudaMemcpy(         temp_levels[level].grids,domain->subdomains[box].levels[level].grids,domain->subdomains[box].levels[level].numGrids*sizeof(double*),cudaMemcpyHostToDevice);
    }
    cudaMalloc((void**)&temp_subdomains[box].levels,            domain->numLevels*sizeof(box_type));
    cudaMemcpy(         temp_subdomains[box].levels,temp_levels,domain->numLevels*sizeof(box_type),cudaMemcpyHostToDevice);
    free(temp_levels);
  }
  cudaMalloc((void**)&domain->gpu_subdomains,                domain->numsubdomains*sizeof(subdomain_type));
  cudaMemcpy(         domain->gpu_subdomains,temp_subdomains,domain->numsubdomains*sizeof(subdomain_type),cudaMemcpyHostToDevice);
  free(temp_subdomains);


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // flatten surface/ghost_buf pointers for easy indexing...
  /*
  double ** cpu_pointers_to_gpu_surface_bufs = (double**)malloc(domain->numsubdomains*27*domain->numLevels*sizeof(double*));
  double ** cpu_pointers_to_gpu_ghost_bufs   = (double**)malloc(domain->numsubdomains*27*domain->numLevels*sizeof(double*));
  for(box=0;box<domain->numsubdomains;box++){
    for(level=0;level<domain->numLevels;level++){
      for(n=0;n<27;n++){
        cpu_pointers_to_gpu_surface_bufs[box*domain->numLevels*27 + level*27 + n] = (domain->subdomains[box].levels[level].surface_bufs[n]);
        cpu_pointers_to_gpu_ghost_bufs[  box*domain->numLevels*27 + level*27 + n] = (domain->subdomains[box].levels[level].ghost_bufs[n]  );
  }}}
  cudaMalloc((void**)&domain->gpu_pointers_to_gpu_surface_bufs, domain->numsubdomains*27*domain->numLevels*sizeof(double*));
  cudaMalloc((void**)&domain->gpu_pointers_to_gpu_ghost_bufs,   domain->numsubdomains*27*domain->numLevels*sizeof(double*));
  cudaMemcpy( domain->gpu_pointers_to_gpu_surface_bufs, cpu_pointers_to_gpu_surface_bufs, domain->numsubdomains*27*domain->numLevels*sizeof(double*), cudaMemcpyHostToDevice);
  cudaMemcpy( domain->gpu_pointers_to_gpu_ghost_bufs,   cpu_pointers_to_gpu_ghost_bufs,   domain->numsubdomains*27*domain->numLevels*sizeof(double*), cudaMemcpyHostToDevice);
  */

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // flatten neighbor data for easy indexing...
  /*
  neighbor_type * cpu_pointer_to_cpu_neighbor_data = (neighbor_type*)malloc(domain->numsubdomains*27*sizeof(neighbor_type));
  for(box=0;box<domain->numsubdomains;box++){
    for(n=0;n<27;n++){
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].rank                = domain->subdomains[box].neighbors[n].rank;
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].local_index         = domain->subdomains[box].neighbors[n].local_index;
      #ifdef __MPI
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].send.buf            = domain->subdomains[box].neighbors[n].send.buf;
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].send.offset.faces   = domain->subdomains[box].neighbors[n].send.offset.faces;
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].send.offset.edges   = domain->subdomains[box].neighbors[n].send.offset.edges;
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].send.offset.corners = domain->subdomains[box].neighbors[n].send.offset.corners;
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].recv.buf            = domain->subdomains[box].neighbors[n].recv.buf;
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].recv.offset.faces   = domain->subdomains[box].neighbors[n].recv.offset.faces;
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].recv.offset.edges   = domain->subdomains[box].neighbors[n].recv.offset.edges;
      cpu_pointer_to_cpu_neighbor_data[box*27 + n].recv.offset.corners = domain->subdomains[box].neighbors[n].recv.offset.corners;
      #endif
  }}
  cudaMalloc((void**)&domain->gpu_pointer_to_gpu_neighbor_data, domain->numsubdomains*27*sizeof(neighbor_type));
  cudaMemcpy( domain->gpu_pointer_to_gpu_neighbor_data, cpu_pointer_to_cpu_neighbor_data, domain->numsubdomains*27*sizeof(neighbor_type), cudaMemcpyHostToDevice);
  */


  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(domain->rank==0){
  printf("done\n");fflush(stdout); 
  printf("  %d x %d x %d (per subdomain)\n",subdomain_dim_i,subdomain_dim_j,subdomain_dim_k);
  printf("  %d x %d x %d (per process)\n",subdomains_per_rank_in_i*subdomain_dim_i,subdomains_per_rank_in_j*subdomain_dim_j,subdomains_per_rank_in_k*subdomain_dim_k);
  printf("  %d x %d x %d (overall)\n",ranks_in_i*subdomains_per_rank_in_i*subdomain_dim_i,ranks_in_j*subdomains_per_rank_in_j*subdomain_dim_j,ranks_in_k*subdomains_per_rank_in_k*subdomain_dim_k);
  printf("  %d-deep ghost zones\n",ghosts);
  printf("  allocated %llu MB\n",memory_allocated>>20);
  fflush(stdout);
  }
  return(memory_allocated);
}

void destroy_domain(domain_type * domain){
  if(domain->rank==0){printf("deallocating domain...   ");fflush(stdout);}
  int box;for(box=0;box<domain->numsubdomains;box++){
    destroy_subdomain(&domain->subdomains[box]);
  }
  free(domain->subdomains);
  cudaFree(domain->gpu_reduction_buffer);
  if(domain->rank==0){printf("done\n");fflush(stdout);}
}

//------------------------------------------------------------------------------------------------------------------------------
void print_timing(domain_type *domain){
  if(domain->rank!=0)return;
  int level,numLevels = domain->numLevels;
  uint64_t _timeStart=CycleTime();sleep(1);uint64_t _timeEnd=CycleTime();
  float SecondsPerCycle = (float)1.0/(float)(_timeEnd-_timeStart);

  cudaDeviceSynchronize();
  for(level=0;level<(numLevels  );level++){
    domain->cycles.Total[level]        =0;
    domain->cycles.smooth[level]       =0;
    domain->cycles.apply_op[level]     =0;
    domain->cycles.residual[level]     =0;
    domain->cycles.restriction[level]  =0;
    domain->cycles.interpolation[level]=0;
    domain->cycles.blas1[level]        =0;
    domain->cycles.blas3[level]        =0;
  //domain->cycles.communication[level]=0;
    domain->cycles.pack[level]         =0;
    domain->cycles.bufcopy[level]      =0;
    domain->cycles.unpack[level]       =0;
    domain->cycles.s2buf[level]        =0;
    domain->cycles.buf2g[level]        =0;
  }
  int ev;
  for(ev=0;ev<domain->num_cudaEvents;ev+=2){
    float deltaTime;cudaEventElapsedTime(&deltaTime,domain->cudaEvents[ev].event,domain->cudaEvents[ev+1].event); // assumes events are queued back-to-back
    if(domain->cudaEvents[ev].type == cudaEvent_smooth       )domain->cycles.smooth[       domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_apply_op     )domain->cycles.apply_op[     domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_residual     )domain->cycles.residual[     domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_restriction  )domain->cycles.restriction[  domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_interpolation)domain->cycles.interpolation[domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_blas1        )domain->cycles.blas1[        domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_blas3        )domain->cycles.blas3[        domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_pack         )domain->cycles.pack[         domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_bufcopy      )domain->cycles.bufcopy[      domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_unpack       )domain->cycles.unpack[       domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_s2buf        )domain->cycles.s2buf[        domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_buf2g        )domain->cycles.buf2g[        domain->cudaEvents[ev].level]+=deltaTime;
    if(domain->cudaEvents[ev].type == cudaEvent_PCIe         )domain->cycles.PCIe[         domain->cudaEvents[ev].level]+=deltaTime;
//  if(domain->cudaEvents[ev].type == cudaEvent_communication)domain->cycles.communication[domain->cudaEvents[ev].level]+=deltaTime;
  }

  float total,seconds;
          printf("                       ");for(level=0;level<(numLevels  );level++){printf("%12d ",level);}printf("\n");
          printf("                       ");for(level=0;level<(numLevels  );level++){printf("%10d^3 ",domain->subdomains[0].levels[0].dim.i>>level);}printf("       total\n");
  total=0;printf("smooth                 ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(       domain->cycles.smooth[level]);printf("%12.6f ",seconds);domain->cycles.Total[level]+=seconds;total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("apply_op               ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(     domain->cycles.apply_op[level]);printf("%12.6f ",seconds);domain->cycles.Total[level]+=seconds;total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("residual               ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(     domain->cycles.residual[level]);printf("%12.6f ",seconds);domain->cycles.Total[level]+=seconds;total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("restriction            ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(  domain->cycles.restriction[level]);printf("%12.6f ",seconds);domain->cycles.Total[level]+=seconds;total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("interpolation          ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(domain->cycles.interpolation[level]);printf("%12.6f ",seconds);domain->cycles.Total[level]+=seconds;total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("BLAS1                  ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(        domain->cycles.blas1[level]);printf("%12.6f ",seconds);domain->cycles.Total[level]+=seconds;total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("communication          ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(domain->cycles.communication[level])+
                                                                                             SecondsPerCycle*(float)(  domain->cycles.collectives[level]);printf("%12.6f ",seconds);domain->cycles.Total[level]+=seconds;total+=seconds;}printf("%12.6f\n",total);
//total=0;printf("communication          ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(domain->cycles.communication[level]);printf("%12.6f ",seconds);domain->cycles.Total[level]+=seconds;total+=seconds;}printf("%12.6f\n",total);
//total=0;printf("communication          ");for(level=0;level<(numLevels  );level++){seconds=SecondsPerCycle*(float)(domain->cycles.communication[level]);printf("%12.6f ",seconds);domain->cycles.Total[level]+=seconds;total+=seconds;}printf("%12.6f\n",total);
//#ifdef __PRINT_COMMUNICATION_BREAKDOWN
  total=0;printf("  surface to buffers   ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(        domain->cycles.s2buf[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("  exchange buffers     ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(      domain->cycles.bufcopy[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("  buffers to ghosts    ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(        domain->cycles.buf2g[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("  PCIe transfers       ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(         domain->cycles.PCIe[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  #ifdef __MPI
  total=0;printf("  pack MPI buffers     ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(         domain->cycles.pack[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("  unpack MPI buffers   ");for(level=0;level<(numLevels  );level++){seconds=          0.001*(float)(       domain->cycles.unpack[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("  MPI_Isend            ");for(level=0;level<(numLevels  );level++){seconds=SecondsPerCycle*(float)(         domain->cycles.send[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("  MPI_Irecv            ");for(level=0;level<(numLevels  );level++){seconds=SecondsPerCycle*(float)(         domain->cycles.recv[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("  MPI_Waitall          ");for(level=0;level<(numLevels  );level++){seconds=SecondsPerCycle*(float)(         domain->cycles.wait[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  total=0;printf("  MPI_collectives      ");for(level=0;level<(numLevels  );level++){seconds=SecondsPerCycle*(float)(  domain->cycles.collectives[level]);printf("%12.6f ",seconds);                                     total+=seconds;}printf("%12.6f\n",total);
  #endif
//#endif
  total=0;printf("--------------         ");for(level=0;level<(numLevels+1);level++){printf("------------ ");}printf("\n");
  total=0;printf("Total by level         ");for(level=0;level<(numLevels  );level++){printf("%12.6f " ,domain->cycles.Total[level]);total+=domain->cycles.Total[level];}
                                                                                     printf("%12.6f\n",total);
  printf("\n");
  printf( "Total time in MGSolve   %12.6f\n",SecondsPerCycle*(float)domain->cycles.MGSolve);
  printf("            \" MGBuild   %12.6f\n",SecondsPerCycle*(float)domain->cycles.build);
  printf("            \" v-cycles  %12.6f\n",SecondsPerCycle*(float)domain->cycles.vcycles);
  printf( "    number of v-cycles  %12d\n"  ,domain->vcycles_performed);
  printf( "   BiCGStab iterations  %12d\n"  ,domain->BiCGStab_iterations);
  printf("\n\n");fflush(stdout);
}


//------------------------------------------------------------------------------------------------------------------------------
void exchange_boundary(domain_type *domain, int level, int grid_id, int exchange_faces, int exchange_edges, int exchange_corners){
//uint64_t _timeStartCommunication = CycleTime();
  cudaEventRecord(domain->cycles.communicationStart,0);

  uint64_t _timeStart,_timeEnd;
  int n;
  if(domain->ghosts>1){printf("ERROR, CUDA version of exchange_boundary is hard-coded to only operate on faces\n");exit(0);}
  int    faces[27] = {0,0,0,0,1,0,0,0,0,  0,1,0,1,0,1,0,1,0,  0,0,0,0,1,0,0,0,0};
  int    edges[27] = {0,1,0,1,0,1,0,1,0,  1,0,1,0,0,0,1,0,1,  0,1,0,1,0,1,0,1,0};
  int  corners[27] = {1,0,1,0,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,  1,0,1,0,0,0,1,0,1};
  int exchange[27] = {0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0};
  #warning Sam, change the code so that faces/edges/corners are mandates and not suggestions...
  for(n=0;n<27;n++){
    if(                       exchange_faces   )exchange[n] |=   faces[n];
    if( (domain->ghosts>1) && exchange_edges   )exchange[n] |=   edges[n];
    if( (domain->ghosts>1) && exchange_corners )exchange[n] |= corners[n];
  }

  #ifdef __MPI
  int   FaceSizeAtLevel = domain->subdomains[0].levels[level].dim.i*domain->subdomains[0].levels[level].dim.i*domain->ghosts;
  int   EdgeSizeAtLevel = domain->subdomains[0].levels[level].dim.i*domain->ghosts*domain->ghosts;
  int CornerSizeAtLevel = domain->ghosts*domain->ghosts*domain->ghosts;

  if(domain->ghosts==1){
      EdgeSizeAtLevel=0;
    CornerSizeAtLevel=0;
  }

  MPI_Request requests[54];
  MPI_Status    status[54];
  int nMessages=0;
  #endif

  // loop through bufs, prepost Irecv's
  #ifdef __MPI
  _timeStart = CycleTime();
  for(n=0;n<27;n++){if(exchange[26-n] && (domain->rank_of_neighbor[26-n] != domain->rank) ){
    int size = FaceSizeAtLevel*domain->buffer_size[26-n].faces + EdgeSizeAtLevel*domain->buffer_size[26-n].edges + CornerSizeAtLevel*domain->buffer_size[26-n].corners;
    #ifdef __GPU_DIRECT
    MPI_Irecv(domain->gpu_recv_buffer[26-n],size,MPI_DOUBLE,domain->rank_of_neighbor[26-n],n,MPI_COMM_WORLD,&requests[nMessages]);
    #else
    MPI_Irecv(    domain->recv_buffer[26-n],size,MPI_DOUBLE,domain->rank_of_neighbor[26-n],n,MPI_COMM_WORLD,&requests[nMessages]);
    #endif
    nMessages++;
  }}
  _timeEnd = CycleTime();
  domain->cycles.recv[level] += (uint64_t)(_timeEnd-_timeStart);
  #endif

  // extract surface, pack into surface_bufs
  grid_to_surface_buffers(domain,level,grid_id);


  // pack domain buffers
  #ifdef __MPI
  //GPU copies faces into gpu_send_buffers[] with offset FaceSizeAtLevel*domain->subdomains[sendBox].neighbors[n].send.offset.faces
  surface_buffers_to_send_buffer(domain,level,grid_id);
  #ifndef __GPU_DIRECT   // without GPU_DIRECT, you have to DMA the MPI send buffers to the host over PCIe
  enqueueEvent(cudaEvent_PCIe);
  for(n=0;n<27;n++){if(exchange[n] && (domain->rank_of_neighbor[n] != domain->rank) ){
    int size = FaceSizeAtLevel*domain->buffer_size[n].faces + EdgeSizeAtLevel*domain->buffer_size[n].edges + CornerSizeAtLevel*domain->buffer_size[n].corners;
    cudaMemcpyAsync(domain->send_buffer[n],domain->gpu_send_buffer[n],size*sizeof(double),cudaMemcpyDeviceToHost,0);
  }}
  enqueueEvent(cudaEvent_PCIe);
  cudaDeviceSynchronize(); // wait for Memcpy's to finish before initiating Isend's...
  #endif
  #endif



  // loop through bufs, post Isend's
  #ifdef __MPI
  _timeStart = CycleTime();
  for(n=0;n<27;n++){if(exchange[n] && (domain->rank_of_neighbor[n] != domain->rank) ){
    int size = FaceSizeAtLevel*domain->buffer_size[n].faces + EdgeSizeAtLevel*domain->buffer_size[n].edges + CornerSizeAtLevel*domain->buffer_size[n].corners;
    #ifdef __GPU_DIRECT
    MPI_Isend(domain->gpu_send_buffer[n],size,MPI_DOUBLE,domain->rank_of_neighbor[n],n,MPI_COMM_WORLD,&requests[nMessages]);
    #else
    MPI_Isend(    domain->send_buffer[n],size,MPI_DOUBLE,domain->rank_of_neighbor[n],n,MPI_COMM_WORLD,&requests[nMessages]);
    #endif
    nMessages++;
  }}
  _timeEnd = CycleTime();
  domain->cycles.send[level] += (uint64_t)(_timeEnd-_timeStart);
  #endif


  // exchange locally... try and hide within Isend latency... 
  surface_buffers_to_ghost_buffers(domain,level,grid_id);


  #ifdef __MPI
  // loop through bufs, MPI_Wait on recvs
  _timeStart = CycleTime();
  MPI_Waitall(nMessages,requests,status);
  _timeEnd = CycleTime();
  domain->cycles.wait[level] += (uint64_t)(_timeEnd-_timeStart);
  #endif


  // unpack domain buffers 
  #ifdef __MPI
  #ifndef __GPU_DIRECT   // without GPU_DIRECT, you have to DMA the MPI receive buffers to the device over PCIe
  enqueueEvent(cudaEvent_PCIe);
  for(n=0;n<27;n++){if(exchange[n] && (domain->rank_of_neighbor[n] != domain->rank) ){
    int size = FaceSizeAtLevel*domain->buffer_size[n].faces + EdgeSizeAtLevel*domain->buffer_size[n].edges + CornerSizeAtLevel*domain->buffer_size[n].corners;
    cudaMemcpyAsync(domain->gpu_recv_buffer[n],domain->recv_buffer[n],size*sizeof(double),cudaMemcpyHostToDevice,0);
  }}
  enqueueEvent(cudaEvent_PCIe);
  #endif
  //cudaDeviceSynchronize(); // FIX - necessary ???
  //GPU extracts faces from gpu_recv_buffers[] with offset FaceSizeAtLevel*domain->subdomains[recvBox].neighbors[n].recv.offset.faces
  recv_buffer_to_ghost_buffers(domain,level,grid_id);
  #endif


  // unpack ghost_bufs, insert into grid
  ghost_buffers_to_grid(domain,level,grid_id);


  cudaEventRecord(domain->cycles.communicationEnd,0);
  cudaEventSynchronize(domain->cycles.communicationEnd);
  float deltaTime;cudaEventElapsedTime(&deltaTime,domain->cycles.communicationStart,domain->cycles.communicationEnd);
  domain->cycles.communication[level] += deltaTime;

//cudaDeviceSynchronize(); // mainly for timing purposes...
//uint64_t _timeEndCommunication = CycleTime();
//domain->cycles.communication[level] += (uint64_t)(_timeEndCommunication-_timeStartCommunication);
}


//------------------------------------------------------------------------------------------------------------------------------
double norm(domain_type * domain, int level, int grid_id){
  double cpu_norm = 0.0;
  enqueueEvent(cudaEvent_PCIe);
  cudaMemcpy(domain->gpu_reduction_buffer,&cpu_norm, 1*sizeof(double), cudaMemcpyHostToDevice);
  enqueueEvent(cudaEvent_PCIe);
  norm_on_gpu(domain,level,grid_id,domain->gpu_reduction_buffer);
  enqueueEvent(cudaEvent_PCIe);
  cudaMemcpy(&cpu_norm,domain->gpu_reduction_buffer, 1*sizeof(double), cudaMemcpyDeviceToHost);
  enqueueEvent(cudaEvent_PCIe);
  #ifdef __MPI
  uint64_t _timeStartAllReduce = CycleTime();
  double send = cpu_norm;
  MPI_Allreduce(&send,&cpu_norm,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  uint64_t _timeEndAllReduce = CycleTime();
  domain->cycles.collectives[level]   += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce);
//domain->cycles.communication[level] += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce); // unlike the CPU version, 'communication' doesn't include 'collectives'
  #endif
  return(cpu_norm);
}

//------------------------------------------------------------------------------------------------------------------------------
double dot(domain_type * domain, int level, int id_a, int id_b){
  double cpu_dot = 0.0;
  enqueueEvent(cudaEvent_PCIe);
  cudaMemcpy(domain->gpu_reduction_buffer,&cpu_dot, 1*sizeof(double), cudaMemcpyHostToDevice);
  enqueueEvent(cudaEvent_PCIe);
  dot_on_gpu(domain,level,id_a,id_b,domain->gpu_reduction_buffer);
  enqueueEvent(cudaEvent_PCIe);
  cudaMemcpy(&cpu_dot,domain->gpu_reduction_buffer, 1*sizeof(double), cudaMemcpyDeviceToHost);
  enqueueEvent(cudaEvent_PCIe);
  #ifdef __MPI
  uint64_t _timeStartAllReduce = CycleTime();
  double send = cpu_dot;
  MPI_Allreduce(&send,&cpu_dot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  uint64_t _timeEndAllReduce = CycleTime();
  domain->cycles.collectives[level]   += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce);
//domain->cycles.communication[level] += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce); // unlike the CPU version, 'communication' doesn't include 'collectives'
  #endif
  return(cpu_dot);
}

//------------------------------------------------------------------------------------------------------------------------------
double mean(domain_type * domain, int level, int id_a){
  double sum_domain = 0.0;
  enqueueEvent(cudaEvent_PCIe);
  cudaMemcpy(domain->gpu_reduction_buffer,&sum_domain, 1*sizeof(double), cudaMemcpyHostToDevice);
  enqueueEvent(cudaEvent_PCIe);
  sum_on_gpu(domain,level,id_a,domain->gpu_reduction_buffer);
  enqueueEvent(cudaEvent_PCIe);
  cudaMemcpy(&sum_domain,domain->gpu_reduction_buffer, 1*sizeof(double), cudaMemcpyDeviceToHost);
  enqueueEvent(cudaEvent_PCIe);

  double ncells_domain = (double)domain->dim.i*(double)domain->dim.j*(double)domain->dim.k;

  #ifdef __MPI
  uint64_t _timeStartAllReduce = CycleTime();
  double send = sum_domain;
  MPI_Allreduce(&send,&sum_domain,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  uint64_t _timeEndAllReduce = CycleTime();
  domain->cycles.collectives[level]   += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce);
//domain->cycles.communication[level] += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce); // unlike the CPU version, 'communication' doesn't include 'collectives'
  #endif

  double mean_domain = sum_domain / ncells_domain;
  if(domain->rank==0){printf("mean_domain = %20.12e  ncells_domain=%20.12e\n",mean_domain,ncells_domain);fflush(stdout);}
  return(mean_domain);
}


//------------------------------------------------------------------------------------------------------------------------------
void MGSolve(domain_type * domain, int u_id, int F_id, double a, double b, double h0, double desired_mg_norm){ 
  int level;
  double hLevel;
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  domain->num_cudaEvents=0;
  for(level=0;level<10;level++){
  domain->cycles.communication[level] = 0;
  domain->cycles.recv[level]          = 0;
  domain->cycles.send[level]          = 0;
  domain->cycles.wait[level]          = 0;
  domain->cycles.collectives[level]   = 0;
  }
//domain->cycles.build                = 0;
  domain->cycles.vcycles              = 0;
  domain->cycles.MGSolve              = 0;
  domain->vcycles_performed           = 0;
  domain->BiCGStab_iterations         = 0;
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ConfigureGPU();
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int e_id = __u;
  int R_id = __f_minus_Av;
  int box,v,s;
  int ghosts = domain->ghosts;
  int numLevels = domain->numLevels;
  int numSmooths       = 4; 
  int numSmoothsBottom = 48;
  #ifdef __TEST_MG_CONVERGENCE
  int maxVCycles       = 20;
  #else
  int maxVCycles       = 10;
  #endif
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(domain->rank==0){printf("MGSolve...");fflush(stdout);}
  uint64_t _timeStartMGSolve = CycleTime();
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // make initial guess for e (=0) and setup the RHS
  #if 0
   zero_grid(domain,0,e_id);                // ee = 0
  scale_grid(domain,0,R_id,1.0,F_id);       // R_id = F_id
  #else
   mul_grids(domain,0,e_id,1.0,F_id,__lambda);  // e_id = F_id*lambda
  scale_grid(domain,0,R_id,1.0,F_id);           // R_id = F_id
  #endif
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  uint64_t _timeStartVCycle = CycleTime();
  for(v=0;v<maxVCycles;v++){
    domain->vcycles_performed++;
    // down the v-cycle.........................................................................................................
    for(level=0;level<(numLevels-1);level++){
      hLevel = h0 * (double)(1<<level);
      if(ghosts>1)exchange_boundary(domain,level,R_id,1,1,1); // need to get RHS from neighbors if duplicating their work (need ghosts-1)
      for(s=0;s<numSmooths;s+=ghosts){
                    exchange_boundary(domain,level,e_id,1,1,1);
                               smooth(domain,level,e_id,R_id,a,b,hLevel,s);
      } // relax
                    exchange_boundary(domain,level,e_id,1,0,0); // technically only needs to be a 1-deep ghost zone & faces only
                             residual(domain,level,__temp,e_id,R_id,a,b,hLevel);
                          restriction(domain,level,R_id,__temp);
                            zero_grid(domain,level+1,e_id);
    } // down
  
    // bottom solve... (exchange phi,rhs) ......................................................................................
    level = numLevels-1;
    hLevel = h0 * (double)(1<<level);
    #ifdef __USE_BICGSTAB
      #warning Using BiCGStab bottom solver...
      // based on scanned page sent to me by Erin/Nick
      // Algorithm 7.7 in Iterative Methods for Sparse Linear Systems(Yousef Saad)???
      int DiagonallyPrecondition=1;
      double desired_reduction_in_norm = 1e-3;
      int jMax=200;
      int j=0;
      int BiCGStabFailed    = 0;
      int BiCGStabConverged = 0;
      exchange_boundary(domain,level,e_id,1,0,0);                       // exchange_boundary(e_id)
      residual(domain,level,__r0,e_id,R_id,a,b,hLevel);                 // r0[] = R_id[] - A(e_id)
      scale_grid(domain,level,__r,1.0,__r0);                            // r[] = r0[]
      scale_grid(domain,level,__p,1.0,__r0);                            // p[] = r0[]
      double r_dot_r0 = dot(domain,level,__r,__r0);                     // r_dot_r0 = dot(r,r0)
      double norm_of_r0 = norm(domain,level,__r);                       // the norm of the initial residual... if L2 norm, then isn't this sqrt(r_dot_r0) ?
      if(norm_of_r0 == 0.0){BiCGStabConverged=1;}                       // entered BiCGStab with exact solutio
      while( (j<jMax) && (!BiCGStabFailed) && (!BiCGStabConverged) ){   // while(not done){
        j++;domain->BiCGStab_iterations++;
        if(DiagonallyPrecondition){
        mul_grids(domain,level,__temp,1.0,__lambda,__p);                //   temp[] = lambda[]*p[]
        exchange_boundary(domain,level,__temp,1,0,0);                   //   exchange_boundary(p)
        apply_op(domain,level,__Ap,__temp,a,b,hLevel);                  //   Ap = AD^{-1}(p)
        }else{
        exchange_boundary(domain,level,__p,1,0,0);                      //   exchange_boundary(p)
        apply_op(domain,level,__Ap,__p,a,b,hLevel);                     //   Ap = A(p)
        }
        //apply_op(domain,level,__Ap,__p,a,b,hLevel);                   //   Ap = A(p)
        double Ap_dot_r0 = dot(domain,level,__Ap,__r0);                 //   Ap_dot_r0 = dot(Ap,r0)
        if(Ap_dot_r0 == 0.0){BiCGStabFailed=1;break;}                   //   pivot breakdown ???
        double alpha = r_dot_r0 / Ap_dot_r0;                            //   alpha = r_dot_r0 / Ap_dot_r0
        add_grids(domain,level,e_id,1.0,e_id, alpha,__p );              //   e_id[] = e_id[] + alpha*p[]
        add_grids(domain,level,__s ,1.0,__r ,-alpha,__Ap);              //   s[]    = r[]    - alpha*Ap[]   (intermediate residual?)
        double norm_of_s = norm(domain,level,__s);                      //   FIX - redundant??  norm of intermediate residual
        if(norm_of_s == 0.0){BiCGStabConverged=1;break;}                //   FIX - redundant??  if As_dot_As==0, then As must be 0 which implies s==0
        if(norm_of_s < desired_reduction_in_norm*norm_of_r0){BiCGStabConverged=1;break;}
        if(DiagonallyPrecondition){
        mul_grids(domain,level,__temp,1.0,__lambda,__s);                //   temp[] = lambda[]*s[]
        exchange_boundary(domain,level,__temp,1,0,0);                   //   exchange_boundary(s)
        apply_op(domain,level,__As,__temp,a,b,hLevel);                  //   As = AD^{-1}(s)
        }else{
        exchange_boundary(domain,level,__s,1,0,0);                      //   exchange_boundary(s)
        apply_op(domain,level,__As,__s,a,b,hLevel);                     //   As = A(s)
        }
      //apply_op(domain,level,__As,__s,a,b,hLevel);                     //   As = A(s)
        double As_dot_As = dot(domain,level,__As,__As);                 //   As_dot_As = dot(As,As)
        double As_dot_s  = dot(domain,level,__As, __s);                 //   As_dot_s  = dot(As, s)
        if(As_dot_As == 0.0){BiCGStabConverged=1;break;}                //   converged ?
        double omega = As_dot_s / As_dot_As;                            //   omega = As_dot_s / As_dot_As
        if(omega == 0.0){BiCGStabFailed=1;break;}                       //   stabilization breakdown ???
        add_grids(domain,level,  e_id,  1.0,e_id, omega,__s   );        //   e_id[] = e_id[] + omega*s[]
        add_grids(domain,level,__r   ,  1.0,__s ,-omega,__As  );        //   r[]    = s[]    - omega*As[]  (recursively computed / updated residual)
        double norm_of_r = norm(domain,level,__r);                      //   norm of recursively computed residual (good enough??)
        if(norm_of_r == 0.0){BiCGStabConverged=1;break;}
        if(norm_of_r < desired_reduction_in_norm*norm_of_r0){BiCGStabConverged=1;break;}
        double r_dot_r0_new = dot(domain,level,__r,__r0);               //   r_dot_r0_new = dot(r,r0)
        if(r_dot_r0_new == 0.0){BiCGStabFailed=1;break;}                //   Lanczos breakdown ???
        double beta = (r_dot_r0_new/r_dot_r0) * (alpha/omega);          //   beta = (r_dot_r0_new/r_dot_r0) * (alpha/omega)
        add_grids(domain,level,__temp,1.0,__p,-omega,__Ap  );           //   __temp =         (p[]-omega*Ap[])
        add_grids(domain,level,__p   ,1.0,__r,  beta,__temp);           //   p[] = r[] + beta*(p[]-omega*Ap[])
        r_dot_r0 = r_dot_r0_new;                                        //   r_dot_r0 = r_dot_r0_new   (save old r_dot_r0)
      }                                                                 // }
      if(DiagonallyPrecondition){
      mul_grids(domain,level,e_id,1.0,__lambda,e_id);                   //   e_id[] = lambda[]*e_id[] // i.e. x = D^{-1}x'
      }
    #else // just multiple GSRB's
      #warning Defaulting to simple GSRB bottom solver with fixed number of iterations...
      // relax...
      if(ghosts>1)exchange_boundary(domain,level,R_id,1,1,1);
      for(s=0;s<numSmoothsBottom;s+=ghosts){
                    exchange_boundary(domain,level,e_id,1,1,1);
                               smooth(domain,level,e_id,R_id,a,b,hLevel,s);
      } // relax
    #endif
  
    // back up the v-cycle......................................................................................................
    for(level=(numLevels-2);level>=0;level--){
      hLevel = h0 * (double)(1<<level);
      interpolation(domain,level,e_id,e_id);
      for(s=0;s<numSmooths;s+=ghosts){
        exchange_boundary(domain,level,e_id,1,1,1); // always communicate
                   smooth(domain,level,e_id,R_id,a,b,hLevel,s);
      }
    } // up

    #if defined(__PRINT_NORM) || defined(__TEST_MG_CONVERGENCE)
    // check norm of residual...................................................................................................
    double norm_of_residual = 0.0;
    exchange_boundary(domain,0,e_id,1,0,0);
    residual(domain,0,__temp,e_id,R_id,a,b,h0);
    #if 1
    #warning Using ||D^{-1}(b-Ax)||_{inf} as convergence criteria...
    mul_grids(domain,0,__temp,1.0,__temp,__lambda); /// <<< precondition by D^{-1} ???
    #endif
    norm_of_residual = norm(domain,0,__temp);
    #if defined(__PRINT_NORM)
    if(domain->rank==0){
      if(v==0)printf("\n");
              printf("v-cycle=%2d, norm=%22.20f (%12.6e)\n",v+1,norm_of_residual,norm_of_residual);
      fflush(stdout);
    }
    #endif
    #if defined(__TEST_MG_CONVERGENCE)
    if(norm_of_residual<desired_mg_norm)break;
    #endif
    #endif
  } // maxVCycles
  cudaDeviceSynchronize(); // mainly for timing purposes...
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  domain->cycles.vcycles += (uint64_t)(CycleTime()-_timeStartVCycle);
  domain->cycles.MGSolve += (uint64_t)(CycleTime()-_timeStartMGSolve);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(domain->rank==0){printf("done\n");fflush(stdout);}
}

//------------------------------------------------------------------------------------------------------------------------------
