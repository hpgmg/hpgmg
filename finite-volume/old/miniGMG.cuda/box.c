//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>  
#include <stdlib.h>  
#include <stdint.h>  
#include <string.h>  
//------------------------------------------------------------------------------------------------------------------------------
#include "cuda_runtime.h"
#include "cuda.h"
//------------------------------------------------------------------------------------------------------------------------------
#include "defines.h"
#include "box.h"
//------------------------------------------------------------------------------------------------------------------------------
int create_box(box_type *box, int numGrids, int low_i, int low_j, int low_k, int dim_i, int dim_j, int dim_k, int ghosts){
  uint64_t memory_allocated = 0;
  box->numGrids = numGrids;
  box->low.i = low_i;
  box->low.j = low_j;
  box->low.k = low_k;
  box->dim.i = dim_i;
  box->dim.j = dim_j;
  box->dim.k = dim_k;
  box->dim_with_ghosts.i = dim_i+2*ghosts;
  box->dim_with_ghosts.j = dim_j+2*ghosts;
  box->dim_with_ghosts.k = dim_k+2*ghosts;
  box->ghosts = ghosts;
  box->pencil = (dim_i+2*ghosts);
  #if 1
  box->plane  =( ((dim_j+2*ghosts)*box->pencil)+0xF) & ~0xF; // pad plane to be a multiple of 128 bytes
  #else
  box->plane  = (dim_j+2*ghosts)*box->pencil;
  #endif
  box->volume = (dim_k+2*ghosts)*box->plane;


  // bufsizes represent the 26 neighboring boxes of the ghostzone
  //    faces are ghosts*dim*dim
  //    edges are ghosts*ghosts*dim
  // vertices are ghosts*ghosts*ghosts
  // buffer 13 (offset = 0,0,0) is the core of the box and is not communicated.  As such, its size is 0
  int di,dj,dk;
  for(dk=-1;dk<=1;dk++){
  for(dj=-1;dj<=1;dj++){
  for(di=-1;di<=1;di++){
    int n=13+di+3*dj+9*dk;
    box->bufsizes[n]=1;
    if(di==0)box->bufsizes[n]*=dim_i;else box->bufsizes[n]*=ghosts;
    if(dj==0)box->bufsizes[n]*=dim_j;else box->bufsizes[n]*=ghosts;
    if(dk==0)box->bufsizes[n]*=dim_k;else box->bufsizes[n]*=ghosts;
  }}}box->bufsizes[13]=0;

  // allocate buffers in one pass
  int n,total_bufsize = 0;
  for(n=0;n<27;n++)total_bufsize+=box->bufsizes[n];
  
  cudaMalloc( (void**) &(box->surface_bufs[0]), total_bufsize*sizeof(double)); cudaMemset( box->surface_bufs[0], 0, total_bufsize*sizeof(double));
  cudaMalloc( (void**) &(  box->ghost_bufs[0]), total_bufsize*sizeof(double)); cudaMemset(   box->ghost_bufs[0], 0, total_bufsize*sizeof(double));
  //cudaDeviceSynchronize();
  
  memory_allocated += total_bufsize*sizeof(double);
  memory_allocated += total_bufsize*sizeof(double);
  double *base;
  base=box->surface_bufs[0];for(n=0;n<27;n++){box->surface_bufs[n]=base;base+=box->bufsizes[n];}   // pointer to surfaces
  base=  box->ghost_bufs[0];for(n=0;n<27;n++){  box->ghost_bufs[n]=base;base+=box->bufsizes[n];}

  // allocate pointers to grids and grids themselves
  posix_memalign((void**)&(box->grids),64,box->numGrids*sizeof(double*));  
  memory_allocated += box->numGrids*sizeof(double*);
  int g;
  for(g=0;g<box->numGrids;g++){
    cudaMalloc( (void**) &(box->grids[g]), box->volume*sizeof(double)); 
    cudaMemset( box->grids[g], 0, box->volume*sizeof(double));
    memory_allocated += box->volume*sizeof(double);
  }

  return(memory_allocated);
}

void destroy_box(box_type *box){
  int g;for(g=0;g<box->numGrids;g++){
    cudaFree(box->grids[g]);
  }
  free(box->grids);
  cudaFree(box->ghost_bufs[0]);
  cudaFree(box->surface_bufs[0]);
}


