//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdint.h>
#include "defines.h"
#include "box.h"
#include "mg.h"

//------------------------------------------------------------------------------------------------------------------------------
#define enqueueEvent(id) domain->cudaEvents[domain->num_cudaEvents].level=level;domain->cudaEvents[domain->num_cudaEvents].type=id;cudaEventRecord(domain->cudaEvents[domain->num_cudaEvents].event,0);domain->num_cudaEvents++;
//------------------------------------------------------------------------------------------------------------------------------
__constant__ int    faces[27] = {0,0,0,0,1,0,0,0,0,  0,1,0,1,0,1,0,1,0,  0,0,0,0,1,0,0,0,0};
//------------------------------------------------------------------------------------------------------------------------------


#warning  !!! Remember to try out __ldg()
#warning  !!! Everything is hard-coded for ghosts=1

#ifdef VL
#warning using vectorized smooth()
__global__ void __smooth_once_GSRB(subdomain_type * gpu_subdomains, int phi_id, int rhs_id, double a, double b, double h, int sweep, int level){

  double h2inv = 1.0/(h*h);
  int box = blockIdx.z;  // CUDA 4 !!!

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  // |--pencil--|--IJStride--|--pencil--|--??--|
  int IJStride = (VL - pencil - pencil) & ~0x0F;
  int ij = ((IJStride*blockIdx.x + pencil - pencil) & ~0xF) + threadIdx.x; // i.e. shift vector so that thread0 is 128-byte aligned
//int i = ij % pencil; // pencil is even and I will AND off all but the LSB
  int i = ij & 0x1; // pencil is even and I will AND off all but the LSB
  int j = ij / pencil;
  int k; 
  
  #ifdef __POINTERS_IN_SHARED
  __shared__ double * phi;
  __shared__ double * rhs;
  __shared__ double * alpha;
  __shared__ double * beta_i;
  __shared__ double * beta_j;
  __shared__ double * beta_k;
  __shared__ double * lambda;
  if(threadIdx.x==0){
       phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane);
       rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane);
     alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane);
    beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane);
    beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane);
    beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane);
    lambda = gpu_subdomains[box].levels[level].grids[__lambda] + (plane);
  }
  __syncthreads();
  #else
  double *         phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane);
  double *         rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane);
  double *       alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane);
  double *      beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane);
  double *      beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane);
  double *      beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane);
  double *      lambda = gpu_subdomains[box].levels[level].grids[__lambda] + (plane);
  #endif


  __shared__ double beta_i_ijk[VL];
  __shared__ double beta_j_ijk[VL];
  __shared__ double       temp[VL];
  double beta_k_ijk,beta_k_ijkp1;
  double    phi_ijk,   phi_ijkp1, phi_ijkm1;

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                                   int withinBounds = 1;
  if( threadIdx.x <  pencil           )withinBounds = 0;
  if( threadIdx.x >= pencil+IJStride  )withinBounds = 0;
  if( ij < pencil+1                   )withinBounds = 0;
  if( ij >=pencil*(subdomain_dim+1)-1 )withinBounds = 0;


  k = 0;
  int ijk = ij + k*plane;

     phi_ijkm1 =    phi[ijk-plane];
     phi_ijk   =    phi[ijk      ];
     phi_ijkp1 =    phi[ijk+plane];
  beta_k_ijk   = beta_k[ijk      ];
  beta_k_ijkp1 = beta_k[ijk+plane];

  int RedBlackUpdate = (i^j^sweep)&0x1;


  for(k=0;k<subdomain_dim;k++){

    __syncthreads(); // WAR guard
                                                                      temp[threadIdx.x] = phi_ijk;
    if( (threadIdx.x>=pencil)&&(threadIdx.x<pencil+IJStride+1) )beta_i_ijk[threadIdx.x] = beta_i[ijk];
    if( (threadIdx.x>=pencil)                                  )beta_j_ijk[threadIdx.x] = beta_j[ijk];
    __syncthreads(); // RAW guard

    if(withinBounds){
      double helmholtz_ijk = a*alpha[ijk]*phi_ijk - b*h2inv*(
        beta_i_ijk[threadIdx.x+     1] * ( temp[threadIdx.x+     1] - phi_ijk                  ) -
        beta_i_ijk[threadIdx.x       ] * ( phi_ijk                  - temp[threadIdx.x-     1] ) +
        beta_j_ijk[threadIdx.x+pencil] * ( temp[threadIdx.x+pencil] - phi_ijk                  ) -
        beta_j_ijk[threadIdx.x       ] * ( phi_ijk                  - temp[threadIdx.x-pencil] ) +
        beta_k_ijkp1                   * ( phi_ijkp1                - phi_ijk                  ) -
        beta_k_ijk                     * ( phi_ijk                  - phi_ijkm1                )
      );
      // GSRB
      double new_phi = phi_ijk - lambda[ijk]*(helmholtz_ijk - rhs[ijk]);
      phi[ijk] = (RedBlackUpdate) ? new_phi : phi_ijk;
    }

    RedBlackUpdate=RedBlackUpdate^1;
    ijk+=plane;
    // rotate register pipeline...
       phi_ijkm1 =    phi_ijk;
       phi_ijk   =    phi_ijkp1;
       phi_ijkp1 =    phi[ijk+plane];
    beta_k_ijk   = beta_k_ijkp1;
    beta_k_ijkp1 = beta_k[ijk+plane];

  } // for k
} // GSRB kernel

#else
#ifdef __LOCALITY_VIA_SHARED
//==============================================================================================================================================================
// shared+L1 version (relies on shared memory but you must still favorL1)
//==============================================================================================================================================================
__global__ void __smooth_once_GSRB(subdomain_type * gpu_subdomains, int phi_id, int rhs_id, double a, double b, double h, int sweep, int level){

  double h2inv = 1.0/(h*h);
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------------------
  #ifdef __POINTERS_IN_SHARED
  __shared__ double * phi;
  __shared__ double * rhs;
  __shared__ double * alpha;
  __shared__ double * beta_i;
  __shared__ double * beta_j;
  __shared__ double * beta_k;
  __shared__ double * lambda;
  if((threadIdx.x==0)&&(threadIdx.y==0)){
       phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane+pencil+1);
       rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane+pencil+1);
     alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
    beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
    beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
    beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
    lambda = gpu_subdomains[box].levels[level].grids[__lambda] + (plane+pencil+1);
  }
  __syncthreads();
  #else
  double *         phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane+pencil+1);
  double *         rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane+pencil+1);
  double *       alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
  double *      beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
  double *      beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
  double *      beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
  double *      lambda = gpu_subdomains[box].levels[level].grids[__lambda] + (plane+pencil+1);
  #endif


  __shared__ double beta_i_ijk[TBDIMY  ][TBDIMX+1];
  __shared__ double beta_j_ijk[TBDIMY+1][TBDIMX  ];
  __shared__ double       temp[TBDIMY+2][TBDIMX+2]; // always index with threadIdx.y+1,x+1
  double beta_k_ijk,beta_k_ijkp1;
  double    phi_ijk,   phi_ijkp1, phi_ijkm1;
  

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

     k = 0; 
     int ijk = k*plane + j*pencil + i;

     phi_ijkm1 = phi[ijk-plane];
     phi_ijk   = phi[ijk      ];
     phi_ijkp1 = phi[ijk+plane];
  beta_k_ijk   = beta_k[ijk      ];
  beta_k_ijkp1 = beta_k[ijk+plane];

  int RedBlackUpdate = (i^j^sweep) & 0x1;

  for(k=0;k<subdomain_dim;k++){

    __syncthreads(); // WAR guard
    // copy phi and beta_i/j to shared temp array...
                            temp[threadIdx.y+1][threadIdx.x+1] = phi_ijk;
    if(threadIdx.x==0)      temp[threadIdx.y+1][            0] = phi[ijk-1];
    if(threadIdx.y==0)      temp[            0][threadIdx.x+1] = phi[ijk-pencil];
    if(threadIdx.x==0)      temp[threadIdx.y+1][     TBDIMX+1] = phi[ijk+TBDIMX];
    if(threadIdx.y==0)      temp[     TBDIMY+1][threadIdx.x+1] = phi[ijk+TBDIMY*pencil];
                      beta_i_ijk[threadIdx.y  ][threadIdx.x  ] = beta_i[ijk];
    if(threadIdx.x==0)beta_i_ijk[threadIdx.y  ][     TBDIMX  ] = beta_i[ijk+TBDIMX];
                      beta_j_ijk[threadIdx.y  ][threadIdx.x  ] = beta_j[ijk];
    if(threadIdx.y==0)beta_j_ijk[     TBDIMY  ][threadIdx.x  ] = beta_j[ijk+TBDIMY*pencil];
    __syncthreads(); // RAW guard

    if(withinBounds){
      double helmholtz_ijk = a*alpha[ijk]*phi_ijk - b*h2inv*( 
        beta_i_ijk[threadIdx.y  ][threadIdx.x+1] * ( temp[threadIdx.y+1][threadIdx.x+2] - phi_ijk                            ) - 
        beta_i_ijk[threadIdx.y  ][threadIdx.x  ] * ( phi_ijk                            - temp[threadIdx.y+1][threadIdx.x  ] ) +
        beta_j_ijk[threadIdx.y+1][threadIdx.x  ] * ( temp[threadIdx.y+2][threadIdx.x+1] - phi_ijk                            ) -
        beta_j_ijk[threadIdx.y  ][threadIdx.x  ] * ( phi_ijk                            - temp[threadIdx.y  ][threadIdx.x+1] ) +
        beta_k_ijkp1                             * ( phi_ijkp1                          - phi_ijk                            ) -
        beta_k_ijk                               * ( phi_ijk                            - phi_ijkm1                          )
      );
      // GSRB
      double new_phi = phi_ijk - lambda[ijk]*(helmholtz_ijk - rhs[ijk]);
      phi[ijk] = (RedBlackUpdate) ? new_phi : phi_ijk;
    }
 
    //if(withinBounds && ((RedBlackUpdate^k)&0x1)){ 
    //  phi[ijk] = phi_ijk - lambda[ijk]*(helmholtz_ijk - rhs[ijk]); // GSRB
    //}
    //double new_phi = phi_ijk - lambda[ijk]*(helmholtz_ijk - rhs[ijk]);
    //if(withinBounds){ 
    //  phi[ijk] = ((RedBlackUpdate^k)&0x1) ? new_phi : phi_ijk;
    //}

    RedBlackUpdate=RedBlackUpdate^1;
    ijk+=plane;
    // rotate register pipeline...
    phi_ijkm1    = phi_ijk;
    phi_ijk      = phi_ijkp1;
    phi_ijkp1    = phi[ijk+plane];
    beta_k_ijk   = beta_k_ijkp1;
    beta_k_ijkp1 = beta_k[ijk+plane];

  } // for k
} // GSRB kernel
#else
//==============================================================================================================================================================
// cache version (relies solely on the L1/L2 cache hierarchy)
//==============================================================================================================================================================
__global__ void __smooth_once_GSRB(subdomain_type * gpu_subdomains, int phi_id, int rhs_id, double a, double b, double h, int sweep, int level){

  double h2inv = 1.0/(h*h);
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------------------
  #ifdef __POINTERS_IN_SHARED
  __shared__ double * phi;
  __shared__ double * rhs;
  __shared__ double * alpha;
  __shared__ double * beta_i;
  __shared__ double * beta_j;
  __shared__ double * beta_k;
  __shared__ double * lambda;
  if((threadIdx.x==0)&&(threadIdx.y==0)){
       phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane+pencil+1);
       rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane+pencil+1);
     alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
    beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
    beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
    beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
    lambda = gpu_subdomains[box].levels[level].grids[__lambda] + (plane+pencil+1);
  }
  __syncthreads();
  #else
  double *         phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane+pencil+1);
  double *         rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane+pencil+1);
  double *       alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
  double *      beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
  double *      beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
  double *      beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
  double *      lambda = gpu_subdomains[box].levels[level].grids[__lambda] + (plane+pencil+1);
  #endif


  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  int RedBlackUpdate = (i^j^sweep) & 0x1;

  for(k=0;k<subdomain_dim;k++){
    int ijk = k*plane + j*pencil + i;

    if(withinBounds){
      double phi_ijk = phi[ijk];
      double helmholtz_ijk = a*alpha[ijk]*phi[ijk] - b*h2inv*( 
        beta_i[ijk     +1] * ( phi[ijk     +1] - phi[ijk       ] ) -
        beta_i[ijk       ] * ( phi[ijk       ] - phi[ijk     -1] ) + 
        beta_j[ijk+pencil] * ( phi[ijk+pencil] - phi[ijk       ] ) -
        beta_j[ijk       ] * ( phi[ijk       ] - phi[ijk-pencil] ) + 
        beta_k[ijk +plane] * ( phi[ijk +plane] - phi[ijk       ] ) -
        beta_k[ijk       ] * ( phi[ijk       ] - phi[ijk -plane] )
      );
      // GSRB
      double new_phi = phi_ijk - lambda[ijk]*(helmholtz_ijk - rhs[ijk]);
      phi[ijk] = (RedBlackUpdate) ? new_phi : phi_ijk;
    }
    //if(withinBounds && ((RedBlackUpdate^k)&0x1)){
    //  phi[ijk] = phi_ijk - lambda[ijk]*(helmholtz_ijk - rhs[ijk]); // GSRB
    //}
    //double new_phi = phi_ijk - lambda[ijk]*(helmholtz_ijk - rhs[ijk]);
    //if(withinBounds){ 
    //  phi[ijk] = ((RedBlackUpdate^k)&0x1) ? new_phi : phi_ijk;
    //}
    RedBlackUpdate=RedBlackUpdate^1;
  } // for k
} // GSRB kernel
#endif
#endif // 1D vs 2D versions



//==================================================================================================================================================
__global__ void __dot(subdomain_type * gpu_subdomains, int id_a, int id_b, double *global_dot, int level){

  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;


  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------
  double * grid_a = gpu_subdomains[box].levels[level].grids[id_a] + (plane+pencil+1);
  double * grid_b = gpu_subdomains[box].levels[level].grids[id_b] + (plane+pencil+1);


  // register/shared pipelines ---------------------------------------------------------------------------------------------------------------------
  __shared__ double local_dot[TBDIMY][TBDIMX];
  local_dot[threadIdx.y][threadIdx.x] = 0.0;
  __syncthreads(); 

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  // sweep through k dimension ---------------------------------------------------------------------------------------------------------------------
  for(k=0;k<subdomain_dim;k++){ // reduction in z
    int ijk = k*plane + j*pencil + i;
    if(withinBounds){
      local_dot[threadIdx.y][threadIdx.x]+=grid_a[ijk]*grid_b[ijk]; // dot product
    }
  } // for k
  __syncthreads();  // ensure everyone has found the max for their column
  int y;
  for(y=1;y<TBDIMY;y++){ // reduction in y
    if(threadIdx.y==0)local_dot[0][threadIdx.x]+=local_dot[y][threadIdx.x];
  }
  __syncthreads();  // ensure all threads[0][x] have have their sum over y
  int x;
  for(x=1;x<TBDIMX;x++){ // reduction in x
    if((threadIdx.x==0)&&(threadIdx.y==0))local_dot[0][0]+=local_dot[0][x];
  }
  __syncthreads();  // ensure thread 0 has thread block's sum;

  // now try and do an atomic add with local_dot[0][0] and *global_dot ------------------------------------------------------------------------------
  if((threadIdx.x==0)&&(threadIdx.y==0)){
    double oldd,newd;
    unsigned long long oldull, newull, CASrv;
    oldd = *global_dot;
    oldull = __double_as_longlong(oldd);
    newd = oldd;
    newd+=local_dot[0][0];
    newull = __double_as_longlong(newd);
    while ((CASrv=atomicCAS((unsigned long long *)global_dot, oldull, newull)) != oldull){
      oldull = CASrv;
      newd = __longlong_as_double(oldull);
      newd+=local_dot[0][0];
      newull = __double_as_longlong(newd);
    }
    // when complete, host should cudaMemcpy(&cpu_dot, gpu_dot, sizeof(double), cudaMemcpyDeviceToHost); 
  }

} // dot product kernel


//==================================================================================================================================================
__global__ void __sum(subdomain_type * gpu_subdomains, int id_a, double *global_sum, int level){

  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;


  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------
  double * grid_a = gpu_subdomains[box].levels[level].grids[id_a] + (plane+pencil+1);


  // register/shared pipelines ---------------------------------------------------------------------------------------------------------------------
  __shared__ double local_sum[TBDIMY][TBDIMX];
  local_sum[threadIdx.y][threadIdx.x] = 0.0;
  __syncthreads(); 

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  // sweep through k dimension ---------------------------------------------------------------------------------------------------------------------
  for(k=0;k<subdomain_dim;k++){ // reduction in z
    int ijk = k*plane + j*pencil + i;
    if(withinBounds){
      local_sum[threadIdx.y][threadIdx.x]+=grid_a[ijk];
    }
  } // for k
  __syncthreads();  // ensure everyone has found the max for their column
  int y;
  for(y=1;y<TBDIMY;y++){ // reduction in y
    if(threadIdx.y==0)local_sum[0][threadIdx.x]+=local_sum[y][threadIdx.x];
  }
  __syncthreads();  // ensure all threads[0][x] have have their sum over y
  int x;
  for(x=1;x<TBDIMX;x++){ // reduction in x
    if((threadIdx.x==0)&&(threadIdx.y==0))local_sum[0][0]+=local_sum[0][x];
  }
  __syncthreads();  // ensure thread 0 has thread block's sum;

  // now try and do an atomic add with local_sum[0][0] and *global_sum ------------------------------------------------------------------------------
  if((threadIdx.x==0)&&(threadIdx.y==0)){
    double oldd,newd;
    unsigned long long oldull, newull, CASrv;
    oldd = *global_sum;
    oldull = __double_as_longlong(oldd);
    newd = oldd;
    newd+=local_sum[0][0];
    newull = __double_as_longlong(newd);
    while ((CASrv=atomicCAS((unsigned long long *)global_sum, oldull, newull)) != oldull){
      oldull = CASrv;
      newd = __longlong_as_double(oldull);
      newd+=local_sum[0][0];
      newull = __double_as_longlong(newd);
    }
    // when complete, host should cudaMemcpy(&cpu_sum, gpu_sum, sizeof(double), cudaMemcpyDeviceToHost); 
  }

} // sum product kernel


//==================================================================================================================================================
__global__ void __norm(subdomain_type * gpu_subdomains, int grid_id, double *global_norm, int level){

  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;


  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------
  double *   grid = gpu_subdomains[box].levels[level].grids[grid_id] + (plane+pencil+1);


  // register/shared pipelines ---------------------------------------------------------------------------------------------------------------------
  __shared__ double max_norm[TBDIMY][TBDIMX];
  max_norm[threadIdx.y][threadIdx.x] = 0.0;
  __syncthreads(); 

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  // sweep through k dimension ---------------------------------------------------------------------------------------------------------------------
  for(k=0;k<subdomain_dim;k++){ // reduction in z
    int ijk = k*plane + j*pencil + i;
    double fabs_grid_ijk = fabs(grid[ijk]);
    if(withinBounds){
      if(fabs_grid_ijk>max_norm[threadIdx.y][threadIdx.x]){max_norm[threadIdx.y][threadIdx.x]=fabs_grid_ijk;} // max norm
    }
  } // for k
  __syncthreads();  // ensure everyone has found the max for their column
  int y;
  if(threadIdx.y==0)for(y=1;y<TBDIMY;y++){ // reduction in y
    if(max_norm[y][threadIdx.x]>max_norm[0][threadIdx.x]){max_norm[0][threadIdx.x]=max_norm[y][threadIdx.x];} // max norm
  }
  __syncthreads();  // ensure all threads in x have found the max for all their corresponding y's
  int x;
  if((threadIdx.x==0)&&(threadIdx.y==0))for(x=1;x<TBDIMX;x++){ // reduction in x
    if(max_norm[0][x]>max_norm[0][0]){max_norm[0][0]=max_norm[0][x];} // max norm
  }
  __syncthreads();  // ensure thread 0 has found the global max;

  // now try and do an atomic max with max_norm[0][0] and*global_norm ------------------------------------------------------------------------------
  if((threadIdx.x==0)&&(threadIdx.y==0)){
    double oldd,newd;
    unsigned long long oldull, newull, CASrv;
    oldd = *global_norm;
    oldull = __double_as_longlong(oldd);
    newd = oldd;if(max_norm[0][0]>newd)newd=max_norm[0][0];
    newull = __double_as_longlong(newd);
    while ((CASrv=atomicCAS((unsigned long long *)global_norm, oldull, newull)) != oldull){
      oldull = CASrv;
      newd = __longlong_as_double(oldull);
      if(max_norm[0][0]>newd)newd=max_norm[0][0];
      newull = __double_as_longlong(newd);
    }
    // when complete, host should cudaMemcpy(&cpu_norm, gpu_norm, sizeof(double), cudaMemcpyDeviceToHost); 
  }

} // norm kernel



//==================================================================================================================================================
#ifdef VL
__global__ void __residual(subdomain_type * gpu_subdomains, int res_id, int phi_id, int rhs_id, double a, double b, double h, int level){

  double h2inv = 1.0/(h*h);
  int box = blockIdx.z;  // CUDA 4 !!!

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  // |--pencil--|--IJStride--|--pencil--|--??--|
  int IJStride = (VL - pencil - pencil) & ~0x0F;
  int ij = ((IJStride*blockIdx.x + pencil - pencil) & ~0xF) + threadIdx.x; // i.e. shift vector so that thread0 is 128-byte aligned
  int k; 
 
  #ifdef __POINTERS_IN_SHARED
  __shared__ double * phi;
  __shared__ double * rhs;
  __shared__ double * alpha;
  __shared__ double * beta_i;
  __shared__ double * beta_j;
  __shared__ double * beta_k;
  __shared__ double * res;
  if(threadIdx.x==0){
       phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane);
       rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane);
     alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane);
    beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane);
    beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane);
    beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane);
       res = gpu_subdomains[box].levels[level].grids[  res_id] + (plane);
  }
  __syncthreads();
  #else
  double *         phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane);
  double *         rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane);
  double *       alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane);
  double *      beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane);
  double *      beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane);
  double *      beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane);
  double *         res = gpu_subdomains[box].levels[level].grids[  res_id] + (plane);
  #endif

 
  __shared__ double beta_i_ijk[VL];
  __shared__ double beta_j_ijk[VL];
  __shared__ double       temp[VL];
  double beta_k_ijk,beta_k_ijkp1;
  double    phi_ijk,   phi_ijkp1, phi_ijkm1;

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                                   int withinBounds = 1;
  if( threadIdx.x <  pencil           )withinBounds = 0;
  if( threadIdx.x >= pencil+IJStride  )withinBounds = 0;
  if( ij < pencil+1                   )withinBounds = 0;
  if( ij >=pencil*(subdomain_dim+1)-1 )withinBounds = 0;

  k = 0;
  int ijk = ij + k*plane;

     phi_ijkm1 =    phi[ijk-plane];
     phi_ijk   =    phi[ijk      ];
     phi_ijkp1 =    phi[ijk+plane];
  beta_k_ijk   = beta_k[ijk      ];
  beta_k_ijkp1 = beta_k[ijk+plane];


  for(k=0;k<subdomain_dim;k++){

    __syncthreads(); // WAR guard
                                                                      temp[threadIdx.x] = phi_ijk;
    if( (threadIdx.x>=pencil)&&(threadIdx.x<pencil+IJStride+1) )beta_i_ijk[threadIdx.x] = beta_i[ijk];
    if( (threadIdx.x>=pencil)                                  )beta_j_ijk[threadIdx.x] = beta_j[ijk];
    __syncthreads(); // RAW guard

    double helmholtz_ijk;
    if(withinBounds)
    helmholtz_ijk = a*alpha[ijk]*phi_ijk - b*h2inv*(
      beta_i_ijk[threadIdx.x+     1] * ( temp[threadIdx.x+     1] - phi_ijk                  ) -
      beta_i_ijk[threadIdx.x       ] * ( phi_ijk                  - temp[threadIdx.x-     1] ) +
      beta_j_ijk[threadIdx.x+pencil] * ( temp[threadIdx.x+pencil] - phi_ijk                  ) -
      beta_j_ijk[threadIdx.x       ] * ( phi_ijk                  - temp[threadIdx.x-pencil] ) +
      beta_k_ijkp1                   * ( phi_ijkp1                - phi_ijk                  ) -
      beta_k_ijk                     * ( phi_ijk                  - phi_ijkm1                )
    );

    if(withinBounds)res[ijk] = rhs[ijk] - helmholtz_ijk;

    ijk+=plane;
    // rotate register pipeline...
    phi_ijkm1    = phi_ijk;
    phi_ijk      = phi_ijkp1;
    phi_ijkp1    = phi[ijk+plane];
    beta_k_ijk   = beta_k_ijkp1;
    if(withinBounds)
    beta_k_ijkp1 = beta_k[ijk+plane];

  } // for k
} // residual kernel

#else
#ifdef __LOCALITY_VIA_SHARED
//==============================================================================================================================================================
// shared+L1 version (relies on shared memory but you must still favorL1)
//==============================================================================================================================================================
__global__ void __residual(subdomain_type * gpu_subdomains, int res_id, int phi_id, int rhs_id, double a, double b, double h, int level){

  double h2inv = 1.0/(h*h);
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------------------
  #ifdef __POINTERS_IN_SHARED
  __shared__ double * phi;
  __shared__ double * rhs;
  __shared__ double * alpha;
  __shared__ double * beta_i;
  __shared__ double * beta_j;
  __shared__ double * beta_k;
  __shared__ double * res;
  if((threadIdx.x==0)&&(threadIdx.y==0)){
       phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane+pencil+1);
       rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane+pencil+1);
     alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
    beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
    beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
    beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
       res = gpu_subdomains[box].levels[level].grids[  res_id] + (plane+pencil+1);
  }
  __syncthreads();
  #else
  double *         phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane+pencil+1);
  double *         rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane+pencil+1);
  double *       alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
  double *      beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
  double *      beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
  double *      beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
  double *         res = gpu_subdomains[box].levels[level].grids[  res_id] + (plane+pencil+1);
  #endif


  __shared__ double beta_i_ijk[TBDIMY  ][TBDIMX+1];
  __shared__ double beta_j_ijk[TBDIMY+1][TBDIMX  ];
  __shared__ double       temp[TBDIMY+2][TBDIMX+2]; // always index with threadIdx.y+1,x+1
  double beta_k_ijk,beta_k_ijkp1;
  double    phi_ijk,   phi_ijkp1, phi_ijkm1;
  

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

     k = 0; 
     int ijk = k*plane + j*pencil + i;

     phi_ijkm1 = phi[ijk-plane];
     phi_ijk   = phi[ijk      ];
     phi_ijkp1 = phi[ijk+plane];
  beta_k_ijk   = beta_k[ijk      ];
  beta_k_ijkp1 = beta_k[ijk+plane];

  for(k=0;k<subdomain_dim;k++){

    __syncthreads(); // WAR guard
    // copy phi and beta_i/j to shared temp array...
                            temp[threadIdx.y+1][threadIdx.x+1] = phi_ijk;
    if(threadIdx.x==0)      temp[threadIdx.y+1][            0] = phi[ijk-1];
    if(threadIdx.y==0)      temp[            0][threadIdx.x+1] = phi[ijk-pencil];
    if(threadIdx.x==0)      temp[threadIdx.y+1][     TBDIMX+1] = phi[ijk+TBDIMX];
    if(threadIdx.y==0)      temp[     TBDIMY+1][threadIdx.x+1] = phi[ijk+TBDIMY*pencil];
                      beta_i_ijk[threadIdx.y  ][threadIdx.x  ] = beta_i[ijk];
    if(threadIdx.x==0)beta_i_ijk[threadIdx.y  ][     TBDIMX  ] = beta_i[ijk+TBDIMX];
                      beta_j_ijk[threadIdx.y  ][threadIdx.x  ] = beta_j[ijk];
    if(threadIdx.y==0)beta_j_ijk[     TBDIMY  ][threadIdx.x  ] = beta_j[ijk+TBDIMY*pencil];
    __syncthreads(); // RAW guard

    double helmholtz_ijk = a*alpha[ijk]*phi_ijk - b*h2inv*( 
      beta_i_ijk[threadIdx.y  ][threadIdx.x+1] * ( temp[threadIdx.y+1][threadIdx.x+2] - phi_ijk                            ) - 
      beta_i_ijk[threadIdx.y  ][threadIdx.x  ] * ( phi_ijk                            - temp[threadIdx.y+1][threadIdx.x  ] ) +
      beta_j_ijk[threadIdx.y+1][threadIdx.x  ] * ( temp[threadIdx.y+2][threadIdx.x+1] - phi_ijk                            ) -
      beta_j_ijk[threadIdx.y  ][threadIdx.x  ] * ( phi_ijk                            - temp[threadIdx.y  ][threadIdx.x+1] ) +
      beta_k_ijkp1                             * ( phi_ijkp1                          - phi_ijk                            ) -
      beta_k_ijk                               * ( phi_ijk                            - phi_ijkm1                          )
    );
 
    if(withinBounds)res[ijk] = rhs[ijk] - helmholtz_ijk;

    ijk+=plane;
    // rotate register pipeline...
    phi_ijkm1    = phi_ijk;
    phi_ijk      = phi_ijkp1;
    phi_ijkp1    = phi[ijk+plane];
    beta_k_ijk   = beta_k_ijkp1;
    beta_k_ijkp1 = beta_k[ijk+plane];

  } // for k
} // residual kernel
#else
//==============================================================================================================================================================
// cache version (relies solely on the L1/L2 cache hierarchy)
//==============================================================================================================================================================
__global__ void __residual(subdomain_type * gpu_subdomains, int res_id, int phi_id, int rhs_id, double a, double b, double h, int level){

  double h2inv = 1.0/(h*h);
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;


  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------------------
  #ifdef __POINTERS_IN_SHARED
  __shared__ double * phi;
  __shared__ double * rhs;
  __shared__ double * alpha;
  __shared__ double * beta_i;
  __shared__ double * beta_j;
  __shared__ double * beta_k;
  __shared__ double * res;
  if((threadIdx.x==0)&&(threadIdx.y==0)){
       phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane+pencil+1);
       rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane+pencil+1);
     alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
    beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
    beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
    beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
       res = gpu_subdomains[box].levels[level].grids[  res_id] + (plane+pencil+1);
  }
  __syncthreads();
  #else
  double *         phi = gpu_subdomains[box].levels[level].grids[  phi_id] + (plane+pencil+1);
  double *         rhs = gpu_subdomains[box].levels[level].grids[  rhs_id] + (plane+pencil+1);
  double *       alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
  double *      beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
  double *      beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
  double *      beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
  double *         res = gpu_subdomains[box].levels[level].grids[  res_id] + (plane+pencil+1);
  #endif



  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  for(k=0;k<subdomain_dim;k++){
    int ijk = k*plane + j*pencil + i;
    double helmholtz_ijk = a*alpha[ijk]*phi[ijk] - b*h2inv*(
      beta_i[ijk     +1] * ( phi[ijk     +1] - phi[ijk       ] ) -
      beta_i[ijk       ] * ( phi[ijk       ] - phi[ijk     -1] ) +
      beta_j[ijk+pencil] * ( phi[ijk+pencil] - phi[ijk       ] ) -
      beta_j[ijk       ] * ( phi[ijk       ] - phi[ijk-pencil] ) +
      beta_k[ijk +plane] * ( phi[ijk +plane] - phi[ijk       ] ) -
      beta_k[ijk       ] * ( phi[ijk       ] - phi[ijk -plane] )
    );
    if(withinBounds)res[ijk] = rhs[ijk] - helmholtz_ijk;
  } // for k
} // residual kernel
#endif
#endif // VL

//=============================================================================================================================================================
#ifdef VL
__global__ void __apply_op(subdomain_type * gpu_subdomains, int  Ax_id, int x_id, double a, double b, double h, int level){

  double h2inv = 1.0/(h*h);
  int box = blockIdx.z;  // CUDA 4 !!!

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  // |--pencil--|--IJStride--|--pencil--|--??--|
  int IJStride = (VL - pencil - pencil) & ~0x0F;
  int ij = ((IJStride*blockIdx.x + pencil - pencil) & ~0xF) + threadIdx.x; // i.e. shift vector so that thread0 is 128-byte aligned
  int k; 
 
  #ifdef __POINTERS_IN_SHARED
  __shared__ double * x;
  __shared__ double * alpha;
  __shared__ double * beta_i;
  __shared__ double * beta_j;
  __shared__ double * beta_k;
  __shared__ double *  Ax;
  if(threadIdx.x==0){
         x = gpu_subdomains[box].levels[level].grids[    x_id] + (plane);
     alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane);
    beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane);
    beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane);
    beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane);
        Ax = gpu_subdomains[box].levels[level].grids[   Ax_id] + (plane);
  }
  __syncthreads();
  #else
  double *           x = gpu_subdomains[box].levels[level].grids[    x_id] + (plane);
  double *       alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane);
  double *      beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane);
  double *      beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane);
  double *      beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane);
  double *          Ax = gpu_subdomains[box].levels[level].grids[   Ax_id] + (plane);
  #endif

 
  __shared__ double beta_i_ijk[VL];
  __shared__ double beta_j_ijk[VL];
  __shared__ double       temp[VL];
  double beta_k_ijk,beta_k_ijkp1;
  double      x_ijk,     x_ijkp1,   x_ijkm1;

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                                   int withinBounds = 1;
  if( threadIdx.x <  pencil           )withinBounds = 0;
  if( threadIdx.x >= pencil+IJStride  )withinBounds = 0;
  if( ij < pencil+1                   )withinBounds = 0;
  if( ij >=pencil*(subdomain_dim+1)-1 )withinBounds = 0;

  k = 0;
  int ijk = ij + k*plane;

       x_ijkm1 =      x[ijk-plane];
       x_ijk   =      x[ijk      ];
       x_ijkp1 =      x[ijk+plane];
  beta_k_ijk   = beta_k[ijk      ];
  beta_k_ijkp1 = beta_k[ijk+plane];


  for(k=0;k<subdomain_dim;k++){

    __syncthreads(); // WAR guard
                                                                      temp[threadIdx.x] =   x_ijk;
    if( (threadIdx.x>=pencil)&&(threadIdx.x<pencil+IJStride+1) )beta_i_ijk[threadIdx.x] = beta_i[ijk];
    if( (threadIdx.x>=pencil)                                  )beta_j_ijk[threadIdx.x] = beta_j[ijk];
    __syncthreads(); // RAW guard

    double helmholtz_ijk;
    if(withinBounds)
    helmholtz_ijk = a*alpha[ijk]*  x_ijk - b*h2inv*(
      beta_i_ijk[threadIdx.x+     1] * ( temp[threadIdx.x+     1] -   x_ijk                  ) -
      beta_i_ijk[threadIdx.x       ] * (   x_ijk                  - temp[threadIdx.x-     1] ) +
      beta_j_ijk[threadIdx.x+pencil] * ( temp[threadIdx.x+pencil] -   x_ijk                  ) -
      beta_j_ijk[threadIdx.x       ] * (   x_ijk                  - temp[threadIdx.x-pencil] ) +
      beta_k_ijkp1                   * (   x_ijkp1                -   x_ijk                  ) -
      beta_k_ijk                     * (   x_ijk                  -   x_ijkm1                )
    );

    if(withinBounds) Ax[ijk] = helmholtz_ijk;

    ijk+=plane;
    // rotate register pipeline...
      x_ijkm1    =   x_ijk;
      x_ijk      =   x_ijkp1;
      x_ijkp1    =   x[ijk+plane];
    beta_k_ijk   = beta_k_ijkp1;
    if(withinBounds)
    beta_k_ijkp1 = beta_k[ijk+plane];

  } // for k
} // residual kernel

#else
#ifdef __LOCALITY_VIA_SHARED
//==============================================================================================================================================================
// shared+L1 version (relies on shared memory but you must still favorL1)
//==============================================================================================================================================================
__global__ void __apply_op(subdomain_type * gpu_subdomains, int  Ax_id, int   x_id, double a, double b, double h, int level){

  double h2inv = 1.0/(h*h);
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------------------
  #ifdef __POINTERS_IN_SHARED
  __shared__ double *   x;
  __shared__ double * alpha;
  __shared__ double * beta_i;
  __shared__ double * beta_j;
  __shared__ double * beta_k;
  __shared__ double *  Ax;
  if((threadIdx.x==0)&&(threadIdx.y==0)){
         x = gpu_subdomains[box].levels[level].grids[    x_id] + (plane+pencil+1);
     alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
    beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
    beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
    beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
        Ax = gpu_subdomains[box].levels[level].grids[   Ax_id] + (plane+pencil+1);
  }
  __syncthreads();
  #else
  double *           x = gpu_subdomains[box].levels[level].grids[    x_id] + (plane+pencil+1);
  double *       alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
  double *      beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
  double *      beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
  double *      beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
  double *          Ax = gpu_subdomains[box].levels[level].grids[   Ax_id] + (plane+pencil+1);
  #endif


  __shared__ double beta_i_ijk[TBDIMY  ][TBDIMX+1];
  __shared__ double beta_j_ijk[TBDIMY+1][TBDIMX  ];
  __shared__ double       temp[TBDIMY+2][TBDIMX+2]; // always index with threadIdx.y+1,x+1
  double beta_k_ijk,beta_k_ijkp1;
  double      x_ijk,     x_ijkp1,   x_ijkm1;
  

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

     k = 0; 
     int ijk = k*plane + j*pencil + i;

       x_ijkm1 =   x[ijk-plane];
       x_ijk   =   x[ijk      ];
       x_ijkp1 =   x[ijk+plane];
  beta_k_ijk   = beta_k[ijk      ];
  beta_k_ijkp1 = beta_k[ijk+plane];

  for(k=0;k<subdomain_dim;k++){

    __syncthreads(); // WAR guard
    // copy   x and beta_i/j to shared temp array...
                            temp[threadIdx.y+1][threadIdx.x+1] =   x_ijk;
    if(threadIdx.x==0)      temp[threadIdx.y+1][            0] =   x[ijk-1];
    if(threadIdx.y==0)      temp[            0][threadIdx.x+1] =   x[ijk-pencil];
    if(threadIdx.x==0)      temp[threadIdx.y+1][     TBDIMX+1] =   x[ijk+TBDIMX];
    if(threadIdx.y==0)      temp[     TBDIMY+1][threadIdx.x+1] =   x[ijk+TBDIMY*pencil];
                      beta_i_ijk[threadIdx.y  ][threadIdx.x  ] = beta_i[ijk];
    if(threadIdx.x==0)beta_i_ijk[threadIdx.y  ][     TBDIMX  ] = beta_i[ijk+TBDIMX];
                      beta_j_ijk[threadIdx.y  ][threadIdx.x  ] = beta_j[ijk];
    if(threadIdx.y==0)beta_j_ijk[     TBDIMY  ][threadIdx.x  ] = beta_j[ijk+TBDIMY*pencil];
    __syncthreads(); // RAW guard

    double helmholtz_ijk = a*alpha[ijk]*  x_ijk - b*h2inv*( 
      beta_i_ijk[threadIdx.y  ][threadIdx.x+1] * ( temp[threadIdx.y+1][threadIdx.x+2] -   x_ijk                            ) - 
      beta_i_ijk[threadIdx.y  ][threadIdx.x  ] * (   x_ijk                            - temp[threadIdx.y+1][threadIdx.x  ] ) +
      beta_j_ijk[threadIdx.y+1][threadIdx.x  ] * ( temp[threadIdx.y+2][threadIdx.x+1] -   x_ijk                            ) -
      beta_j_ijk[threadIdx.y  ][threadIdx.x  ] * (   x_ijk                            - temp[threadIdx.y  ][threadIdx.x+1] ) +
      beta_k_ijkp1                             * (   x_ijkp1                          -   x_ijk                            ) -
      beta_k_ijk                               * (   x_ijk                            -   x_ijkm1                          )
    );
 
    if(withinBounds) Ax[ijk] = helmholtz_ijk;

    ijk+=plane;
    // rotate register pipeline...
      x_ijkm1    =   x_ijk;
      x_ijk      =   x_ijkp1;
      x_ijkp1    =   x[ijk+plane];
    beta_k_ijk   = beta_k_ijkp1;
    beta_k_ijkp1 = beta_k[ijk+plane];

  } // for k
} // residual kernel
#else
//==============================================================================================================================================================
// cache version (relies solely on the L1/L2 cache hierarchy)
//==============================================================================================================================================================
__global__ void __apply_op(subdomain_type * gpu_subdomains, int  Ax_id, int   x_id, double a, double b, double h, int level){

  double h2inv = 1.0/(h*h);
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;


  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------------------
  #ifdef __POINTERS_IN_SHARED
  __shared__ double *   x;
  __shared__ double * alpha;
  __shared__ double * beta_i;
  __shared__ double * beta_j;
  __shared__ double * beta_k;
  __shared__ double *  Ax;
  if((threadIdx.x==0)&&(threadIdx.y==0)){
         x = gpu_subdomains[box].levels[level].grids[    x_id] + (plane+pencil+1);
     alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
    beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
    beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
    beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
        Ax = gpu_subdomains[box].levels[level].grids[   Ax_id] + (plane+pencil+1);
  }
  __syncthreads();
  #else
  double *           x = gpu_subdomains[box].levels[level].grids[    x_id] + (plane+pencil+1);
  double *       alpha = gpu_subdomains[box].levels[level].grids[ __alpha] + (plane+pencil+1);
  double *      beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (plane+pencil+1);
  double *      beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (plane+pencil+1);
  double *      beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (plane+pencil+1);
  double *          Ax = gpu_subdomains[box].levels[level].grids[   Ax_id] + (plane+pencil+1);
  #endif



  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  for(k=0;k<subdomain_dim;k++){
    int ijk = k*plane + j*pencil + i;
    double helmholtz_ijk = a*alpha[ijk]*  x[ijk] - b*h2inv*(
      beta_i[ijk     +1] * (   x[ijk     +1] -   x[ijk       ] ) -
      beta_i[ijk       ] * (   x[ijk       ] -   x[ijk     -1] ) +
      beta_j[ijk+pencil] * (   x[ijk+pencil] -   x[ijk       ] ) -
      beta_j[ijk       ] * (   x[ijk       ] -   x[ijk-pencil] ) +
      beta_k[ijk +plane] * (   x[ijk +plane] -   x[ijk       ] ) -
      beta_k[ijk       ] * (   x[ijk       ] -   x[ijk -plane] )
    );
    if(withinBounds) Ax[ijk] = helmholtz_ijk;
  } // for k
} // residual kernel
#endif
#endif // VL

//=============================================================================================================================================================
__global__ void __restriction_betas(subdomain_type * gpu_subdomains, int fine_level){

  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int   fine_dim    = gpu_subdomains[box].levels[fine_level  ].dim.i;
  int   fine_pencil = gpu_subdomains[box].levels[fine_level  ].pencil;
  int   fine_plane  = gpu_subdomains[box].levels[fine_level  ].plane;
  int coarse_dim    = gpu_subdomains[box].levels[fine_level+1].dim.i;
  int coarse_pencil = gpu_subdomains[box].levels[fine_level+1].pencil;
  int coarse_plane  = gpu_subdomains[box].levels[fine_level+1].plane;

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                    int withinCoarseBounds = 1;
  if((i>>1)>=coarse_dim)withinCoarseBounds = 0;
  if((j>>1)>=coarse_dim)withinCoarseBounds = 0;

  // restrict beta_i  (== face in jk)
  double * beta_f = gpu_subdomains[box].levels[fine_level  ].grids[__beta_i] + (  fine_plane+  fine_pencil+1);
  double * beta_c = gpu_subdomains[box].levels[fine_level+1].grids[__beta_i] + (coarse_plane+coarse_pencil+1);
  for(k=0;k<fine_dim;k+=2){
    int   fine_ijk = (k   )*  fine_plane + (j   )*  fine_pencil + (i   );
    int coarse_ijk = (k>>1)*coarse_plane + (j>>1)*coarse_pencil + (i>>1);
    if(withinCoarseBounds){
    if(!((threadIdx.x|threadIdx.y)&0x1)){ // i.e. x and y LSB are 0
      beta_c[coarse_ijk] = ( beta_f[fine_ijk            ]+beta_f[fine_ijk+fine_pencil           ] +
                             beta_f[fine_ijk+fine_plane ]+beta_f[fine_ijk+fine_pencil+fine_plane] ) * 0.25;
    }}
  } // for k
  // restrict beta_j  (== face in ik)
  beta_f = gpu_subdomains[box].levels[fine_level  ].grids[__beta_j] + (  fine_plane+  fine_pencil+1);
  beta_c = gpu_subdomains[box].levels[fine_level+1].grids[__beta_j] + (coarse_plane+coarse_pencil+1);
  for(k=0;k<fine_dim;k+=2){
    int   fine_ijk = (k   )*  fine_plane + (j   )*  fine_pencil + (i   );
    int coarse_ijk = (k>>1)*coarse_plane + (j>>1)*coarse_pencil + (i>>1);
    if(withinCoarseBounds){
    if(!((threadIdx.x|threadIdx.y)&0x1)){ // i.e. x and y LSB are 0
      beta_c[coarse_ijk] = ( beta_f[fine_ijk            ]+beta_f[fine_ijk+1                     ] +
                             beta_f[fine_ijk+fine_plane ]+beta_f[fine_ijk+1+fine_plane          ] ) * 0.25;
    }}
  } // for k
  // restrict beta_k  (== face in ij)
  beta_f = gpu_subdomains[box].levels[fine_level  ].grids[__beta_k] + (  fine_plane+  fine_pencil+1);
  beta_c = gpu_subdomains[box].levels[fine_level+1].grids[__beta_k] + (coarse_plane+coarse_pencil+1);
  for(k=0;k<fine_dim;k+=2){
    int   fine_ijk = (k   )*  fine_plane + (j   )*  fine_pencil + (i   );
    int coarse_ijk = (k>>1)*coarse_plane + (j>>1)*coarse_pencil + (i>>1);
    if(withinCoarseBounds){
    if(!((threadIdx.x|threadIdx.y)&0x1)){ // i.e. x and y LSB are 0
      beta_c[coarse_ijk] = ( beta_f[fine_ijk            ]+beta_f[fine_ijk+1                     ] +
                             beta_f[fine_ijk+fine_pencil]+beta_f[fine_ijk+1+fine_pencil         ] ) * 0.25;
    }}
  } // for k
} // restriction of face-centered coefficients kernel

//=============================================================================================================================================================
__global__ void __restriction(subdomain_type * gpu_subdomains, int fine_id, int coarse_id, int fine_level){

  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int   fine_dim    = gpu_subdomains[box].levels[fine_level  ].dim.i;
  int   fine_pencil = gpu_subdomains[box].levels[fine_level  ].pencil;
  int   fine_plane  = gpu_subdomains[box].levels[fine_level  ].plane;
  int coarse_dim    = gpu_subdomains[box].levels[fine_level+1].dim.i;
  int coarse_pencil = gpu_subdomains[box].levels[fine_level+1].pencil;
  int coarse_plane  = gpu_subdomains[box].levels[fine_level+1].plane;

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                    int withinCoarseBounds = 1;
  if((i>>1)>=coarse_dim)withinCoarseBounds = 0;
  if((j>>1)>=coarse_dim)withinCoarseBounds = 0;

  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------------------
  double *   fine = gpu_subdomains[box].levels[fine_level  ].grids[  fine_id] + (  fine_plane+  fine_pencil+1);
  double * coarse = gpu_subdomains[box].levels[fine_level+1].grids[coarse_id] + (coarse_plane+coarse_pencil+1);

  // FIX - make coarse-oriented, not fine oriented...
  for(k=0;k<fine_dim;k+=2){
    int   fine_ijk = (k   )*  fine_plane + (j   )*  fine_pencil + (i   );
    int coarse_ijk = (k>>1)*coarse_plane + (j>>1)*coarse_pencil + (i>>1);

    if(withinCoarseBounds){
    if(!((threadIdx.x|threadIdx.y)&0x1)){ // i.e. x and y LSB are 0
      coarse[coarse_ijk] = 0.125 * (
                                   fine[fine_ijk                         ] +
                                   fine[fine_ijk                       +1] +
                                   fine[fine_ijk           +fine_pencil  ] +
                                   fine[fine_ijk           +fine_pencil+1] +
                                   fine[fine_ijk+fine_plane              ] +
                                   fine[fine_ijk+fine_plane            +1] +
                                   fine[fine_ijk+fine_plane+fine_pencil  ] +
                                   fine[fine_ijk+fine_plane+fine_pencil+1]
                                   );
    }}
  } // for k
} // restriction kernel

//=============================================================================================================================================================
__global__ void __interpolation(subdomain_type * gpu_subdomains, int coarse_id, int fine_id, int fine_level){

  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int   fine_dim    = gpu_subdomains[box].levels[fine_level  ].dim.i;
  int   fine_pencil = gpu_subdomains[box].levels[fine_level  ].pencil;
  int   fine_plane  = gpu_subdomains[box].levels[fine_level  ].plane;
  int coarse_pencil = gpu_subdomains[box].levels[fine_level+1].pencil;
  int coarse_plane  = gpu_subdomains[box].levels[fine_level+1].plane;

  // construct pointers to element (0,0,0) in each array for the current subdomain -----------------------------------------------------------------------------
  double *   fine = gpu_subdomains[box].levels[fine_level  ].grids[  fine_id] + (  fine_plane+  fine_pencil+1);
  double * coarse = gpu_subdomains[box].levels[fine_level+1].grids[coarse_id] + (coarse_plane+coarse_pencil+1);

  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                    int   withinFineBounds = 1;
  if( i    >=  fine_dim)  withinFineBounds = 0;
  if( j    >=  fine_dim)  withinFineBounds = 0;
//                  int withinCoarseBounds = 1;
//if((i>>1)>=coarse_dim)withinCoarseBounds = 0;
//if((j>>1)>=coarse_dim)withinCoarseBounds = 0;

  for(k=0;k<fine_dim;k+=2){
    int   fine_ijk = (k   )*  fine_plane + (j   )*  fine_pencil + (i   );
    int coarse_ijk = (k>>1)*coarse_plane + (j>>1)*coarse_pencil + (i>>1);

    if(withinFineBounds){
      fine[fine_ijk           ] += coarse[coarse_ijk];
      fine[fine_ijk+fine_plane] += coarse[coarse_ijk];
    }
   
  } // for k
} // interpolation kernel

//==============================================================================================================================================================
__global__ void __grid_to_surface_buffers(subdomain_type * gpu_subdomains, int grid_id, int level){

  // ASSUME: ThreadBlock=(Dim x 1 x 1), Grid=(1,Dim,subdomains)

  //  |\
  //  |  \
  //  |\   \
  //  |  \   \
  //  |\   \  |   
  //  |  \   \|   ^ blockIdx.y (blockIdx.x is unused for grid dims<512)
  //  |\   \  |   |
  //  |  \   \|   |
  //   \   \  |   |
  //     \   \|
  //       \  |
  //         \|
  //   ---->
  //  threadIdx.x (threadIdx.y is unused)

  // FIX, what if ghosts > 1 ???

  int box = blockIdx.z;  // CUDA 4 !!!


  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  // do the 6 faces...
  int di,dj,dk;
  for(dk=-1;dk<=1;dk++){
  for(dj=-1;dj<=1;dj++){
  for(di=-1;di<=1;di++){int n=13+di+3*dj+9*dk;if(faces[n]){
    int low_i,low_j,low_k;
    int elementStride;
    int  pencilStride;
    switch(di){
      case -1:low_i=  1;break;
      case  0:low_i=  1;break;
      case  1:low_i=subdomain_dim;break;
    };
    switch(dj){
      case -1:low_j=  1;break;
      case  0:low_j=  1;break;
      case  1:low_j=subdomain_dim;break;
    };
    switch(dk){
      case -1:low_k=  1;break;
      case  0:low_k=  1;break;
      case  1:low_k=subdomain_dim;break;
    };
    if(di==0)elementStride =      1; // ij and ik faces
        else elementStride = pencil; // jk faces
    if(dk==0) pencilStride =  plane; // ik and jk faces
        else  pencilStride = pencil; // ij faces
    int ijk = low_k*plane + low_j*pencil + low_i + pencilStride*blockIdx.y + elementStride*threadIdx.x;
    int b = subdomain_dim*blockIdx.y + threadIdx.x;

    double * surface_buf = gpu_subdomains[box].levels[level].surface_bufs[n];
    double * grid        = gpu_subdomains[box].levels[level].grids[grid_id];
 
    #ifdef __USE_LDG
    surface_buf[b] = __ldg(grid+ijk);
    #else
    surface_buf[b] = grid[ijk];
    #endif
  }}}}
}
//==============================================================================================================================================================
__global__ void __ghost_buffers_to_grid(subdomain_type * gpu_subdomains, int grid_id, int level){


  int box = blockIdx.z;  // CUDA 4 !!!


  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  // do the 6 faces...
  int di,dj,dk;
  for(dk=-1;dk<=1;dk++){
  for(dj=-1;dj<=1;dj++){
  for(di=-1;di<=1;di++){int n=13+di+3*dj+9*dk;if(faces[n]){
    int low_i,low_j,low_k;
    int elementStride;
    int  pencilStride;
    switch(di){
      case -1:low_i=              0;break;
      case  0:low_i=              1;break;
      case  1:low_i=subdomain_dim+1;break;
    };
    switch(dj){
      case -1:low_j=              0;break;
      case  0:low_j=              1;break;
      case  1:low_j=subdomain_dim+1;break;
    };
    switch(dk){
      case -1:low_k=              0;break;
      case  0:low_k=              1;break;
      case  1:low_k=subdomain_dim+1;break;
    };
    if(di==0)elementStride =      1; // ij and ik faces
        else elementStride = pencil; // jk faces
    if(dk==0) pencilStride =  plane; // ik and jk faces
        else  pencilStride = pencil; // ij faces
    int ijk = low_k*plane + low_j*pencil + low_i + pencilStride*blockIdx.y + elementStride*threadIdx.x;
    int b = subdomain_dim*blockIdx.y + threadIdx.x;

    double * ghost_buf = gpu_subdomains[box].levels[level].ghost_bufs[n];
    double * grid      = gpu_subdomains[box].levels[level].grids[grid_id];
 
    grid[ijk] = ghost_buf[b];
  }}}}
}


//==============================================================================================================================================================
__global__ void __surface_buffers_to_ghost_buffers(subdomain_type * gpu_subdomains, int grid_id, int level, int MyRank){
  int recvBox = blockIdx.z;  // CUDA 4 !!!
  int subdomain_dim = gpu_subdomains[recvBox].levels[level].dim.i;
  int n;
  for(n=0;n<27;n++)
  if( faces[n] && (gpu_subdomains[recvBox].neighbors[n].rank == MyRank) ){
    int sendBox = gpu_subdomains[recvBox].neighbors[n].local_index;
    double *   ghost_buf = gpu_subdomains[recvBox].levels[level].ghost_bufs[n];
    double * surface_buf = gpu_subdomains[sendBox].levels[level].surface_bufs[26-n];
    int b = subdomain_dim*blockIdx.y + threadIdx.x;
    ghost_buf[b] = surface_buf[b];
  }
}

//==============================================================================================================================================================
#ifdef __MPI
__global__ void __surface_buffers_to_send_buffer(subdomain_type * gpu_subdomains, double ** gpu_send_buffer, double ** gpu_recv_buffer, int grid_id, int level, int MyRank){
  int sendBox = blockIdx.z;  // CUDA 4 !!!
  int subdomain_dim = gpu_subdomains[sendBox].levels[level].dim.i;
  int n;
  int FaceSizeAtLevel = subdomain_dim*subdomain_dim;
  for(n=0;n<27;n++)
  if( faces[n] &&   (gpu_subdomains[sendBox].neighbors[n].rank != MyRank) ){
    int        buf = gpu_subdomains[sendBox].neighbors[n].send.buf;
    int FaceOffset = gpu_subdomains[sendBox].neighbors[n].send.offset.faces;
    double * surface_buf = gpu_subdomains[sendBox].levels[level].surface_bufs[n];
    double *    send_buf = gpu_send_buffer[buf] + FaceSizeAtLevel*FaceOffset;
    int b = subdomain_dim*blockIdx.y + threadIdx.x;
    send_buf[b] = surface_buf[b];
  }
}

__global__ void __recv_buffer_to_ghost_buffers(subdomain_type * gpu_subdomains, double ** gpu_send_buffer, double ** gpu_recv_buffer, int grid_id, int level, int MyRank){
  int recvBox = blockIdx.z;  // CUDA 4 !!!
  int subdomain_dim = gpu_subdomains[recvBox].levels[level].dim.i;
  int n;
  int FaceSizeAtLevel = subdomain_dim*subdomain_dim;
  for(n=0;n<27;n++)
  if( faces[n] &&   (gpu_subdomains[recvBox].neighbors[n].rank != MyRank) ){
    int        buf = gpu_subdomains[recvBox].neighbors[n].recv.buf;
    int FaceOffset = gpu_subdomains[recvBox].neighbors[n].recv.offset.faces;
    double *   ghost_buf = gpu_subdomains[recvBox].levels[level].ghost_bufs[n];
    double *    recv_buf = gpu_recv_buffer[buf] + FaceSizeAtLevel*FaceOffset;
    int b = subdomain_dim*blockIdx.y + threadIdx.x;
    ghost_buf[b] = recv_buf[b];
  }
}
#endif

 
//==============================================================================================================================================================
__global__ void __initialize_grid_to_scalar(subdomain_type * gpu_subdomains, int grid_id, double h, double value, int level){
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  double * grid = gpu_subdomains[box].levels[level].grids[grid_id] + (plane+pencil+1);
  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  k = 0; 
  int ijk = k*plane + j*pencil + i;

  for(k=0;k<subdomain_dim;k++){
    if(withinBounds)grid[ijk] = value;
    ijk+=plane;
  } // for k
} // iniitalize kernel


//==============================================================================================================================================================
__global__ void __scale_grid(subdomain_type * gpu_subdomains, int id_c, double scale_a, int id_a, int level){ // c=scale_a*id_a
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  double * grid_a = gpu_subdomains[box].levels[level].grids[id_a] + (plane+pencil+1);
  double * grid_c = gpu_subdomains[box].levels[level].grids[id_c] + (plane+pencil+1);
  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  k = 0; 
  int ijk = k*plane + j*pencil + i;

  for(k=0;k<subdomain_dim;k++){
    if(withinBounds)grid_c[ijk] = scale_a*grid_a[ijk];
    ijk+=plane;
  } // for k
} // scale_grid kernel


//==============================================================================================================================================================
__global__ void __mul_grids(subdomain_type * gpu_subdomains, int id_c, double scale, int id_a, int id_b, int level){ // c=scale*id_a*id_b
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  double * grid_a = gpu_subdomains[box].levels[level].grids[id_a] + (plane+pencil+1);
  double * grid_b = gpu_subdomains[box].levels[level].grids[id_b] + (plane+pencil+1);
  double * grid_c = gpu_subdomains[box].levels[level].grids[id_c] + (plane+pencil+1);
  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  k = 0; 
  int ijk = k*plane + j*pencil + i;

  for(k=0;k<subdomain_dim;k++){
    if(withinBounds)grid_c[ijk] = scale*grid_a[ijk]*grid_b[ijk];
    ijk+=plane;
  } // for k
} // mul_grids kernel


//==============================================================================================================================================================
__global__ void __shift_grid(subdomain_type * gpu_subdomains, int id_c, int id_a, double shift_a, int level){ // id_c=id_a + shift_a
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  double * grid_c = gpu_subdomains[box].levels[level].grids[id_c] + (plane+pencil+1);
  double * grid_a = gpu_subdomains[box].levels[level].grids[id_a] + (plane+pencil+1);
  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  k = 0; 
  int ijk = k*plane + j*pencil + i;

  for(k=0;k<subdomain_dim;k++){
    if(withinBounds)grid_c[ijk] = grid_a[ijk] + shift_a;
    ijk+=plane;
  } // for k
} // shift_grid kernel


//==============================================================================================================================================================
__global__ void __add_grids(subdomain_type * gpu_subdomains, int id_c, double scale_a, int id_a, double scale_b, int id_b, int level){ // c=scale_a*id_a + scale_b*id_b
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  double * grid_a = gpu_subdomains[box].levels[level].grids[id_a] + (plane+pencil+1);
  double * grid_b = gpu_subdomains[box].levels[level].grids[id_b] + (plane+pencil+1);
  double * grid_c = gpu_subdomains[box].levels[level].grids[id_c] + (plane+pencil+1);
  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  k = 0; 
  int ijk = k*plane + j*pencil + i;

  for(k=0;k<subdomain_dim;k++){
    if(withinBounds)grid_c[ijk] = scale_a*grid_a[ijk] + scale_b*grid_b[ijk];
    ijk+=plane;
  } // for k
} // add_grids kernel


//==============================================================================================================================================================
__global__ void __rebuild_lambda(subdomain_type * gpu_subdomains, double a, double b, double h, int level){ 
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  double h2inv = 1.0/(h*h);
  double * alpha  = gpu_subdomains[box].levels[level].grids[__alpha ] + (1+pencil+plane);
  double * beta_i = gpu_subdomains[box].levels[level].grids[__beta_i] + (1+pencil+plane);
  double * beta_j = gpu_subdomains[box].levels[level].grids[__beta_j] + (1+pencil+plane);
  double * beta_k = gpu_subdomains[box].levels[level].grids[__beta_k] + (1+pencil+plane);
  double * lambda = gpu_subdomains[box].levels[level].grids[__lambda] + (1+pencil+plane);
  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  k = 0; 
  int ijk = k*plane + j*pencil + i;

  for(k=0;k<subdomain_dim;k++){
    if(withinBounds){
      // centr of Gershgorin disc is the diagonal element...
      double    Aii = a*alpha[ijk] - b*h2inv*( -beta_i[ijk]-beta_i[ijk+     1]
                                               -beta_j[ijk]-beta_j[ijk+pencil]
                                               -beta_k[ijk]-beta_k[ijk+ plane] );
      lambda[ijk] = 1.0/Aii; // inverse of the diagonal Aii
    }
    ijk+=plane;
  } // for k
} // rebuild_lambda kernel


//==============================================================================================================================================================
__global__ void __initialize_exact(subdomain_type * gpu_subdomains, int level, double hLevel, double a, double b){
  double NPi = 2.0*M_PI;
  double Bmin =  1.0;
  double Bmax = 10.0;
  double c2 = (Bmax-Bmin)/2;
  double c1 = (Bmax+Bmin)/2;
  double c3=10.0; // how sharply (B)eta transitions
  double c4 = -5.0/0.25;

  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int   subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int          pencil = gpu_subdomains[box].levels[level].pencil;
  int           plane = gpu_subdomains[box].levels[level].plane;
  int           low_i = gpu_subdomains[box].levels[level].low.i;
  int           low_j = gpu_subdomains[box].levels[level].low.j;
  int           low_k = gpu_subdomains[box].levels[level].low.k;
  double * grid_u     = gpu_subdomains[box].levels[level].grids[__u_exact] + (plane+pencil+1);
  double * grid_f     = gpu_subdomains[box].levels[level].grids[__f      ] + (plane+pencil+1);
  double * grid_alpha = gpu_subdomains[box].levels[level].grids[__alpha  ] + (plane+pencil+1);
  double * grid_beta  = gpu_subdomains[box].levels[level].grids[__beta   ] + (plane+pencil+1);
  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  k = 0;
  int ijk = k*plane + j*pencil + i;

  for(k=0;k<subdomain_dim;k++){
      double x = hLevel*((double)(i+low_i)+0.5);
      double y = hLevel*((double)(j+low_j)+0.5);
      double z = hLevel*((double)(k+low_k)+0.5);
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      double r2   = pow((x-0.50),2) +  pow((y-0.50),2) +  pow((z-0.50),2); // distance from center squared
      double r2x  = 2.0*(x-0.50);
      double r2y  = 2.0*(y-0.50);
      double r2z  = 2.0*(z-0.50);
      double r2xx = 2.0;
      double r2yy = 2.0;
      double r2zz = 2.0;
      double r    = pow(r2,0.5);
      double rx   = 0.5*r2x*pow(r2,-0.5);
      double ry   = 0.5*r2y*pow(r2,-0.5);
      double rz   = 0.5*r2z*pow(r2,-0.5);
      double rxx  = 0.5*r2xx*pow(r2,-0.5) - 0.25*r2x*r2x*pow(r2,-1.5);
      double ryy  = 0.5*r2yy*pow(r2,-0.5) - 0.25*r2y*r2y*pow(r2,-1.5);
      double rzz  = 0.5*r2zz*pow(r2,-0.5) - 0.25*r2z*r2z*pow(r2,-1.5);
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      double A  = 1.0;
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      #if 1
      double B  =           c1+c2*tanh( c3*(r-0.25) );
      double Bx = c2*c3*rx*(1-pow(tanh( c3*(r-0.25) ),2));
      double By = c2*c3*ry*(1-pow(tanh( c3*(r-0.25) ),2));
      double Bz = c2*c3*rz*(1-pow(tanh( c3*(r-0.25) ),2));
      #else
      double B  = 1.0;
      double Bx = 0.0;
      double By = 0.0;
      double Bz = 0.0;
      #endif
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      #if 1
      double u   =                exp(c4*r2)*sin(NPi*x)*sin(NPi*y)*sin(NPi*z);
      double ux  = c4*r2x*u + NPi*exp(c4*r2)*cos(NPi*x)*sin(NPi*y)*sin(NPi*z);
      double uy  = c4*r2y*u + NPi*exp(c4*r2)*sin(NPi*x)*cos(NPi*y)*sin(NPi*z);
      double uz  = c4*r2z*u + NPi*exp(c4*r2)*sin(NPi*x)*sin(NPi*y)*cos(NPi*z);
      double uxx = c4*r2xx*u + c4*r2x*ux + c4*r2x*NPi*exp(c4*r2)*cos(NPi*x)*sin(NPi*y)*sin(NPi*z) - NPi*NPi*exp(c4*r2)*sin(NPi*x)*sin(NPi*y)*sin(NPi*z);
      double uyy = c4*r2yy*u + c4*r2y*uy + c4*r2y*NPi*exp(c4*r2)*sin(NPi*x)*cos(NPi*y)*sin(NPi*z) - NPi*NPi*exp(c4*r2)*sin(NPi*x)*sin(NPi*y)*sin(NPi*z);
      double uzz = c4*r2zz*u + c4*r2z*uz + c4*r2z*NPi*exp(c4*r2)*sin(NPi*x)*sin(NPi*y)*cos(NPi*z) - NPi*NPi*exp(c4*r2)*sin(NPi*x)*sin(NPi*y)*sin(NPi*z);
      double f = a*A*u - b*( (Bx*ux + By*uy + Bz*uz)  +  B*(uxx + uyy + uzz) );
      #else
      // should be continuous in u, u', and u''
      // v(w) = w^4 - 2w^3 + w^2
      // u(x,y,z) = v(x)v(y)v(z)
      double X   =  1.0*pow(x,4) -  2.0*pow(x,3) + 1.0*pow(x,2);
      double Y   =  1.0*pow(y,4) -  2.0*pow(y,3) + 1.0*pow(y,2);
      double Z   =  1.0*pow(z,4) -  2.0*pow(z,3) + 1.0*pow(z,2);
      double Xx  =  4.0*pow(x,3) -  6.0*pow(x,2) + 2.0*x;
      double Yy  =  4.0*pow(y,3) -  6.0*pow(y,2) + 2.0*y;
      double Zz  =  4.0*pow(z,3) -  6.0*pow(z,2) + 2.0*z;
      double Xxx = 12.0*pow(x,2) - 12.0*x        + 2.0;
      double Yyy = 12.0*pow(y,2) - 12.0*y        + 2.0;
      double Zzz = 12.0*pow(z,2) - 12.0*z        + 2.0;
      double u   = X*Y*Z;
      double ux  = Xx*Y*Z;
      double uy  = X*Yy*Z;
      double uz  = X*Y*Zz;
      double uxx = Xxx*Y*Z;
      double uyy = X*Yyy*Z;
      double uzz = X*Y*Zzz;
      double f = a*A*u - b*( (Bx*ux + By*uy + Bz*uz)  +  B*(uxx + uyy + uzz) );
      #endif
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      if(withinBounds){
        grid_alpha[ijk] = A;
        grid_beta[ijk]  = B;
        grid_u[ijk]     = u;
        grid_f[ijk]     = f;
      }
      ijk+=plane;
  } // for k
} // initialize exact on gpu kernel


//=============================================================================================================================================================
__global__ void __project_cell_to_face(subdomain_type * gpu_subdomains, int level, int id_cell, int id_face, int dir){
  int box = blockIdx.z;  // CUDA 4 !!!
  int i,j,k; // (0,0,0) = first non ghost zone element
  i = TBDIMX*blockIdx.x + threadIdx.x;
  j = TBDIMY*blockIdx.y + threadIdx.y;

  int subdomain_dim = gpu_subdomains[box].levels[level].dim.i;
  int pencil = gpu_subdomains[box].levels[level].pencil;
  int  plane = gpu_subdomains[box].levels[level].plane;

  double * grid_cell = gpu_subdomains[box].levels[level].grids[id_cell] + (1+pencil+plane);
  double * grid_face = gpu_subdomains[box].levels[level].grids[id_face] + (1+pencil+plane);
  // when the subdomain is coarsened to smaller than a thread block, certain threads won't commit their result -------------------------------------------------
                  int withinBounds = 1;
  if(i>=subdomain_dim)withinBounds = 0;
  if(j>=subdomain_dim)withinBounds = 0;

  int stride;
  switch(dir){
    case 0: stride =      1;break;//i-direction
    case 1: stride = pencil;break;//j-direction
    case 2: stride =  plane;break;//k-direction
  }

  k = 0;
  int ijk = k*plane + j*pencil + i;

  for(k=0;k<subdomain_dim;k++){
    if(withinBounds){
      grid_face[ijk] = 0.5*(grid_cell[ijk-stride] + grid_cell[ijk]); // simple linear interpolation
    }
    ijk+=plane;
  } // for k
} // projection of cell-centered coefficients to face-centered 


//==============================================================================================================================================================
// wrappers...
//==============================================================================================================================================================
extern "C" void zero_grid(domain_type *domain, int level, int grid_id){
  // zeros the grid and the ghost zones
  enqueueEvent(cudaEvent_blas1);
  int box;
  for(box=0;box<domain->numsubdomains;box++){
    cudaMemsetAsync(domain->subdomains[box].levels[level].grids[grid_id], 0, domain->subdomains[0].levels[level].volume*sizeof(double),0);
  }
  enqueueEvent(cudaEvent_blas1);
}

extern "C" void initialize_grid_to_scalar(domain_type *domain, int level, int grid_id, double h, double scalar){
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __initialize_grid_to_scalar<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,grid_id,h,scalar,level);
  enqueueEvent(cudaEvent_blas1);
}


extern "C" void scale_grid(domain_type *domain, int level, int id_c, double scale_a, int id_a){
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __scale_grid<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,id_c,scale_a,id_a,level);
  enqueueEvent(cudaEvent_blas1);
}


extern "C" void mul_grids(domain_type *domain, int level, int id_c, double scale, int id_a, int id_b){
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __mul_grids<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,id_c,scale,id_a,id_b,level);
  enqueueEvent(cudaEvent_blas1);
}


extern "C" void shift_grid(domain_type *domain, int level, int id_c, int id_a, double scale_a){
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __shift_grid<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,id_c,id_a,scale_a,level);
  enqueueEvent(cudaEvent_blas1);
}

extern "C" void project_cell_to_face(domain_type *domain, int level, int id_cell, int id_face, int dir){
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __project_cell_to_face<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,level,id_cell,id_face,dir);
  enqueueEvent(cudaEvent_blas1);
}


extern "C" void add_grids(domain_type *domain, int level, int id_c, double scale_a, int id_a, double scale_b, int id_b){
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __add_grids<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,id_c,scale_a,id_a,scale_b,id_b,level);
  enqueueEvent(cudaEvent_blas1);
}


extern "C" void rebuild_lambda(domain_type *domain, int level, double a, double b, double hLevel){
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __rebuild_lambda<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,a,b,hLevel,level);
  enqueueEvent(cudaEvent_blas1);
}


extern "C" void grid_to_surface_buffers(domain_type *domain, int level, int grid_id){
  enqueueEvent(cudaEvent_s2buf);
  // each thread block copies all minimum stride points for all faces for a given least unit stride offset of one box
  dim3 dimBlock(domain->subdomains[0].levels[level].dim.i,1,1);                    // a stanza of dim points
  dim3 dimGrid(1,domain->subdomains[0].levels[level].dim.i,domain->numsubdomains); // all offsets x all boxes
  __grid_to_surface_buffers<<<dimGrid, dimBlock>>>(domain->gpu_subdomains, grid_id,level);
  enqueueEvent(cudaEvent_s2buf);
}


extern "C" void surface_buffers_to_ghost_buffers(domain_type *domain, int level, int grid_id){
  enqueueEvent(cudaEvent_bufcopy);
  // each thread block copies dim.i points for all faces for a given offset within a face of one box
  dim3 dimBlock(domain->subdomains[0].levels[level].dim.i,1,1);                    // a stanza of dim points
  dim3 dimGrid(1,domain->subdomains[0].levels[level].dim.i,domain->numsubdomains); // all offsets x all boxes
  __surface_buffers_to_ghost_buffers<<<dimGrid, dimBlock>>>(domain->gpu_subdomains, grid_id,level,domain->rank);
  enqueueEvent(cudaEvent_bufcopy);
}

#ifdef __MPI
extern "C" void surface_buffers_to_send_buffer(domain_type *domain, int level, int grid_id){
  enqueueEvent(cudaEvent_pack);
  // each thread block copies dim.i points for all faces for a given offset within a face of one box
  dim3 dimBlock(domain->subdomains[0].levels[level].dim.i,1,1);                    // a stanza of dim points
  dim3 dimGrid(1,domain->subdomains[0].levels[level].dim.i,domain->numsubdomains); // all offsets x all boxes
  __surface_buffers_to_send_buffer<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,domain->gpu_pointers_to_gpu_send_buffer,domain->gpu_pointers_to_gpu_recv_buffer,grid_id,level,domain->rank);
  enqueueEvent(cudaEvent_pack);
}
extern "C" void recv_buffer_to_ghost_buffers(domain_type *domain, int level, int grid_id){
  enqueueEvent(cudaEvent_unpack);
  // each thread block copies dim.i points for all faces for a given offset within a face of one box
  dim3 dimBlock(domain->subdomains[0].levels[level].dim.i,1,1);                    // a stanza of dim points
  dim3 dimGrid(1,domain->subdomains[0].levels[level].dim.i,domain->numsubdomains); // all offsets x all boxes
  __recv_buffer_to_ghost_buffers<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,domain->gpu_pointers_to_gpu_send_buffer,domain->gpu_pointers_to_gpu_recv_buffer,grid_id,level,domain->rank);
  enqueueEvent(cudaEvent_unpack);
}
#endif

extern "C" void ghost_buffers_to_grid(domain_type *domain, int level, int grid_id){
  enqueueEvent(cudaEvent_buf2g);
  // each thread block copies all minimum stride points for all faces for a given least unit stride offset of one box
  dim3 dimBlock(domain->subdomains[0].levels[level].dim.i,1,1);                    // a stanza of dim points
  dim3 dimGrid(1,domain->subdomains[0].levels[level].dim.i,domain->numsubdomains); // all offsets x all boxes
  __ghost_buffers_to_grid<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,grid_id,level);
  enqueueEvent(cudaEvent_buf2g);
}

extern "C" void smooth(domain_type *domain, int level, int phi_id, int rhs_id, double a, double b, double hLevel, int s){
  enqueueEvent(cudaEvent_smooth);
  #ifdef VL
  int IJStride = (VL - 2*(domain->subdomains[0].levels[level].dim.i+2)) & ~0x0F; // i.e. VL minus halos rounded down to a multiple of 16
  if(IJStride==0){printf("Error, IJStride==0 for VL=%d, dim=%d\n",VL,domain->subdomains[0].levels[level].dim.i);exit(0);}
  int NumVectors = ( (domain->subdomains[0].levels[level].dim.i+2)*(domain->subdomains[0].levels[level].dim.j)+IJStride+16-1)/IJStride; // i.e. size of a plane less the first/last halo + 16 for rounding
  // FIX, round VLSat up to a multiple of 16
  int VLSat = VL;
  if(VLSat>(domain->subdomains[0].levels[level].dim.i+2)*(domain->subdomains[0].levels[level].dim.j+2)+16)
     VLSat=(domain->subdomains[0].levels[level].dim.i+2)*(domain->subdomains[0].levels[level].dim.j+2)+16;
  dim3 dimBlock(VLSat,1);
  dim3 dimGrid(NumVectors,1,domain->numsubdomains);
  #else
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1 
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  #endif
  __smooth_once_GSRB<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,phi_id,rhs_id,a,b,hLevel,s,level);
  enqueueEvent(cudaEvent_smooth);
}

extern "C" void restriction(domain_type *domain, int level, int coarse_id, int fine_id){
  enqueueEvent(cudaEvent_restriction);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __restriction<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,fine_id,coarse_id,level); // fine_id@level -> coarse_id@level+1
  enqueueEvent(cudaEvent_restriction);
}

extern "C" void restriction_betas(domain_type *domain, int level){
  enqueueEvent(cudaEvent_restriction);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __restriction_betas<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,level); // update betas at level+1
  enqueueEvent(cudaEvent_restriction);
}

extern "C" void residual(domain_type * domain, int level, int res_id, int phi_id, int rhs_id, double a, double b, double hLevel){
  enqueueEvent(cudaEvent_residual);
  #ifdef VL   
  int IJStride = (VL - 2*(domain->subdomains[0].levels[level].dim.i+2)) & ~0x0F; // i.e. VL minus halos rounded down to a multiple of 16
  if(IJStride==0){printf("Error, IJStride==0 for VL=%d, dim=%d\n",VL,domain->subdomains[0].levels[level].dim.i);exit(0);}
  int NumVectors = ( (domain->subdomains[0].levels[level].dim.i+2)*(domain->subdomains[0].levels[level].dim.j)+IJStride+16-1)/IJStride; // i.e. size of a plane less the first/last halo + 16 for rounding
  // FIX, round VLSat up to a multiple of 16
  int VLSat = VL;
  if(VLSat>(domain->subdomains[0].levels[level].dim.i+2)*(domain->subdomains[0].levels[level].dim.j+2)+16)
     VLSat=(domain->subdomains[0].levels[level].dim.i+2)*(domain->subdomains[0].levels[level].dim.j+2)+16;
  dim3 dimBlock(VLSat,1);
  dim3 dimGrid(NumVectors,1,domain->numsubdomains);
  #else
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  #endif
  __residual<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,res_id,phi_id,rhs_id,a,b,hLevel,level);
  enqueueEvent(cudaEvent_residual);
}

extern "C" void apply_op(domain_type * domain, int level, int     Ax_id, int      x_id, double a, double b, double hLevel){
  enqueueEvent(cudaEvent_apply_op);
  #ifdef VL   
  int IJStride = (VL - 2*(domain->subdomains[0].levels[level].dim.i+2)) & ~0x0F; // i.e. VL minus halos rounded down to a multiple of 16
  if(IJStride==0){printf("Error, IJStride==0 for VL=%d, dim=%d\n",VL,domain->subdomains[0].levels[level].dim.i);exit(0);}
  int NumVectors = ( (domain->subdomains[0].levels[level].dim.i+2)*(domain->subdomains[0].levels[level].dim.j)+IJStride+16-1)/IJStride; // i.e. size of a plane less the first/last halo + 16 for rounding
  // FIX, round VLSat up to a multiple of 16
  int VLSat = VL;
  if(VLSat>(domain->subdomains[0].levels[level].dim.i+2)*(domain->subdomains[0].levels[level].dim.j+2)+16)
     VLSat=(domain->subdomains[0].levels[level].dim.i+2)*(domain->subdomains[0].levels[level].dim.j+2)+16;
  dim3 dimBlock(VLSat,1);
  dim3 dimGrid(NumVectors,1,domain->numsubdomains);
  #else
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  #endif
  __apply_op<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,Ax_id,x_id,a,b,hLevel,level);
  enqueueEvent(cudaEvent_apply_op);
}


extern "C" void interpolation(domain_type * domain, int level, int fine_id, int coarse_id){ // interpolate from level+1 to level
  enqueueEvent(cudaEvent_interpolation);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __interpolation<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,coarse_id,fine_id,level); // interpolate from level+1 onto level
  enqueueEvent(cudaEvent_interpolation);
}


extern "C" void norm_on_gpu(domain_type * domain, int level, int grid_id, double *gpu_norm){ // norm returned in *gpu_norm.  However, must be cudaMemcpy'd to host
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1 
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __norm<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,grid_id,gpu_norm,level);
  enqueueEvent(cudaEvent_blas1);
}

extern "C" void dot_on_gpu(domain_type * domain, int level, int id_a, int id_b, double *gpu_dot){ // dot returned in *gpu_dot.  However, must be cudaMemcpy'd to host
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1 
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __dot<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,id_a,id_b,gpu_dot,level);
  enqueueEvent(cudaEvent_blas1);
}

extern "C" void sum_on_gpu(domain_type * domain, int level, int id_a, double *gpu_sum){ // sum returned in *gpu_sum.  However, must be cudaMemcpy'd to host
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1 
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __sum<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,id_a,gpu_sum,level);
  enqueueEvent(cudaEvent_blas1);
}

//=======================================================================================================================================================================
extern "C" void initialize_exact_on_gpu(domain_type * domain, int level, double hLevel, double a, double b){
  enqueueEvent(cudaEvent_blas1);
  int GDIMX = (domain->subdomains[0].levels[level].dim.i+TBDIMX-1)/TBDIMX; // i.e. what happens when TBDIMX==32 on 4^3, GDIMX should always be >=1 
  int GDIMY = (domain->subdomains[0].levels[level].dim.j+TBDIMY-1)/TBDIMY; // i.e. what happens when TBDIMY== 8 on 4^3, GDIMY should always be >=1
  dim3 dimBlock(TBDIMX,TBDIMY);
  dim3 dimGrid(GDIMX,GDIMY,domain->numsubdomains);
  __initialize_exact<<<dimGrid, dimBlock>>>(domain->gpu_subdomains,level,hLevel,a,b);
  enqueueEvent(cudaEvent_blas1);
}
//=======================================================================================================================================================================
extern "C" void ConfigureGPU(){
  cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);  // K20x shared memory banking optimized for 64-bit accesses
#ifdef __PREFER_SHARED
  cudaFuncSetCacheConfig(__smooth_once_GSRB, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(        __residual, cudaFuncCachePreferShared);
  cudaFuncSetCacheConfig(        __apply_op, cudaFuncCachePreferShared);
#else
  cudaFuncSetCacheConfig(__smooth_once_GSRB, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(        __residual, cudaFuncCachePreferL1);
  cudaFuncSetCacheConfig(        __apply_op, cudaFuncCachePreferL1);
#endif
}
