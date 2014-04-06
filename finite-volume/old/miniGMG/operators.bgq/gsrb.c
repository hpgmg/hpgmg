//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include "../timer.h"
#include <spi/include/l1p/sprefetch.h>
//------------------------------------------------------------------------------------------------------------------------------
void __box_smooth_GSRB_multiple(box_type *box, int phi_id, int rhs_id, double a, double b, int sweep){
    int pencil = box->pencil;
    int plane = box->plane;
    int plane8   = plane<<3;
    int pencil8  = pencil<<3;
    int ghosts = box->ghosts;
    int DimI = box->dim.i;
    int DimJ = box->dim.j;
    int DimK = box->dim.k;
    double h2inv = 1.0/(box->h*box->h);
    double  * __restrict__ phi          = box->grids[  phi_id] + ghosts*plane;
    double  * __restrict__ rhs          = box->grids[  rhs_id] + ghosts*plane;
    double  * __restrict__ alpha        = box->grids[__alpha ] + ghosts*plane;
    double  * __restrict__ beta_i       = box->grids[__beta_i] + ghosts*plane;
    double  * __restrict__ beta_j       = box->grids[__beta_j] + ghosts*plane;
    double  * __restrict__ beta_k       = box->grids[__beta_k] + ghosts*plane;
    double  * __restrict__ lambda       = box->grids[__lambda] + ghosts*plane;
    uint64_t* __restrict__ RedBlackMask = box->RedBlack_64bMask;
        vector4double       a_splat4 = vec_splats(a);
        vector4double b_h2inv_splat4 = vec_mul(vec_splats(b),vec_splats(h2inv));
    int smt_id      = 0;
    int smt_leader  = 0;
    int global_ij_start[8];
    int global_ij_end[8];
    int ij_start[8];
    int ij_end[8];
    int planeInWavefront;for(planeInWavefront=0;planeInWavefront<ghosts;planeInWavefront++){
      global_ij_start[planeInWavefront] = (                   (1+planeInWavefront)*pencil)&~3;
      global_ij_end[planeInWavefront]   = ((ghosts+DimJ+ghosts-1-planeInWavefront)*pencil);
      ij_start[planeInWavefront] = global_ij_start[planeInWavefront];
      ij_end[planeInWavefront]   = global_ij_end[planeInWavefront];
    }
    int leadingK;
    int kLow  =     -(ghosts-1);
    int kHigh = DimK+(ghosts-1);
    for(leadingK=kLow-1;leadingK<kHigh;leadingK++){
      {
        int j,k;
        for(planeInWavefront=0;planeInWavefront<ghosts;planeInWavefront++){
          k=(leadingK-planeInWavefront);if((k>=kLow)&&(k<kHigh)){
        uint64_t invertMask = 0-((1^sweep^k^planeInWavefront)&0x1);
        vector4double invertMask4 = vec_splats(*((double*)&invertMask));
        int ij           = ij_start[planeInWavefront];
        int ijk_start    = ij_start[planeInWavefront]+k*plane;
        int ijk          = ijk_start;
        int _ij_end      = ij_end[planeInWavefront];
        while(ij<_ij_end){ // smooth an interleaved vector...
                                                         __dcbtst((void*)(   phi+ijk+4));
                vector4double                 phi_m4   = vec_lda(-32,(double*)(   phi+ijk));
                vector4double                 phi_00 = vec_lda(  0,(double*)(   phi+ijk));
                vector4double                 phi_04 = vec_lda( 32,(double*)(   phi+ijk));
                vector4double              beta_i_00 = vec_lda(  0,(double*)(beta_i+ijk));
                vector4double              beta_i_04 = vec_lda( 32,(double*)(beta_i+ijk));
                vector4double                 phi_m1   = vec_sldw(   phi_m4  ,   phi_00,3);
                vector4double                 phi_01 = vec_sldw(   phi_00,   phi_04,1);
                vector4double              beta_i_01 = vec_sldw(beta_i_00,beta_i_04,1);
                vector4double        phi_jm1_temp_00 = vec_ld(   0-pencil8,(double*)(   phi+ijk));
                vector4double        phi_jm1_temp_04 = vec_ld(  32-pencil8,(double*)(   phi+ijk));
                vector4double        phi_jp1_temp_00 = vec_ld(   0+pencil8,(double*)(   phi+ijk));
                vector4double        phi_jp1_temp_04 = vec_ld(  32+pencil8,(double*)(   phi+ijk));
                vector4double     beta_j_jp1_temp_00 = vec_ld(   0+pencil8,(double*)(beta_j+ijk));
                vector4double     beta_j_jp1_temp_04 = vec_ld(  32+pencil8,(double*)(beta_j+ijk));
                vector4double              beta_j_00 = vec_lda(  0        ,(double*)(beta_j+ijk));
                vector4double          phi_jm1_permute = vec_lvsl( 0-pencil8,(double*)(   phi+ijk));
                vector4double          phi_jp1_permute = vec_lvsl( 0+pencil8,(double*)(   phi+ijk));
                vector4double       beta_j_jp1_permute = vec_lvsl( 0+pencil8,(double*)(beta_j+ijk));
                vector4double             phi_jm1_00 = vec_perm(   phi_jm1_temp_00,   phi_jm1_temp_04,   phi_jm1_permute);
                vector4double             phi_jp1_00 = vec_perm(   phi_jp1_temp_00,   phi_jp1_temp_04,   phi_jp1_permute);
                vector4double          beta_j_jp1_00 = vec_perm(beta_j_jp1_temp_00,beta_j_jp1_temp_04,beta_j_jp1_permute);
                vector4double              beta_k_00 = vec_lda(  0       ,(double*)(beta_k+ijk));
                vector4double          beta_k_kp1_00 = vec_lda(  0+plane8,(double*)(beta_k+ijk));
                vector4double             phi_km1_00 = vec_lda(  0-plane8,(double*)(   phi+ijk));
                vector4double             phi_kp1_00 = vec_lda(  0+plane8,(double*)(   phi+ijk));
                vector4double           laplacian_00;
                vector4double         laplacian_i_00;
                vector4double         laplacian_k_00;
                vector4double         laplacian_j_00;
                                      laplacian_i_00 = vec_mul(  beta_i_00    ,phi_m1                       );
                                      laplacian_i_00 = vec_nmsub(beta_i_00    ,phi_00    ,laplacian_i_00);
                                      laplacian_i_00 = vec_nmsub(beta_i_01    ,phi_00    ,laplacian_i_00);
                                      laplacian_i_00 = vec_madd( beta_i_01    ,phi_01    ,laplacian_i_00);
                                      laplacian_j_00 = vec_mul(  beta_j_00    ,phi_jm1_00                 );
                                      laplacian_j_00 = vec_nmsub(beta_j_00    ,phi_00    ,laplacian_j_00);
                                      laplacian_j_00 = vec_nmsub(beta_j_jp1_00,phi_00    ,laplacian_j_00);
                                      laplacian_j_00 = vec_madd( beta_j_jp1_00,phi_jp1_00,laplacian_j_00);
                                      laplacian_k_00 = vec_mul(  beta_k_00    ,phi_km1_00                 );
                                      laplacian_k_00 = vec_nmsub(beta_k_00    ,phi_00    ,laplacian_k_00);
                                      laplacian_k_00 = vec_nmsub(beta_k_kp1_00,phi_00    ,laplacian_k_00);
                                      laplacian_k_00 = vec_madd( beta_k_kp1_00,phi_kp1_00,laplacian_k_00);
                vector4double               alpha_00 = vec_lda(  0,(double*)( alpha+ijk));
                vector4double         a_alpha_phi_00 = vec_mul(vec_mul(a_splat4,phi_00),alpha_00);
                vector4double           helmholtz_00 = vec_nmsub(b_h2inv_splat4,vec_add(laplacian_i_00,laplacian_j_00),vec_nmsub(b_h2inv_splat4,laplacian_k_00,a_alpha_phi_00));
                vector4double                 rhs_00 = vec_lda(  0,(double*)(   rhs+ijk));
                vector4double              lambda_00 = vec_lda(  0,(double*)(lambda+ijk));
                vector4double            RedBlack_00 = vec_xor(invertMask4,vec_lda( 0,(double*)(RedBlackMask+ij) ));
                vector4double phi_plus_lambda_rhs_00 = vec_madd(lambda_00,rhs_00,phi_00);
                vector4double                 new_00 = vec_nmsub(lambda_00,helmholtz_00,phi_plus_lambda_rhs_00);
                                              new_00 = vec_sel(new_00,phi_00,RedBlack_00);
                                                         vec_sta(new_00, 0,(double*)(phi+ijk));
          ij+=4;ijk+=4;
        }
        }}
      } // if stencil
    } // leadingK
}


void __box_smooth_GSRB_multiple_threaded(box_type *box, int phi_id, int rhs_id, double a, double b, int sweep){
  #pragma omp parallel
  {
    int pencil = box->pencil;
    int plane = box->plane;
    int plane8   = plane<<3;
    int pencil8  = pencil<<3;
    int ghosts = box->ghosts;
    int DimI = box->dim.i;
    int DimJ = box->dim.j;
    int DimK = box->dim.k;
    double h2inv = 1.0/(box->h*box->h);
    double  * __restrict__ phi          = box->grids[  phi_id] + ghosts*plane;
    double  * __restrict__ rhs          = box->grids[  rhs_id] + ghosts*plane;
    double  * __restrict__ alpha        = box->grids[__alpha ] + ghosts*plane;
    double  * __restrict__ beta_i       = box->grids[__beta_i] + ghosts*plane;
    double  * __restrict__ beta_j       = box->grids[__beta_j] + ghosts*plane;
    double  * __restrict__ beta_k       = box->grids[__beta_k] + ghosts*plane;
    double  * __restrict__ lambda       = box->grids[__lambda] + ghosts*plane;
    uint64_t* __restrict__ RedBlackMask = box->RedBlack_64bMask;
        vector4double       a_splat4 = vec_splats(a);
        vector4double b_h2inv_splat4 = vec_mul(vec_splats(b),vec_splats(h2inv));
    int id      = omp_get_thread_num();
    int threads = omp_get_num_threads();
    int DRAM_scout_threads = 0;
    int smt_id      = id &  3;
    int smt_leader  = id & ~3;
    int global_ij_start[4];
    int global_ij_end[4];
    int ij_start[4];
    int ij_end[4];
    int planeInWavefront;for(planeInWavefront=0;planeInWavefront<ghosts;planeInWavefront++){
      global_ij_start[planeInWavefront] = (                   (1)*pencil)&~3;
      global_ij_end[planeInWavefront]   = ((ghosts+DimJ+ghosts-1)*pencil);
      int TotalUnrollings = ((global_ij_end[planeInWavefront]-global_ij_start[planeInWavefront]+4-1)/4);
      ij_start[planeInWavefront] = global_ij_start[planeInWavefront] + 4*( (smt_leader  )*(TotalUnrollings)/(threads-DRAM_scout_threads)) + 4*smt_id;
      ij_end[planeInWavefront]   = global_ij_start[planeInWavefront] + 4*( (smt_leader+4)*(TotalUnrollings)/(threads-DRAM_scout_threads));
      if(id>=(threads-1-DRAM_scout_threads))ij_end[planeInWavefront] = global_ij_end[planeInWavefront];
      if(ij_end[planeInWavefront]>global_ij_end[planeInWavefront])ij_end[planeInWavefront]=global_ij_end[planeInWavefront];
    }
    uint32_t old_stream_depth;L1P_GetStreamDepth(&old_stream_depth);
    if(smt_id==0){
      L1P_SetStreamPolicy(L1P_stream_confirmed);
      if(id < (threads-DRAM_scout_threads))L1P_SetStreamDepth(2);
      if(id >=(threads-DRAM_scout_threads))L1P_SetStreamDepth(6);
    }
    int leadingK;
    int kLow  =     -(ghosts-1);
    int kHigh = DimK+(ghosts-1);
    for(leadingK=kLow-1;leadingK<kHigh;leadingK++){
      if(ghosts>1){
        #pragma omp barrier
      }
      if(smt_id<3){
        int j,k;
        for(planeInWavefront=0;planeInWavefront<ghosts;planeInWavefront++){
          k=(leadingK-planeInWavefront);if((k>=kLow)&&(k<kHigh)){
        uint64_t invertMask = 0-((1^sweep^k^planeInWavefront)&0x1);
        vector4double invertMask4 = vec_splats(*((double*)&invertMask));
        int ij           = ij_start[planeInWavefront];
        int ijk_start    = ij_start[planeInWavefront]+k*plane;
        int ijk          = ijk_start;
        int _ij_end      = ij_end[planeInWavefront];
        while(ij<_ij_end){ // smooth an interleaved vector...
                                                         __dcbtst((void*)(   phi+ijk+12));
                vector4double                 phi_m4   = vec_lda(-32,(double*)(   phi+ijk));
                vector4double                 phi_00 = vec_lda(  0,(double*)(   phi+ijk));
                vector4double                 phi_04 = vec_lda( 32,(double*)(   phi+ijk));
                vector4double              beta_i_00 = vec_lda(  0,(double*)(beta_i+ijk));
                vector4double              beta_i_04 = vec_lda( 32,(double*)(beta_i+ijk));
                vector4double                 phi_m1   = vec_sldw(   phi_m4  ,   phi_00,3);
                vector4double                 phi_01 = vec_sldw(   phi_00,   phi_04,1);
                vector4double              beta_i_01 = vec_sldw(beta_i_00,beta_i_04,1);
                vector4double        phi_jm1_temp_00 = vec_ld(   0-pencil8,(double*)(   phi+ijk));
                vector4double        phi_jm1_temp_04 = vec_ld(  32-pencil8,(double*)(   phi+ijk));
                vector4double        phi_jp1_temp_00 = vec_ld(   0+pencil8,(double*)(   phi+ijk));
                vector4double        phi_jp1_temp_04 = vec_ld(  32+pencil8,(double*)(   phi+ijk));
                vector4double     beta_j_jp1_temp_00 = vec_ld(   0+pencil8,(double*)(beta_j+ijk));
                vector4double     beta_j_jp1_temp_04 = vec_ld(  32+pencil8,(double*)(beta_j+ijk));
                vector4double              beta_j_00 = vec_lda(  0        ,(double*)(beta_j+ijk));
                vector4double          phi_jm1_permute = vec_lvsl( 0-pencil8,(double*)(   phi+ijk));
                vector4double          phi_jp1_permute = vec_lvsl( 0+pencil8,(double*)(   phi+ijk));
                vector4double       beta_j_jp1_permute = vec_lvsl( 0+pencil8,(double*)(beta_j+ijk));
                vector4double             phi_jm1_00 = vec_perm(   phi_jm1_temp_00,   phi_jm1_temp_04,   phi_jm1_permute);
                vector4double             phi_jp1_00 = vec_perm(   phi_jp1_temp_00,   phi_jp1_temp_04,   phi_jp1_permute);
                vector4double          beta_j_jp1_00 = vec_perm(beta_j_jp1_temp_00,beta_j_jp1_temp_04,beta_j_jp1_permute);
                vector4double              beta_k_00 = vec_lda(  0       ,(double*)(beta_k+ijk));
                vector4double          beta_k_kp1_00 = vec_lda(  0+plane8,(double*)(beta_k+ijk));
                vector4double             phi_km1_00 = vec_lda(  0-plane8,(double*)(   phi+ijk));
                vector4double             phi_kp1_00 = vec_lda(  0+plane8,(double*)(   phi+ijk));
                vector4double           laplacian_00;
                vector4double         laplacian_i_00;
                vector4double         laplacian_k_00;
                vector4double         laplacian_j_00;
                                      laplacian_i_00 = vec_mul(  beta_i_00    ,phi_m1                       );
                                      laplacian_i_00 = vec_nmsub(beta_i_00    ,phi_00    ,laplacian_i_00);
                                      laplacian_i_00 = vec_nmsub(beta_i_01    ,phi_00    ,laplacian_i_00);
                                      laplacian_i_00 = vec_madd( beta_i_01    ,phi_01    ,laplacian_i_00);
                                      laplacian_j_00 = vec_mul(  beta_j_00    ,phi_jm1_00                 );
                                      laplacian_j_00 = vec_nmsub(beta_j_00    ,phi_00    ,laplacian_j_00);
                                      laplacian_j_00 = vec_nmsub(beta_j_jp1_00,phi_00    ,laplacian_j_00);
                                      laplacian_j_00 = vec_madd( beta_j_jp1_00,phi_jp1_00,laplacian_j_00);
                                      laplacian_k_00 = vec_mul(  beta_k_00    ,phi_km1_00                 );
                                      laplacian_k_00 = vec_nmsub(beta_k_00    ,phi_00    ,laplacian_k_00);
                                      laplacian_k_00 = vec_nmsub(beta_k_kp1_00,phi_00    ,laplacian_k_00);
                                      laplacian_k_00 = vec_madd( beta_k_kp1_00,phi_kp1_00,laplacian_k_00);
                vector4double               alpha_00 = vec_lda(  0,(double*)( alpha+ijk));
                vector4double         a_alpha_phi_00 = vec_mul(vec_mul(a_splat4,phi_00),alpha_00);
                vector4double           helmholtz_00 = vec_nmsub(b_h2inv_splat4,vec_add(laplacian_i_00,laplacian_j_00),vec_nmsub(b_h2inv_splat4,laplacian_k_00,a_alpha_phi_00));
                vector4double                 rhs_00 = vec_lda(  0,(double*)(   rhs+ijk));
                vector4double              lambda_00 = vec_lda(  0,(double*)(lambda+ijk));
                vector4double            RedBlack_00 = vec_xor(invertMask4,vec_lda( 0,(double*)(RedBlackMask+ij) ));
                vector4double phi_plus_lambda_rhs_00 = vec_madd(lambda_00,rhs_00,phi_00);
                vector4double                 new_00 = vec_nmsub(lambda_00,helmholtz_00,phi_plus_lambda_rhs_00);
                                              new_00 = vec_sel(new_00,phi_00,RedBlack_00);
                                                         vec_sta(new_00, 0,(double*)(phi+ijk));
          ij+=12;ijk+=12;
        }
        }}
      } // if stencil
    } // leadingK
    if(smt_id==0)L1P_SetStreamDepth(old_stream_depth);
    if(smt_id==0)L1P_SetStreamPolicy(L1P_confirmed_or_dcbt);
  } // omp parallel region
}


//==================================================================================================
void smooth(domain_type * domain, int level, int phi_id, int rhs_id, double a, double b){
  int CollaborativeThreadingBoxSize = 100000; // i.e. never
  #ifdef __COLLABORATIVE_THREADING
    CollaborativeThreadingBoxSize = 1 << __COLLABORATIVE_THREADING;
  #endif
  uint64_t _timeStart = CycleTime();
  int box,s;
  int ghosts = domain->ghosts;
  // if communication-avoiding, need RHS for stencils in ghost zones
  if(ghosts>1)exchange_boundary(domain,level,rhs_id,1,1,1);
  for(s=0;s<numSmooths;s+=ghosts){
    exchange_boundary(domain,level,phi_id,1,ghosts>1,ghosts>1);  // corners/edges if doing communication-avoiding...
    if(domain->subdomains[0].levels[level].dim.i >= CollaborativeThreadingBoxSize){
      uint64_t _timeStart = CycleTime();
      for(box=0;box<domain->subdomains_per_rank;box++){__box_smooth_GSRB_multiple_threaded(&domain->subdomains[box].levels[level],phi_id,rhs_id,a,b,s);}
      domain->cycles.smooth[level] += (uint64_t)(CycleTime()-_timeStart);
    }else{
      // now do ghosts communication-avoiding smooths on each box...
      uint64_t _timeStart = CycleTime();
      #pragma omp parallel for private(box)
      for(box=0;box<domain->subdomains_per_rank;box++){__box_smooth_GSRB_multiple(&domain->subdomains[box].levels[level],phi_id,rhs_id,a,b,s);}
      domain->cycles.smooth[level] += (uint64_t)(CycleTime()-_timeStart);
    }
  }
}
//==================================================================================================
