//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include "../timer.h"
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX_THREADS 256
//------------------------------------------------------------------------------------------------------------------------------
void __box_smooth_GSRB_multiple(box_type *box, int phi_id, int rhs_id, double a, double b, int sweep){
    int pencil = box->pencil;
    int plane = box->plane;
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
    const __m128d       a_splat2 =            _mm_loaddup_pd(&a);
    const __m128d b_h2inv_splat2 = _mm_mul_pd(_mm_loaddup_pd(&b),_mm_loaddup_pd(&h2inv));
    int global_ij_start[8];
    int global_ij_end[8];
    int ij_start[8];
    int ij_end[8];
    int planeInWavefront;for(planeInWavefront=0;planeInWavefront<ghosts;planeInWavefront++){
      global_ij_start[planeInWavefront] = (                   (1+planeInWavefront)*pencil)&~1;
      global_ij_end[planeInWavefront]   = ((ghosts+DimJ+ghosts-1-planeInWavefront)*pencil);
      ij_start[planeInWavefront] = global_ij_start[planeInWavefront];
      ij_end[planeInWavefront]   = global_ij_end[planeInWavefront];
    }
    #if defined(__PREFETCH_NEXT_PLANE_FROM_DRAM)
    double * __restrict__ DRAM_PREFETCH_POINTERS[20];
    DRAM_PREFETCH_POINTERS[0] =    phi+plane-pencil;
    DRAM_PREFETCH_POINTERS[1] = beta_k+plane       ;
    DRAM_PREFETCH_POINTERS[2] = beta_j             ;
    DRAM_PREFETCH_POINTERS[3] = beta_i             ;
    DRAM_PREFETCH_POINTERS[4] =  alpha             ;
    DRAM_PREFETCH_POINTERS[5] =    rhs             ;
    DRAM_PREFETCH_POINTERS[6] = lambda             ;
    #endif
    int leadingK;
    int kLow  =     -(ghosts-1);
    int kHigh = DimK+(ghosts-1);
    for(leadingK=kLow;leadingK<kHigh;leadingK++){
      #if defined(__PREFETCH_NEXT_PLANE_FROM_DRAM)
      int DRAM_prefetch_stream=0;
      if(leadingK>=(kHigh-1))DRAM_prefetch_stream=7; // don't prefetch next plane when on last plane
      int DRAM_prefetch_ijk_start = ij_start[0] + (leadingK+1)*plane;
      int DRAM_prefetch_ijk_end   = ij_end[0]   + (leadingK+1)*plane;
      int DRAM_prefetch_ijk       = DRAM_prefetch_ijk_start;
      #endif
      for(planeInWavefront=0;planeInWavefront<ghosts;planeInWavefront++){
        int k=(leadingK-planeInWavefront);
        if(k>=kLow){
        uint64_t invertMask = 0-((k^planeInWavefront^sweep^1)&0x1);
        const __m128d invertMask2 = _mm_loaddup_pd((double*)&invertMask);
        double * __restrict__ RedBlackFP = box->RedBlack_FP[(k^planeInWavefront^sweep^1)&0x1];
        int kplane=k*plane;
        int ij           = ij_start[planeInWavefront];
        int _ij_end      = ij_end[  planeInWavefront];
        int ijk=ij+kplane;
        while(ij<_ij_end){ // smooth a vector...
          #if defined(__PREFETCH_NEXT_PLANE_FROM_DRAM)
          #warning will attempt to prefetch the next plane from DRAM one component at a time
          if(DRAM_prefetch_stream<7){
            double * _base = DRAM_PREFETCH_POINTERS[DRAM_prefetch_stream] + DRAM_prefetch_ijk;
            _mm_prefetch((const char*)(_base+ 0),_MM_HINT_T1);
            _mm_prefetch((const char*)(_base+ 8),_MM_HINT_T1);
            _mm_prefetch((const char*)(_base+16),_MM_HINT_T1);
            DRAM_prefetch_ijk+=17;
            if(DRAM_prefetch_ijk>DRAM_prefetch_ijk_end){DRAM_prefetch_stream++;DRAM_prefetch_ijk=DRAM_prefetch_ijk_start;}
          }
          #endif
          //careful... assumes the compiler maps _mm128_load_pd to unaligned vmovupd and not the aligned version (should be faster when pencil is a multiple of 2 doubles (16 bytes)
          #if 1
                __m128d helmholtz_00;
                __m128d helmholtz_02;
                __m128d helmholtz_04;
                __m128d helmholtz_06;
                __m128d helmholtz_08;
                                         //_mm_prefetch((const char*)(   phi+ijk+1+10),_MM_HINT_T0);
          const __m128d         phi_m2 = _mm_load_pd(phi+ijk+ -2);
          const __m128d         phi_00 = _mm_load_pd(phi+ijk+  0);
          const __m128d         phi_02 = _mm_load_pd(phi+ijk+  2);
          const __m128d         phi_m1 = _mm_shuffle_pd(phi_m2,phi_00,1);
          const __m128d         phi_04 = _mm_load_pd(phi+ijk+  4);
          const __m128d         phi_01 = _mm_shuffle_pd(phi_00,phi_02,1);
          const __m128d         phi_06 = _mm_load_pd(phi+ijk+  6);
          const __m128d         phi_03 = _mm_shuffle_pd(phi_02,phi_04,1);
          const __m128d         phi_08 = _mm_load_pd(phi+ijk+  8);
          const __m128d         phi_05 = _mm_shuffle_pd(phi_04,phi_06,1);
          const __m128d         phi_10 = _mm_load_pd(phi+ijk+ 10);
          const __m128d         phi_07 = _mm_shuffle_pd(phi_06,phi_08,1);
          const __m128d         phi_09 = _mm_shuffle_pd(phi_08,phi_10,1);
                                         //_mm_prefetch((const char*)(beta_i+ijk+1+10),_MM_HINT_T0);
                        helmholtz_00 =                           _mm_mul_pd(_mm_sub_pd(phi_01,phi_00),_mm_loadu_pd(beta_i+ijk+        1)); 
                        helmholtz_02 =                           _mm_mul_pd(_mm_sub_pd(phi_03,phi_02),_mm_loadu_pd(beta_i+ijk+        3)); 
                        helmholtz_04 =                           _mm_mul_pd(_mm_sub_pd(phi_05,phi_04),_mm_loadu_pd(beta_i+ijk+        5)); 
                        helmholtz_06 =                           _mm_mul_pd(_mm_sub_pd(phi_07,phi_06),_mm_loadu_pd(beta_i+ijk+        7)); 
                        helmholtz_08 =                           _mm_mul_pd(_mm_sub_pd(phi_09,phi_08),_mm_loadu_pd(beta_i+ijk+        9)); 
                        helmholtz_00 = _mm_sub_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(phi_00,phi_m1),_mm_load_pd( beta_i+ijk+        0)));
                        helmholtz_02 = _mm_sub_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(phi_02,phi_01),_mm_load_pd( beta_i+ijk+        2)));
                        helmholtz_04 = _mm_sub_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(phi_04,phi_03),_mm_load_pd( beta_i+ijk+        4)));
                        helmholtz_06 = _mm_sub_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(phi_06,phi_05),_mm_load_pd( beta_i+ijk+        6)));
                        helmholtz_08 = _mm_sub_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(phi_08,phi_07),_mm_load_pd( beta_i+ijk+        8)));
          #else
                __m128d helmholtz_00;
                __m128d helmholtz_02;
                __m128d helmholtz_04;
                __m128d helmholtz_06;
                __m128d helmholtz_08;
                                         //_mm_prefetch((const char*)(   phi+ijk+1+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_i+ijk+1+10),_MM_HINT_T0);
          // this version performs unalligned accesses for phi+/-1, betai+1 and phi+/-pencil
          const __m128d       phi_00 = _mm_load_pd(phi+ijk+  0);
          const __m128d       phi_02 = _mm_load_pd(phi+ijk+  2);
          const __m128d       phi_04 = _mm_load_pd(phi+ijk+  4);
          const __m128d       phi_06 = _mm_load_pd(phi+ijk+  6);
          const __m128d       phi_08 = _mm_load_pd(phi+ijk+  8);
                        helmholtz_00 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        1),             phi_00           ),_mm_loadu_pd(beta_i+ijk+        1)); 
                        helmholtz_02 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        3),             phi_02           ),_mm_loadu_pd(beta_i+ijk+        3)); 
                        helmholtz_04 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        5),             phi_04           ),_mm_loadu_pd(beta_i+ijk+        5)); 
                        helmholtz_06 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        7),             phi_06           ),_mm_loadu_pd(beta_i+ijk+        7)); 
                        helmholtz_08 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        9),             phi_08           ),_mm_loadu_pd(beta_i+ijk+        9)); 
                        helmholtz_00 = _mm_sub_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(             phi_00           ,_mm_loadu_pd(phi+ijk+       -1)),_mm_load_pd( beta_i+ijk+        0)));
                        helmholtz_02 = _mm_sub_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(             phi_02           ,_mm_loadu_pd(phi+ijk+        1)),_mm_load_pd( beta_i+ijk+        2)));
                        helmholtz_04 = _mm_sub_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(             phi_04           ,_mm_loadu_pd(phi+ijk+        3)),_mm_load_pd( beta_i+ijk+        4)));
                        helmholtz_06 = _mm_sub_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(             phi_06           ,_mm_loadu_pd(phi+ijk+        5)),_mm_load_pd( beta_i+ijk+        6)));
                        helmholtz_08 = _mm_sub_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(             phi_08           ,_mm_loadu_pd(phi+ijk+        7)),_mm_load_pd( beta_i+ijk+        8)));
          #endif
                                         //_mm_prefetch((const char*)(   phi+ijk-pencil+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(   phi+ijk+pencil+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_j+ijk       +10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_j+ijk+pencil+10),_MM_HINT_T0);
                        helmholtz_00 = _mm_add_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 0),             phi_00           ),_mm_load_pd( beta_j+ijk+pencil+ 0)));
                        helmholtz_02 = _mm_add_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 2),             phi_02           ),_mm_load_pd( beta_j+ijk+pencil+ 2)));
                        helmholtz_04 = _mm_add_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 4),             phi_04           ),_mm_load_pd( beta_j+ijk+pencil+ 4)));
                        helmholtz_06 = _mm_add_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 6),             phi_06           ),_mm_load_pd( beta_j+ijk+pencil+ 6)));
                        helmholtz_08 = _mm_add_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 8),             phi_08           ),_mm_load_pd( beta_j+ijk+pencil+ 8)));
                        helmholtz_00 = _mm_sub_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(             phi_00           ,_mm_load_pd( phi+ijk-pencil+ 0)),_mm_load_pd( beta_j+ijk+        0)));
                        helmholtz_02 = _mm_sub_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(             phi_02           ,_mm_load_pd( phi+ijk-pencil+ 2)),_mm_load_pd( beta_j+ijk+        2)));
                        helmholtz_04 = _mm_sub_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(             phi_04           ,_mm_load_pd( phi+ijk-pencil+ 4)),_mm_load_pd( beta_j+ijk+        4)));
                        helmholtz_06 = _mm_sub_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(             phi_06           ,_mm_load_pd( phi+ijk-pencil+ 6)),_mm_load_pd( beta_j+ijk+        6)));
                        helmholtz_08 = _mm_sub_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(             phi_08           ,_mm_load_pd( phi+ijk-pencil+ 8)),_mm_load_pd( beta_j+ijk+        8)));
                                         //_mm_prefetch((const char*)(   phi+ijk-plane+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(   phi+ijk+plane+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_k+ijk      +10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_k+ijk+plane+10),_MM_HINT_T0);
                        helmholtz_00 = _mm_add_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 0),             phi_00           ),_mm_load_pd( beta_k+ijk+ plane+ 0)));
                        helmholtz_02 = _mm_add_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 2),             phi_02           ),_mm_load_pd( beta_k+ijk+ plane+ 2)));
                        helmholtz_04 = _mm_add_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 4),             phi_04           ),_mm_load_pd( beta_k+ijk+ plane+ 4)));
                        helmholtz_06 = _mm_add_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 6),             phi_06           ),_mm_load_pd( beta_k+ijk+ plane+ 6)));
                        helmholtz_08 = _mm_add_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 8),             phi_08           ),_mm_load_pd( beta_k+ijk+ plane+ 8)));
                        helmholtz_00 = _mm_sub_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(             phi_00           ,_mm_load_pd( phi+ijk- plane+ 0)),_mm_load_pd( beta_k+ijk       + 0)));
                        helmholtz_02 = _mm_sub_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(             phi_02           ,_mm_load_pd( phi+ijk- plane+ 2)),_mm_load_pd( beta_k+ijk       + 2)));
                        helmholtz_04 = _mm_sub_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(             phi_04           ,_mm_load_pd( phi+ijk- plane+ 4)),_mm_load_pd( beta_k+ijk       + 4)));
                        helmholtz_06 = _mm_sub_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(             phi_06           ,_mm_load_pd( phi+ijk- plane+ 6)),_mm_load_pd( beta_k+ijk       + 6)));
                        helmholtz_08 = _mm_sub_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(             phi_08           ,_mm_load_pd( phi+ijk- plane+ 8)),_mm_load_pd( beta_k+ijk       + 8)));
          #ifdef __GSRB_FP
          #warning GSRB using precomputed 64b FP array for Red-Black
                                         //_mm_prefetch((const char*)(      alpha+ijk+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(        rhs+ijk+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(     lambda+ijk+10),_MM_HINT_T0);
                        helmholtz_00 = _mm_mul_pd(helmholtz_00,b_h2inv_splat2);
                        helmholtz_02 = _mm_mul_pd(helmholtz_02,b_h2inv_splat2);
                        helmholtz_04 = _mm_mul_pd(helmholtz_04,b_h2inv_splat2);
                        helmholtz_06 = _mm_mul_pd(helmholtz_06,b_h2inv_splat2);
                        helmholtz_08 = _mm_mul_pd(helmholtz_08,b_h2inv_splat2);
                        helmholtz_00 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 0)),phi_00),helmholtz_00);
                        helmholtz_02 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 2)),phi_02),helmholtz_02);
                        helmholtz_04 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 4)),phi_04),helmholtz_04);
                        helmholtz_06 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 6)),phi_06),helmholtz_06);
                        helmholtz_08 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 8)),phi_08),helmholtz_08);
                __m128d       new_00 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 0),_mm_sub_pd(helmholtz_00,_mm_load_pd(rhs+ijk+ 0)));
                __m128d       new_02 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 2),_mm_sub_pd(helmholtz_02,_mm_load_pd(rhs+ijk+ 2)));
                __m128d       new_04 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 4),_mm_sub_pd(helmholtz_04,_mm_load_pd(rhs+ijk+ 4)));
                __m128d       new_06 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 6),_mm_sub_pd(helmholtz_06,_mm_load_pd(rhs+ijk+ 6)));
                __m128d       new_08 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 8),_mm_sub_pd(helmholtz_08,_mm_load_pd(rhs+ijk+ 8)));
          const __m128d  RedBlack_00 = _mm_load_pd(RedBlackFP+ij+ 0);
          const __m128d  RedBlack_02 = _mm_load_pd(RedBlackFP+ij+ 2);
          const __m128d  RedBlack_04 = _mm_load_pd(RedBlackFP+ij+ 4);
          const __m128d  RedBlack_06 = _mm_load_pd(RedBlackFP+ij+ 6);
          const __m128d  RedBlack_08 = _mm_load_pd(RedBlackFP+ij+ 8);
                              new_00 = _mm_sub_pd(phi_00,_mm_mul_pd(RedBlack_00,new_00));
                              new_02 = _mm_sub_pd(phi_02,_mm_mul_pd(RedBlack_02,new_02));
                              new_04 = _mm_sub_pd(phi_04,_mm_mul_pd(RedBlack_04,new_04));
                              new_06 = _mm_sub_pd(phi_06,_mm_mul_pd(RedBlack_06,new_06));
                              new_08 = _mm_sub_pd(phi_08,_mm_mul_pd(RedBlack_08,new_08));
                                              _mm_store_pd(phi+ijk+ 0,new_00);
                                              _mm_store_pd(phi+ijk+ 2,new_02);
                                              _mm_store_pd(phi+ijk+ 4,new_04);
                                              _mm_store_pd(phi+ijk+ 6,new_06);
                                              _mm_store_pd(phi+ijk+ 8,new_08);
          ij+=10;
          ijk+=10;
          #else
          #warning GSRB using precomputed 64b integer mask array for Red-Black
                                         //_mm_prefetch((const char*)(      alpha+ijk+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(        rhs+ijk+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(     lambda+ijk+10),_MM_HINT_T0);
                        helmholtz_00 = _mm_mul_pd(helmholtz_00,b_h2inv_splat2);
                        helmholtz_02 = _mm_mul_pd(helmholtz_02,b_h2inv_splat2);
                        helmholtz_04 = _mm_mul_pd(helmholtz_04,b_h2inv_splat2);
                        helmholtz_06 = _mm_mul_pd(helmholtz_06,b_h2inv_splat2);
                        helmholtz_08 = _mm_mul_pd(helmholtz_08,b_h2inv_splat2);
                        helmholtz_00 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 0)),phi_00),helmholtz_00);
                        helmholtz_02 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 2)),phi_02),helmholtz_02);
                        helmholtz_04 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 4)),phi_04),helmholtz_04);
                        helmholtz_06 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 6)),phi_06),helmholtz_06);
                        helmholtz_08 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 8)),phi_08),helmholtz_08);
                __m128d       new_00 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 0),_mm_sub_pd(helmholtz_00,_mm_load_pd(rhs+ijk+ 0)));
                __m128d       new_02 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 2),_mm_sub_pd(helmholtz_02,_mm_load_pd(rhs+ijk+ 2)));
                __m128d       new_04 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 4),_mm_sub_pd(helmholtz_04,_mm_load_pd(rhs+ijk+ 4)));
                __m128d       new_06 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 6),_mm_sub_pd(helmholtz_06,_mm_load_pd(rhs+ijk+ 6)));
                __m128d       new_08 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 8),_mm_sub_pd(helmholtz_08,_mm_load_pd(rhs+ijk+ 8)));
                              new_00 = _mm_sub_pd(phi_00,new_00);
                              new_02 = _mm_sub_pd(phi_02,new_02);
                              new_04 = _mm_sub_pd(phi_04,new_04);
                              new_06 = _mm_sub_pd(phi_06,new_06);
                              new_08 = _mm_sub_pd(phi_08,new_08);
          const __m128d RedBlack_00 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 0) ));
          const __m128d RedBlack_02 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 2) ));
          const __m128d RedBlack_04 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 4) ));
          const __m128d RedBlack_06 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 6) ));
          const __m128d RedBlack_08 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 8) ));
                        new_00  = _mm_or_pd(_mm_and_pd(RedBlack_00,new_00),_mm_andnot_pd(RedBlack_00,phi_00));
                        new_02  = _mm_or_pd(_mm_and_pd(RedBlack_02,new_02),_mm_andnot_pd(RedBlack_02,phi_02));
                        new_04  = _mm_or_pd(_mm_and_pd(RedBlack_04,new_04),_mm_andnot_pd(RedBlack_04,phi_04));
                        new_06  = _mm_or_pd(_mm_and_pd(RedBlack_06,new_06),_mm_andnot_pd(RedBlack_06,phi_06));
                        new_08  = _mm_or_pd(_mm_and_pd(RedBlack_08,new_08),_mm_andnot_pd(RedBlack_08,phi_08));
                                              _mm_store_pd(phi+ijk+ 0,new_00);
                                              _mm_store_pd(phi+ijk+ 2,new_02);
                                              _mm_store_pd(phi+ijk+ 4,new_04);
                                              _mm_store_pd(phi+ijk+ 6,new_06);
                                              _mm_store_pd(phi+ijk+ 8,new_08);
          ij+=10;
          ijk+=10;
          #endif
        }
        } // active plane
      }
    } // leadingK
}


void __box_smooth_GSRB_multiple_threaded(box_type *box, int phi_id, int rhs_id, double a, double b, int sweep){
  volatile int64_t KPlaneFinishedByThread[MAX_THREADS];
  #pragma omp parallel shared(KPlaneFinishedByThread)
  {
    int pencil = box->pencil;
    int plane = box->plane;
    int ghosts = box->ghosts;
    int DimI = box->dim.i;
    int DimJ = box->dim.j;
    int DimK = box->dim.k;
    double h2inv = 1.0/(box->h*box->h);
    double * __restrict__ phi    = box->grids[  phi_id] + ghosts*plane;
    double * __restrict__ rhs    = box->grids[  rhs_id] + ghosts*plane;
    double * __restrict__ alpha  = box->grids[__alpha ] + ghosts*plane;
    double * __restrict__ beta_i = box->grids[__beta_i] + ghosts*plane;
    double * __restrict__ beta_j = box->grids[__beta_j] + ghosts*plane;
    double * __restrict__ beta_k = box->grids[__beta_k] + ghosts*plane;
    double * __restrict__ lambda = box->grids[__lambda] + ghosts*plane;
    const __m128d       a_splat2 =            _mm_loaddup_pd(&a);
    const __m128d b_h2inv_splat2 = _mm_mul_pd(_mm_loaddup_pd(&b),_mm_loaddup_pd(&h2inv));
    int id      = omp_get_thread_num();
    int threads = omp_get_num_threads();
    // only works if (ij_end-ij_start)>=pencil;
    int  left = MAX(        0,id-1);
    int right = MIN(threads-1,id+1);
    if(ghosts==1){right=id;left=id;}
    if(ghosts>1){
      KPlaneFinishedByThread[id]=-100;
      #pragma omp barrier
    }
    int global_ij_start[8];
    int global_ij_end[8];
    int ij_start[8];
    int ij_end[8];
    int planeInWavefront=0;
      global_ij_start[planeInWavefront] = (                   (1)*pencil)&~1;
      global_ij_end[  planeInWavefront] = ((ghosts+DimJ+ghosts-1)*pencil);
      int TotalUnrollings = ((global_ij_end[planeInWavefront]-global_ij_start[planeInWavefront]+10-1)/10);
      ij_start[planeInWavefront] = global_ij_start[planeInWavefront] + 10*( (id          )*(TotalUnrollings)/(threads));
      ij_end[  planeInWavefront] = global_ij_start[planeInWavefront] + 10*( (id+1        )*(TotalUnrollings)/(threads));
      if(ij_end[planeInWavefront]>global_ij_end[planeInWavefront])ij_end[planeInWavefront]=global_ij_end[planeInWavefront];
    for(planeInWavefront=1;planeInWavefront<ghosts;planeInWavefront++){
      ij_start[planeInWavefront] = ij_start[0];
      ij_end[  planeInWavefront] = ij_end[0];
    }
    #if defined(__PREFETCH_NEXT_PLANE_FROM_DRAM)
    double * __restrict__ DRAM_PREFETCH_POINTERS[20];
    DRAM_PREFETCH_POINTERS[0] =    phi+plane-pencil;
    DRAM_PREFETCH_POINTERS[1] = beta_k+plane       ;
    DRAM_PREFETCH_POINTERS[2] = beta_j             ;
    DRAM_PREFETCH_POINTERS[3] = beta_i             ;
    DRAM_PREFETCH_POINTERS[4] =  alpha             ;
    DRAM_PREFETCH_POINTERS[5] =    rhs             ;
    DRAM_PREFETCH_POINTERS[6] = lambda             ;
    #endif
    int leadingK;
    int kLow  =     -(ghosts-1);
    int kHigh = DimK+(ghosts-1);
    for(leadingK=kLow;leadingK<kHigh;leadingK++){
      #if defined(__PREFETCH_NEXT_PLANE_FROM_DRAM)
      int DRAM_prefetch_stream=0;
      if(leadingK>=(kHigh-1))DRAM_prefetch_stream=7; // don't prefetch next plane when on last plane
      int DRAM_prefetch_ijk_start = ij_start[0] + (leadingK+1)*plane;
      int DRAM_prefetch_ijk_end   = ij_end[0]   + (leadingK+1)*plane;
      int DRAM_prefetch_ijk       = DRAM_prefetch_ijk_start;
      #endif
      for(planeInWavefront=0;planeInWavefront<ghosts;planeInWavefront++){
        int k=(leadingK-planeInWavefront);
        if(k>=kLow){
        uint64_t invertMask = 0-((k^planeInWavefront^sweep^1)&0x1);
        const __m128d invertMask2 = _mm_loaddup_pd((double*)&invertMask);
        double * __restrict__ RedBlackFP = box->RedBlack_FP[(k^planeInWavefront^sweep^1)&0x1];
        int kplane=k*plane;
        int ij           = ij_start[planeInWavefront];
        int _ij_end      = ij_end[  planeInWavefront];
        int ijk=ij+kplane;
        while(ij<_ij_end){ // smooth a vector...
          #if defined(__PREFETCH_NEXT_PLANE_FROM_DRAM)
          #warning will attempt to prefetch the next plane from DRAM one component at a time
          if(DRAM_prefetch_stream<7){
            double * _base = DRAM_PREFETCH_POINTERS[DRAM_prefetch_stream] + DRAM_prefetch_ijk;
            _mm_prefetch((const char*)(_base+ 0),_MM_HINT_T1);
            _mm_prefetch((const char*)(_base+ 8),_MM_HINT_T1);
            _mm_prefetch((const char*)(_base+16),_MM_HINT_T1);
            DRAM_prefetch_ijk+=17;
            if(DRAM_prefetch_ijk>DRAM_prefetch_ijk_end){DRAM_prefetch_stream++;DRAM_prefetch_ijk=DRAM_prefetch_ijk_start;}
          }
          #endif
          //careful... assumes the compiler maps _mm128_load_pd to unaligned vmovupd and not the aligned version (should be faster when pencil is a multiple of 2 doubles (16 bytes)
          #if 1
                __m128d helmholtz_00;
                __m128d helmholtz_02;
                __m128d helmholtz_04;
                __m128d helmholtz_06;
                __m128d helmholtz_08;
                                         //_mm_prefetch((const char*)(   phi+ijk+1+10),_MM_HINT_T0);
          const __m128d         phi_m2 = _mm_load_pd(phi+ijk+ -2);
          const __m128d         phi_00 = _mm_load_pd(phi+ijk+  0);
          const __m128d         phi_02 = _mm_load_pd(phi+ijk+  2);
          const __m128d         phi_m1 = _mm_shuffle_pd(phi_m2,phi_00,1);
          const __m128d         phi_04 = _mm_load_pd(phi+ijk+  4);
          const __m128d         phi_01 = _mm_shuffle_pd(phi_00,phi_02,1);
          const __m128d         phi_06 = _mm_load_pd(phi+ijk+  6);
          const __m128d         phi_03 = _mm_shuffle_pd(phi_02,phi_04,1);
          const __m128d         phi_08 = _mm_load_pd(phi+ijk+  8);
          const __m128d         phi_05 = _mm_shuffle_pd(phi_04,phi_06,1);
          const __m128d         phi_10 = _mm_load_pd(phi+ijk+ 10);
          const __m128d         phi_07 = _mm_shuffle_pd(phi_06,phi_08,1);
          const __m128d         phi_09 = _mm_shuffle_pd(phi_08,phi_10,1);
                                         //_mm_prefetch((const char*)(beta_i+ijk+1+10),_MM_HINT_T0);
                        helmholtz_00 =                           _mm_mul_pd(_mm_sub_pd(phi_01,phi_00),_mm_loadu_pd(beta_i+ijk+        1)); 
                        helmholtz_02 =                           _mm_mul_pd(_mm_sub_pd(phi_03,phi_02),_mm_loadu_pd(beta_i+ijk+        3)); 
                        helmholtz_04 =                           _mm_mul_pd(_mm_sub_pd(phi_05,phi_04),_mm_loadu_pd(beta_i+ijk+        5)); 
                        helmholtz_06 =                           _mm_mul_pd(_mm_sub_pd(phi_07,phi_06),_mm_loadu_pd(beta_i+ijk+        7)); 
                        helmholtz_08 =                           _mm_mul_pd(_mm_sub_pd(phi_09,phi_08),_mm_loadu_pd(beta_i+ijk+        9)); 
                        helmholtz_00 = _mm_sub_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(phi_00,phi_m1),_mm_load_pd( beta_i+ijk+        0)));
                        helmholtz_02 = _mm_sub_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(phi_02,phi_01),_mm_load_pd( beta_i+ijk+        2)));
                        helmholtz_04 = _mm_sub_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(phi_04,phi_03),_mm_load_pd( beta_i+ijk+        4)));
                        helmholtz_06 = _mm_sub_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(phi_06,phi_05),_mm_load_pd( beta_i+ijk+        6)));
                        helmholtz_08 = _mm_sub_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(phi_08,phi_07),_mm_load_pd( beta_i+ijk+        8)));
          #else
                __m128d helmholtz_00;
                __m128d helmholtz_02;
                __m128d helmholtz_04;
                __m128d helmholtz_06;
                __m128d helmholtz_08;
                                         //_mm_prefetch((const char*)(   phi+ijk+1+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_i+ijk+1+10),_MM_HINT_T0);
          // this version performs unalligned accesses for phi+/-1, betai+1 and phi+/-pencil
          const __m128d       phi_00 = _mm_load_pd(phi+ijk+  0);
          const __m128d       phi_02 = _mm_load_pd(phi+ijk+  2);
          const __m128d       phi_04 = _mm_load_pd(phi+ijk+  4);
          const __m128d       phi_06 = _mm_load_pd(phi+ijk+  6);
          const __m128d       phi_08 = _mm_load_pd(phi+ijk+  8);
                        helmholtz_00 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        1),             phi_00           ),_mm_loadu_pd(beta_i+ijk+        1)); 
                        helmholtz_02 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        3),             phi_02           ),_mm_loadu_pd(beta_i+ijk+        3)); 
                        helmholtz_04 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        5),             phi_04           ),_mm_loadu_pd(beta_i+ijk+        5)); 
                        helmholtz_06 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        7),             phi_06           ),_mm_loadu_pd(beta_i+ijk+        7)); 
                        helmholtz_08 =                           _mm_mul_pd(_mm_sub_pd(_mm_loadu_pd(phi+ijk+        9),             phi_08           ),_mm_loadu_pd(beta_i+ijk+        9)); 
                        helmholtz_00 = _mm_sub_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(             phi_00           ,_mm_loadu_pd(phi+ijk+       -1)),_mm_load_pd( beta_i+ijk+        0)));
                        helmholtz_02 = _mm_sub_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(             phi_02           ,_mm_loadu_pd(phi+ijk+        1)),_mm_load_pd( beta_i+ijk+        2)));
                        helmholtz_04 = _mm_sub_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(             phi_04           ,_mm_loadu_pd(phi+ijk+        3)),_mm_load_pd( beta_i+ijk+        4)));
                        helmholtz_06 = _mm_sub_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(             phi_06           ,_mm_loadu_pd(phi+ijk+        5)),_mm_load_pd( beta_i+ijk+        6)));
                        helmholtz_08 = _mm_sub_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(             phi_08           ,_mm_loadu_pd(phi+ijk+        7)),_mm_load_pd( beta_i+ijk+        8)));
          #endif
                                         //_mm_prefetch((const char*)(   phi+ijk-pencil+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(   phi+ijk+pencil+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_j+ijk       +10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_j+ijk+pencil+10),_MM_HINT_T0);
                        helmholtz_00 = _mm_add_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 0),             phi_00           ),_mm_load_pd( beta_j+ijk+pencil+ 0)));
                        helmholtz_02 = _mm_add_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 2),             phi_02           ),_mm_load_pd( beta_j+ijk+pencil+ 2)));
                        helmholtz_04 = _mm_add_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 4),             phi_04           ),_mm_load_pd( beta_j+ijk+pencil+ 4)));
                        helmholtz_06 = _mm_add_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 6),             phi_06           ),_mm_load_pd( beta_j+ijk+pencil+ 6)));
                        helmholtz_08 = _mm_add_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+pencil+ 8),             phi_08           ),_mm_load_pd( beta_j+ijk+pencil+ 8)));
                        helmholtz_00 = _mm_sub_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(             phi_00           ,_mm_load_pd( phi+ijk-pencil+ 0)),_mm_load_pd( beta_j+ijk+        0)));
                        helmholtz_02 = _mm_sub_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(             phi_02           ,_mm_load_pd( phi+ijk-pencil+ 2)),_mm_load_pd( beta_j+ijk+        2)));
                        helmholtz_04 = _mm_sub_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(             phi_04           ,_mm_load_pd( phi+ijk-pencil+ 4)),_mm_load_pd( beta_j+ijk+        4)));
                        helmholtz_06 = _mm_sub_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(             phi_06           ,_mm_load_pd( phi+ijk-pencil+ 6)),_mm_load_pd( beta_j+ijk+        6)));
                        helmholtz_08 = _mm_sub_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(             phi_08           ,_mm_load_pd( phi+ijk-pencil+ 8)),_mm_load_pd( beta_j+ijk+        8)));
                                         //_mm_prefetch((const char*)(   phi+ijk-plane+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(   phi+ijk+plane+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_k+ijk      +10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(beta_k+ijk+plane+10),_MM_HINT_T0);
                        helmholtz_00 = _mm_add_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 0),             phi_00           ),_mm_load_pd( beta_k+ijk+ plane+ 0)));
                        helmholtz_02 = _mm_add_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 2),             phi_02           ),_mm_load_pd( beta_k+ijk+ plane+ 2)));
                        helmholtz_04 = _mm_add_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 4),             phi_04           ),_mm_load_pd( beta_k+ijk+ plane+ 4)));
                        helmholtz_06 = _mm_add_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 6),             phi_06           ),_mm_load_pd( beta_k+ijk+ plane+ 6)));
                        helmholtz_08 = _mm_add_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(_mm_load_pd( phi+ijk+ plane+ 8),             phi_08           ),_mm_load_pd( beta_k+ijk+ plane+ 8)));
                        helmholtz_00 = _mm_sub_pd(helmholtz_00,_mm_mul_pd(_mm_sub_pd(             phi_00           ,_mm_load_pd( phi+ijk- plane+ 0)),_mm_load_pd( beta_k+ijk       + 0)));
                        helmholtz_02 = _mm_sub_pd(helmholtz_02,_mm_mul_pd(_mm_sub_pd(             phi_02           ,_mm_load_pd( phi+ijk- plane+ 2)),_mm_load_pd( beta_k+ijk       + 2)));
                        helmholtz_04 = _mm_sub_pd(helmholtz_04,_mm_mul_pd(_mm_sub_pd(             phi_04           ,_mm_load_pd( phi+ijk- plane+ 4)),_mm_load_pd( beta_k+ijk       + 4)));
                        helmholtz_06 = _mm_sub_pd(helmholtz_06,_mm_mul_pd(_mm_sub_pd(             phi_06           ,_mm_load_pd( phi+ijk- plane+ 6)),_mm_load_pd( beta_k+ijk       + 6)));
                        helmholtz_08 = _mm_sub_pd(helmholtz_08,_mm_mul_pd(_mm_sub_pd(             phi_08           ,_mm_load_pd( phi+ijk- plane+ 8)),_mm_load_pd( beta_k+ijk       + 8)));
          #ifdef __GSRB_FP
          #warning GSRB using precomputed 64b FP array for Red-Black
                                         //_mm_prefetch((const char*)(      alpha+ijk+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(        rhs+ijk+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(     lambda+ijk+10),_MM_HINT_T0);
                        helmholtz_00 = _mm_mul_pd(helmholtz_00,b_h2inv_splat2);
                        helmholtz_02 = _mm_mul_pd(helmholtz_02,b_h2inv_splat2);
                        helmholtz_04 = _mm_mul_pd(helmholtz_04,b_h2inv_splat2);
                        helmholtz_06 = _mm_mul_pd(helmholtz_06,b_h2inv_splat2);
                        helmholtz_08 = _mm_mul_pd(helmholtz_08,b_h2inv_splat2);
                        helmholtz_00 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 0)),phi_00),helmholtz_00);
                        helmholtz_02 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 2)),phi_02),helmholtz_02);
                        helmholtz_04 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 4)),phi_04),helmholtz_04);
                        helmholtz_06 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 6)),phi_06),helmholtz_06);
                        helmholtz_08 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 8)),phi_08),helmholtz_08);
                __m128d       new_00 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 0),_mm_sub_pd(helmholtz_00,_mm_load_pd(rhs+ijk+ 0)));
                __m128d       new_02 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 2),_mm_sub_pd(helmholtz_02,_mm_load_pd(rhs+ijk+ 2)));
                __m128d       new_04 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 4),_mm_sub_pd(helmholtz_04,_mm_load_pd(rhs+ijk+ 4)));
                __m128d       new_06 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 6),_mm_sub_pd(helmholtz_06,_mm_load_pd(rhs+ijk+ 6)));
                __m128d       new_08 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 8),_mm_sub_pd(helmholtz_08,_mm_load_pd(rhs+ijk+ 8)));
          const __m128d  RedBlack_00 = _mm_load_pd(RedBlackFP+ij+ 0);
          const __m128d  RedBlack_02 = _mm_load_pd(RedBlackFP+ij+ 2);
          const __m128d  RedBlack_04 = _mm_load_pd(RedBlackFP+ij+ 4);
          const __m128d  RedBlack_06 = _mm_load_pd(RedBlackFP+ij+ 6);
          const __m128d  RedBlack_08 = _mm_load_pd(RedBlackFP+ij+ 8);
                              new_00 = _mm_sub_pd(phi_00,_mm_mul_pd(RedBlack_00,new_00));
                              new_02 = _mm_sub_pd(phi_02,_mm_mul_pd(RedBlack_02,new_02));
                              new_04 = _mm_sub_pd(phi_04,_mm_mul_pd(RedBlack_04,new_04));
                              new_06 = _mm_sub_pd(phi_06,_mm_mul_pd(RedBlack_06,new_06));
                              new_08 = _mm_sub_pd(phi_08,_mm_mul_pd(RedBlack_08,new_08));
                                              _mm_store_pd(phi+ijk+ 0,new_00);
                                              _mm_store_pd(phi+ijk+ 2,new_02);
                                              _mm_store_pd(phi+ijk+ 4,new_04);
                                              _mm_store_pd(phi+ijk+ 6,new_06);
                                              _mm_store_pd(phi+ijk+ 8,new_08);
          ij+=10;
          ijk+=10;
          #else
          #warning GSRB using precomputed 64b integer mask array for Red-Black
                                         //_mm_prefetch((const char*)(      alpha+ijk+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(        rhs+ijk+10),_MM_HINT_T0);
                                         //_mm_prefetch((const char*)(     lambda+ijk+10),_MM_HINT_T0);
                        helmholtz_00 = _mm_mul_pd(helmholtz_00,b_h2inv_splat2);
                        helmholtz_02 = _mm_mul_pd(helmholtz_02,b_h2inv_splat2);
                        helmholtz_04 = _mm_mul_pd(helmholtz_04,b_h2inv_splat2);
                        helmholtz_06 = _mm_mul_pd(helmholtz_06,b_h2inv_splat2);
                        helmholtz_08 = _mm_mul_pd(helmholtz_08,b_h2inv_splat2);
                        helmholtz_00 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 0)),phi_00),helmholtz_00);
                        helmholtz_02 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 2)),phi_02),helmholtz_02);
                        helmholtz_04 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 4)),phi_04),helmholtz_04);
                        helmholtz_06 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 6)),phi_06),helmholtz_06);
                        helmholtz_08 = _mm_sub_pd(_mm_mul_pd(_mm_mul_pd(a_splat2,_mm_load_pd(alpha+ijk+ 8)),phi_08),helmholtz_08);
                __m128d       new_00 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 0),_mm_sub_pd(helmholtz_00,_mm_load_pd(rhs+ijk+ 0)));
                __m128d       new_02 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 2),_mm_sub_pd(helmholtz_02,_mm_load_pd(rhs+ijk+ 2)));
                __m128d       new_04 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 4),_mm_sub_pd(helmholtz_04,_mm_load_pd(rhs+ijk+ 4)));
                __m128d       new_06 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 6),_mm_sub_pd(helmholtz_06,_mm_load_pd(rhs+ijk+ 6)));
                __m128d       new_08 = _mm_mul_pd(_mm_load_pd(lambda+ijk+ 8),_mm_sub_pd(helmholtz_08,_mm_load_pd(rhs+ijk+ 8)));
                              new_00 = _mm_sub_pd(phi_00,new_00);
                              new_02 = _mm_sub_pd(phi_02,new_02);
                              new_04 = _mm_sub_pd(phi_04,new_04);
                              new_06 = _mm_sub_pd(phi_06,new_06);
                              new_08 = _mm_sub_pd(phi_08,new_08);
          const __m128d RedBlack_00 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 0) ));
          const __m128d RedBlack_02 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 2) ));
          const __m128d RedBlack_04 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 4) ));
          const __m128d RedBlack_06 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 6) ));
          const __m128d RedBlack_08 = _mm_xor_pd(invertMask2,(__m128d)_mm_load_si128( (__m128i*)(box->RedBlack_64bMask+ij+ 8) ));
                        new_00  = _mm_or_pd(_mm_and_pd(RedBlack_00,new_00),_mm_andnot_pd(RedBlack_00,phi_00));
                        new_02  = _mm_or_pd(_mm_and_pd(RedBlack_02,new_02),_mm_andnot_pd(RedBlack_02,phi_02));
                        new_04  = _mm_or_pd(_mm_and_pd(RedBlack_04,new_04),_mm_andnot_pd(RedBlack_04,phi_04));
                        new_06  = _mm_or_pd(_mm_and_pd(RedBlack_06,new_06),_mm_andnot_pd(RedBlack_06,phi_06));
                        new_08  = _mm_or_pd(_mm_and_pd(RedBlack_08,new_08),_mm_andnot_pd(RedBlack_08,phi_08));
                                              _mm_store_pd(phi+ijk+ 0,new_00);
                                              _mm_store_pd(phi+ijk+ 2,new_02);
                                              _mm_store_pd(phi+ijk+ 4,new_04);
                                              _mm_store_pd(phi+ijk+ 6,new_06);
                                              _mm_store_pd(phi+ijk+ 8,new_08);
          ij+=10;
          ijk+=10;
          #endif
        }
        } // active plane
      }
      if(ghosts>1){
        KPlaneFinishedByThread[id]=leadingK;
        while( (KPlaneFinishedByThread[left ]<leadingK) || (KPlaneFinishedByThread[right]<leadingK) ){_mm_pause();}; // pause() in case HT is in use...
      }
    } // leadingK
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
