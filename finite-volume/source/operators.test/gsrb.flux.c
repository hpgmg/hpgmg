//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#if   defined(GSRB_FP)
  #warning Overriding default GSRB implementation and using pre-computed 1.0/0.0 FP array for Red-Black to facilitate vectorization...
#elif defined(GSRB_STRIDE2)
  #if defined(GSRB_OOP)
  #warning Overriding default GSRB implementation and using out-of-place and stride-2 accesses to minimize the number of flops
  #else
  #warning Overriding default GSRB implementation and using stride-2 accesses to minimize the number of flops
  #endif
#elif defined(GSRB_BRANCH)
  #if defined(GSRB_OOP)
  #warning Overriding default GSRB implementation and using out-of-place implementation with an if-then-else on loop indices...
  #else
  #warning Overriding default GSRB implementation and using if-then-else on loop indices...
  #endif
#else
#define GSRB_STRIDE2 // default implementation
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include <immintrin.h>
void smooth(level_type * level, int x_id, int rhs_id, double a, double b){
  if(level->fluxes==NULL){posix_memalign( (void**)&(level->fluxes), 64, level->num_threads*(level->box_jStride)*(BLOCKCOPY_TILE_J+1)*(4)*sizeof(double) );}

  int s;
  for(s=0;s<2*NUM_SMOOTHS;s++){ // there are two sweeps per GSRB smooth

    // exchange the ghost zone...
    if((s&1)==0){exchange_boundary(level,       x_id,stencil_get_shape());apply_BCs(level,       x_id,stencil_get_shape());}
            else{exchange_boundary(level,VECTOR_TEMP,stencil_get_shape());apply_BCs(level,VECTOR_TEMP,stencil_get_shape());}

    // apply the smoother...
    double _timeStart = getTime();
    double h2inv = 1.0/(level->h*level->h);

    // loop over all block/tiles this process owns...
    #pragma omp parallel if(level->num_my_blocks>1)
    {
      int block;
      int threadID=0;if(level->num_my_blocks>1)threadID = omp_get_thread_num();

      double * __restrict__ flux_i = level->fluxes + (level->box_jStride)*(BLOCKCOPY_TILE_J+1)*( (threadID*4) + 0);
      double * __restrict__ flux_j = level->fluxes + (level->box_jStride)*(BLOCKCOPY_TILE_J+1)*( (threadID*4) + 1);
      double * __restrict__ flux_k = level->fluxes + (level->box_jStride)*(BLOCKCOPY_TILE_J+1)*( (threadID*4) + 2);


      for(block=threadID;block<level->num_my_blocks;block+=level->num_threads){ // i.e. chunksize=1 block per thread
        const int box  = level->my_blocks[block].read.box;
        const int jlo  = level->my_blocks[block].read.j;
        const int klo  = level->my_blocks[block].read.k;
        const int idim = level->my_blocks[block].dim.i;
        const int jdim = level->my_blocks[block].dim.j;
        const int kdim = level->my_blocks[block].dim.k;

        const int ghosts  = level->my_boxes[box].ghosts;
        const int jStride = level->my_boxes[box].jStride;
        const int kStride = level->my_boxes[box].kStride;
        const int flux_kStride = (BLOCKCOPY_TILE_J+1)*level->box_jStride;

        const int color000 = (level->my_boxes[box].low.i^level->my_boxes[box].low.j^level->my_boxes[box].low.k^jlo^klo^s)&1;  // is element 000 red or black on *THIS* sweep

        const double * __restrict__ rhs      = level->my_boxes[box].vectors[       rhs_id] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
        const double * __restrict__ alpha    = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
        const double * __restrict__ beta_i   = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
        const double * __restrict__ beta_j   = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
        const double * __restrict__ beta_k   = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
        const double * __restrict__ Dinv     = level->my_boxes[box].vectors[VECTOR_DINV  ] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
        const double * __restrict__ x_n;
              double * __restrict__ x_np1;
                       if((s&1)==0){x_n      = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
                                    x_np1    = level->my_boxes[box].vectors[VECTOR_TEMP  ] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);}
                               else{x_n      = level->my_boxes[box].vectors[VECTOR_TEMP  ] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
                                    x_np1    = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);}

        #ifdef __INTEL_COMPILER
        __assume_aligned(rhs   ,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(alpha ,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(beta_i,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(beta_j,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(beta_k,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(Dinv  ,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(x_n   ,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(x_np1 ,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(flux_i,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(flux_j,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(flux_k,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume(        (+jStride) % BOX_ALIGN_JSTRIDE == 0); // e.g. jStride%4==0 or jStride%8==0, hence x+jStride is aligned
        __assume(        (-jStride) % BOX_ALIGN_JSTRIDE == 0);
        __assume(        (+kStride) % BOX_ALIGN_KSTRIDE == 0);
        __assume(        (-kStride) % BOX_ALIGN_KSTRIDE == 0);
        __assume(((jdim  )*jStride) % BOX_ALIGN_JSTRIDE == 0);
        __assume(((jdim+1)*jStride) % BOX_ALIGN_JSTRIDE == 0);
        __assume(    (flux_kStride) % BOX_ALIGN_JSTRIDE == 0);
        #elif __xlC__
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), rhs   );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), alpha );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_i);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_j);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_k);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), Dinv  );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), x_n   );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), x_np1 );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_i);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_j);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_k);
        #endif


        int k;
        for(k=0;k<kdim;k++){
          int i,j,ij;
          double * __restrict__ flux_klo = flux_k + ((k  )&0x1)*flux_kStride;
          double * __restrict__ flux_khi = flux_k + ((k+1)&0x1)*flux_kStride;

          #if (BLOCKCOPY_TILE_I != 10000)
          #error operators.flux.c cannot block the unit stride dimension (BLOCKCOPY_TILE_I!=10000).
          #endif

          // calculate fluxes (pipeline flux_k)...
          #if (_OPENMP>=201307)
          #pragma omp simd aligned(beta_i,x_n,flux_i:BOX_ALIGN_JSTRIDE*sizeof(double))
          #endif
          for(ij=0;ij<jdim*jStride;ij++){
            int ijk = ij + (k  )*kStride;
            flux_i[ij] = beta_dxdi(x_n,ijk);
          }
          #if (_OPENMP>=201307)
          #pragma omp simd aligned(beta_j,x_n,flux_j:BOX_ALIGN_JSTRIDE*sizeof(double))
          #endif
          for(ij=0;ij<(jdim+1)*jStride;ij++){
            int ijk = ij + (k  )*kStride;
            flux_j[ij] = beta_dxdj(x_n,ijk);
          }
          if(k==0){ // startup / prolog for flux_k on jdim pencils...
          #if (_OPENMP>=201307)
          #pragma omp simd aligned(beta_k,x_n,flux_klo:BOX_ALIGN_JSTRIDE*sizeof(double))
          #endif
          for(ij=0;ij<jdim*jStride;ij++){
            int ijk = ij + 0;
            flux_klo[ij] = beta_dxdk(x_n,ijk);
          }}
          #if (_OPENMP>=201307)
          #pragma omp simd aligned(beta_k,x_n,flux_khi:BOX_ALIGN_JSTRIDE*sizeof(double))
          #endif
          for(ij=0;ij<jdim*jStride;ij++){
            int ijk   = ij + (k+1)*kStride;
            flux_khi[ij] = beta_dxdk(x_n,ijk);
          }


          #if defined(GSRB_STRIDE2)
          for(j=0;j<jdim;j++){
            if(((j^k^color000)&1))
            //#if (_OPENMP>=201307)
            //#pragma omp simd aligned(flux_i,flux_j,flux_klo,flux_khi,alpha,rhs,Dinv,x_n,x_np1:BOX_ALIGN_JSTRIDE*sizeof(double))
            //#endif
            for(i=0;i<idim;i+=2){ // stride-2 GSRB
              int ijk = i + j*jStride + k*kStride;
              int ij  = i + j*jStride;
              double Lx = (- flux_i[  ij+1] + flux_i[  ij+2        ]
                           - flux_j[  ij+1] + flux_j[  ij+1+jStride]
                           - flux_klo[ij+1] + flux_khi[ij+1        ] );
              #ifdef USE_HELMHOLTZ
              double Ax = a*alpha[ijk+1]*x_n[ijk+1] - b*Lx;
              #else
              double Ax = -b*Lx;
              #endif
              x_np1[ijk  ] = x_n[ijk  ]; // out-of-place copy
              x_np1[ijk+1] = x_n[ijk+1] + Dinv[ijk+1]*(rhs[ijk+1]-Ax);
            }else
            //#if (_OPENMP>=201307)
            //#pragma omp simd aligned(flux_i,flux_j,flux_klo,flux_khi,alpha,rhs,Dinv,x_n,x_np1:BOX_ALIGN_JSTRIDE*sizeof(double))
            //#endif
            for(i=0;i<idim;i+=2){ // stride-2 GSRB
              int ijk = i + j*jStride + k*kStride;
              int ij  = i + j*jStride;
              double Lx = (- flux_i[  ij] + flux_i[  ij+1      ]
                           - flux_j[  ij] + flux_j[  ij+jStride]
                           - flux_klo[ij] + flux_khi[ij        ] );
              #ifdef USE_HELMHOLTZ
              double Ax = a*alpha[ijk]*x_n[ijk] - b*Lx;
              #else
              double Ax = -b*Lx;
              #endif
              x_np1[ijk  ] = x_n[ijk  ] + Dinv[ijk]*(rhs[ijk]-Ax);
              x_np1[ijk+1] = x_n[ijk+1]; // out-of-place copy
            }
          }
          /*
          */

          #elif defined(GSRB_FP)
          const int color000 = (level->my_boxes[box].low.i^level->my_boxes[box].low.j^level->my_boxes[box].low.k^klo^s)&1;  // is *BLOCK* 000 red or black on *THIS* sweep
          const double * __restrict__ RedBlack = level->RedBlack_FP + ghosts*(1+jStride) + (jlo*jStride) + kStride*((k^color000)&0x1);
          //#if (_OPENMP>=201307)
          //#pragma omp simd aligned(flux_i,flux_j,flux_klo,flux_khi,alpha,rhs,Dinv,x_n,x_np1,RedBlack:BOX_ALIGN_JSTRIDE*sizeof(double))
          //#endif
          for(ij=0;ij<(jdim-1)*jStride+idim;ij++){
            int ijk = ij + k*kStride;
            double Lx = (- flux_i[  ij] + flux_i[  ij+      1]  +
                         - flux_j[  ij] + flux_j[  ij+jStride]  +
                         - flux_klo[ij] + flux_khi[ij        ] );
            #ifdef USE_HELMHOLTZ
            double Ax = a*alpha[ijk]*x_n[ijk] - b*Lx;
            #else
            double Ax = -b*Lx;
            #endif
            x_np1[ijk] = x_n[ijk] + RedBlack[ij]*Dinv[ijk]*(rhs[ijk]-Ax);
          }

          #else
          #error no GSRB implementation was specified
          #endif

        } // kdim

      } // boxes
    } // omp parallel for
    level->timers.smooth += (double)(getTime()-_timeStart);
  } // s-loop
}


//------------------------------------------------------------------------------------------------------------------------------
