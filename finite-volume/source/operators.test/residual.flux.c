//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// This version fissions the FV stencil into the 6 fluxes associate with each direction for each cell.  In order to avoid 
// redundant computation, each flux is calculated only once.  However, in order to avoid writing all these fluxes to memory and 
// then rereading them to complete the laplacian, the calculation of fluxes and summation in the laplacian are performed in a 
// pipelined wavefront.  To further enhance performance, the ij loops are fused (ghost zones are clobbered) and OpenMP simd 
// pragmas are utilized.  Finally, compiler specific hints and directives are utilized to facilitate simdization and nontemporal
// stores.
//------------------------------------------------------------------------------------------------------------------------------
#if (BLOCKCOPY_TILE_I != 10000)
#error operators.flux.c cannot block the unit stride dimension (BLOCKCOPY_TILE_I!=10000).
#endif
//------------------------------------------------------------------------------------------------------------------------------
void residual(level_type * level, int res_id, int x_id, int rhs_id, double a, double b){
  // BLOCKCOPY_TILE_J+1 is required for flux_j, but BLOCKCOPY_TILE_J+2 seems to be faster (avoids prefetcher/coherency effects??)
  if(level->fluxes==NULL){posix_memalign( (void**)&(level->fluxes), 512, (4)*(level->num_threads)*(BLOCKCOPY_TILE_J+2)*(level->box_jStride)*sizeof(double) );}

  // exchange the boundary for x in prep for Ax...
  exchange_boundary(level,x_id,stencil_get_shape());
          apply_BCs(level,x_id,stencil_get_shape());

  // now do residual/restriction proper...
  double _timeStart = getTime();
  double h2inv = 1.0/(level->h*level->h);

  // loop over all block/tiles this process owns...
  #pragma omp parallel if(level->num_my_blocks>1)
  {
    int block;
    int    threadID=omp_get_thread_num();
    int num_threads=omp_get_num_threads();

    #if 0 // [flux][thread][ij] layout
    double * __restrict__ flux_i    =  level->fluxes + (0*level->num_threads + threadID)*(BLOCKCOPY_TILE_J+2)*(level->box_jStride);
    double * __restrict__ flux_j    =  level->fluxes + (1*level->num_threads + threadID)*(BLOCKCOPY_TILE_J+2)*(level->box_jStride);
    double * __restrict__ flux_k[2] = {level->fluxes + (2*level->num_threads + threadID)*(BLOCKCOPY_TILE_J+2)*(level->box_jStride),
                                       level->fluxes + (3*level->num_threads + threadID)*(BLOCKCOPY_TILE_J+2)*(level->box_jStride)};
    #else // [thread][flux][ij] layout
    double * __restrict__ flux_i    =  level->fluxes + (4*threadID + 0)*(BLOCKCOPY_TILE_J+2)*(level->box_jStride);
    double * __restrict__ flux_j    =  level->fluxes + (4*threadID + 1)*(BLOCKCOPY_TILE_J+2)*(level->box_jStride);
    double * __restrict__ flux_k[2] = {level->fluxes + (4*threadID + 2)*(BLOCKCOPY_TILE_J+2)*(level->box_jStride),
                                       level->fluxes + (4*threadID + 3)*(BLOCKCOPY_TILE_J+2)*(level->box_jStride)};
    #endif

    #if 1 // static, chunksize == 1... try to exploit constructive locality in the cache
    int blockStart  = threadID;
    int blockEnd    = level->num_my_blocks;
    int blockStride = num_threads;
    #else // static uniform
    int blockStart  = level->num_my_blocks*(threadID  )/num_threads;
    int blockEnd    = level->num_my_blocks*(threadID+1)/num_threads;
    int blockStride = 1;
    #endif

    for(block=blockStart;block<blockEnd;block+=blockStride){
      const int box  = level->my_blocks[block].read.box;
      const int jlo  = level->my_blocks[block].read.j;
      const int klo  = level->my_blocks[block].read.k;
      const int jdim = level->my_blocks[block].dim.j;
      const int kdim = level->my_blocks[block].dim.k;

      const int ghosts  = level->my_boxes[box].ghosts;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;

      const double * __restrict__ rhs    = level->my_boxes[box].vectors[       rhs_id] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
      const double * __restrict__ alpha  = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
      const double * __restrict__ beta_i = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
      const double * __restrict__ beta_j = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
      const double * __restrict__ beta_k = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
      const double * __restrict__ x      = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride); // i.e. [0] = first non ghost zone point
            double * __restrict__ res    = level->my_boxes[box].vectors[       res_id] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);

        #ifdef __INTEL_COMPILER
        // superfluous with OMP4 simd (?)
        //__assume_aligned(x        ,BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(rhs      ,BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(alpha    ,BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(beta_i   ,BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(beta_j   ,BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(beta_k   ,BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(res      ,BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(flux_i   ,BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(flux_j   ,BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(flux_k[0],BOX_ALIGN_JSTRIDE*sizeof(double));
        //__assume_aligned(flux_k[1],BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume(         (jStride) % BOX_ALIGN_JSTRIDE == 0); // e.g. jStride%4==0 or jStride%8==0, hence x+jStride is aligned
        __assume(         (kStride) % BOX_ALIGN_JSTRIDE == 0);
        __assume(               jStride >= BOX_ALIGN_JSTRIDE);
        __assume(               kStride >= BOX_ALIGN_JSTRIDE);
        __assume(                                   jdim > 0);
        __assume(                                   kdim > 0);
        #elif __xlC__
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), x        );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), rhs      );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), alpha    );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_i   );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_j   );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_k   );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), res      );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_i   );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_j   );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_k[0]);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_k[1]);
        #endif


      int ij,k;
      double * __restrict__ flux_klo = flux_k[0];
      // startup / prolog... calculate flux_klo (bottom of cell)...
      #if (_OPENMP>=201307)
      #pragma omp simd aligned(beta_k,x,flux_klo:BOX_ALIGN_JSTRIDE*sizeof(double))
      #endif
      #ifdef __INTEL_COMPILER
      #pragma loop_count min=BOX_ALIGN_JSTRIDE, avg=512
      #endif
      for(ij=0;ij<jdim*jStride;ij++){
        flux_klo[ij] = beta_dxdk(x,ij); // k==0
      }


      // wavefront loop...
      for(k=0;k<kdim;k++){
        double * __restrict__ flux_klo = flux_k[(k  )&0x1];
        double * __restrict__ flux_khi = flux_k[(k+1)&0x1];

        #if 1
        // calculate flux_i and flux_j together
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(beta_i,beta_j,x,flux_i,flux_j:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        #ifdef __INTEL_COMPILER
        #pragma loop_count min=BOX_ALIGN_JSTRIDE, avg=512
        #endif
        for(ij=0;ij<(jdim+1)*jStride;ij++){
          int ijk = ij + k*kStride;
          flux_i[ij] = beta_dxdi(x,ijk);
          flux_j[ij] = beta_dxdj(x,ijk);
        }
        #else
        // calculate flux_i
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(beta_i,x,flux_i:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        #ifdef __INTEL_COMPILER
        #pragma loop_count min=BOX_ALIGN_JSTRIDE, avg=512
        #endif
        for(ij=0;ij<jdim*jStride;ij++){
          int ijk = ij + k*kStride;
          flux_i[ij] = beta_dxdi(x,ijk);
        }
        // calculate flux_j
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(beta_j,x,flux_j:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        #ifdef __INTEL_COMPILER
        #pragma loop_count min=BOX_ALIGN_JSTRIDE, avg=512
        #endif
        for(ij=0;ij<(jdim+1)*jStride;ij++){
          int ijk = ij + k*kStride;
          flux_j[ij] = beta_dxdj(x,ijk);
        }
        #endif


        // calculate flux_khi (top of cell)
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(beta_k,x,flux_khi:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        #ifdef __INTEL_COMPILER
        #pragma loop_count min=BOX_ALIGN_JSTRIDE, avg=512
        #endif
        for(ij=0;ij<jdim*jStride;ij++){
          int ijk = ij + k*kStride;
          flux_khi[ij] = beta_dxdk(x,ijk+kStride);
        }


        // residual...
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(flux_i,flux_j,flux_klo,flux_khi,alpha,rhs,x,res:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        #ifdef __INTEL_COMPILER
        #pragma loop_count min=BOX_ALIGN_JSTRIDE, avg=512
        #pragma vector nontemporal // generally, we don't expect to reuse res
        #endif
        for(ij=0;ij<jdim*jStride;ij++){
          int ijk = ij + k*kStride;
          double Lx = - flux_i[  ij] + flux_i[  ij+      1]
                      - flux_j[  ij] + flux_j[  ij+jStride]
                      - flux_klo[ij] + flux_khi[ij        ];
          #ifdef USE_HELMHOLTZ
          double Ax = a*alpha[ijk]*x[ijk] - b*Lx;
          #else
          double Ax = -b*Lx;
          #endif
          res[ijk] = rhs[ijk]-Ax;
        }

      } // kdim

    } // block
  } // omp
  level->timers.residual += (double)(getTime()-_timeStart);
}

