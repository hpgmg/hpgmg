//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// This routines calculates the residual (res=rhs-Ax) using the linear operator specified in the apply_op_ijk macro
// This requires exchanging a ghost zone and/or enforcing a boundary condition.
// NOTE, x_id must be distinct from rhs_id and res_id
void residual(level_type * level, int res_id, int x_id, int rhs_id, double a, double b){
  if(level->fluxes==NULL){posix_memalign( (void**)&(level->fluxes), 64, level->num_threads*(level->box_jStride)*(BLOCKCOPY_TILE_J+1)*(4)*sizeof(double) );}

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
    int threadID=0;if(level->num_my_blocks>1)threadID = omp_get_thread_num();

    double * __restrict__ flux_i = level->fluxes + (level->box_jStride)*(BLOCKCOPY_TILE_J+1)*( (threadID*4) + 0);
    double * __restrict__ flux_j = level->fluxes + (level->box_jStride)*(BLOCKCOPY_TILE_J+1)*( (threadID*4) + 1);
    double * __restrict__ flux_k = level->fluxes + (level->box_jStride)*(BLOCKCOPY_TILE_J+1)*( (threadID*4) + 2);


    for(block=threadID;block<level->num_my_blocks;block+=level->num_threads){
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

      const double * __restrict__ x      = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride); // i.e. [0] = first non ghost zone point
      const double * __restrict__ rhs    = level->my_boxes[box].vectors[       rhs_id] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
      const double * __restrict__ alpha  = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
      const double * __restrict__ beta_i = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
      const double * __restrict__ beta_j = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
      const double * __restrict__ beta_k = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);
            double * __restrict__ res    = level->my_boxes[box].vectors[       res_id] + ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride);

        #ifdef __INTEL_COMPILER
        __assume_aligned(x     ,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(rhs   ,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(alpha ,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(beta_i,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(beta_j,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(beta_k,BOX_ALIGN_JSTRIDE*sizeof(double));
        __assume_aligned(res   ,BOX_ALIGN_JSTRIDE*sizeof(double));
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
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), x     );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), rhs   );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), alpha );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_i);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_j);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_k);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), res   );
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_i);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_j);
        __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), flux_k);
        #endif

      int i,j,k,ij;
      for(k=0;k<kdim;k++){
        double * __restrict__ flux_klo = flux_k + ((k  )&0x1)*flux_kStride;
        double * __restrict__ flux_khi = flux_k + ((k+1)&0x1)*flux_kStride;

        #if (BLOCKCOPY_TILE_I != 10000)
        #error operators.flux.c cannot block the unit stride dimension (BLOCKCOPY_TILE_I!=10000).
        #endif

        // calculate fluxes (pipeline flux_k)...
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(beta_i,x,flux_i:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        for(ij=0;ij<jdim*jStride;ij++){ // flux_i for jdim pencils...
          int ijk   = ij + (k  )*kStride;
          flux_i[  ij] = beta_dxdi(x,ijk  );
        }
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(beta_j,x,flux_j:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        for(ij=0;ij<(jdim+1)*jStride;ij++){ // flux_j for jdim+1 pencils...
          int ijk   = ij + (k  )*kStride;
          flux_j[  ij] = beta_dxdj(x,ijk  );
        }
        if(k==0){ // startup / prolog for flux_k on jdim pencils...
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(beta_k,x,flux_klo:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        for(ij=0;ij<jdim*jStride;ij++){
          int ijk   = ij + 0;
          flux_klo[ij] = beta_dxdk(x,ijk);
        }}
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(beta_k,x,flux_khi:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        for(ij=0;ij<jdim*jStride;ij++){ // for flux_k on jdim pencils...
          int ijk   = ij + (k+1)*kStride;
          flux_khi[ij] = beta_dxdk(x,ijk); // flux_k needs k+1
        }


        // residual...
        #if (_OPENMP>=201307)
        #pragma omp simd aligned(flux_i,flux_j,flux_klo,flux_khi,alpha,rhs,x,res:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        for(ij=0;ij<(jdim-1)*jStride+idim;ij++){
          int     ijk = ij + k*kStride;
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

