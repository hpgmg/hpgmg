//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#if (BLOCKCOPY_TILE_I != 10000)
#error gsrb.simd.c does not permit blocking in the unit stride (BLOCKCOPY_TILE_I != 10000)
#endif
//------------------------------------------------------------------------------------------------------------------------------
void residual(level_type * level, int res_id, int x_id, int rhs_id, double a, double b){
  int block;

  // exchange the boundary for x in prep for Ax...
  exchange_boundary(level,x_id,stencil_get_shape());
          apply_BCs(level,x_id,stencil_get_shape());

  // now do residual/restriction proper...
  double _timeStart = getTime();
  double h2inv = 1.0/(level->h*level->h);

  // loop over all block/tiles this process owns...
  PRAGMA_THREAD_ACROSS_BLOCKS(level,block,level->num_my_blocks)
  for(block=0;block<level->num_my_blocks;block++){
    const int box  = level->my_blocks[block].read.box;
    const int jlo  = level->my_blocks[block].read.j;
    const int klo  = level->my_blocks[block].read.k;
    const int idim = level->my_blocks[block].dim.i;
    const int jdim = level->my_blocks[block].dim.j;
    const int kdim = level->my_blocks[block].dim.k;

    const double h2inv = 1.0/(level->h*level->h);
    const int ghosts =  level->box_ghosts;
    const int jStride = level->box_jStride;
    const int kStride = level->box_kStride;
    //const int jStride = level->my_boxes[box].jStride;
    //const int kStride = level->my_boxes[box].kStride;

    const int offset = ghosts*(1+jStride+kStride) + (jlo*jStride + klo*kStride); // offset from first ghost zone to first element this block operates on
    const double * __restrict__ rhs      = level->my_boxes[box].vectors[       rhs_id] + offset;
    #ifdef VECTOR_ALPHA
    const double * __restrict__ alpha    = level->my_boxes[box].vectors[VECTOR_ALPHA ] + offset;
    #endif
    const double * __restrict__ beta_i   = level->my_boxes[box].vectors[VECTOR_BETA_I] + offset;
    const double * __restrict__ beta_j   = level->my_boxes[box].vectors[VECTOR_BETA_J] + offset;
    const double * __restrict__ beta_k   = level->my_boxes[box].vectors[VECTOR_BETA_K] + offset;
    const double * __restrict__ x        = level->my_boxes[box].vectors[         x_id] + offset;
          double * __restrict__ res      = level->my_boxes[box].vectors[       res_id] + offset;

    #ifdef __INTEL_COMPILER
    //__assume_aligned(rhs   ,BOX_ALIGN_JSTRIDE*sizeof(double));
    //__assume_aligned(alpha ,BOX_ALIGN_JSTRIDE*sizeof(double));
    //__assume_aligned(beta_i,BOX_ALIGN_JSTRIDE*sizeof(double));
    //__assume_aligned(beta_j,BOX_ALIGN_JSTRIDE*sizeof(double));
    //__assume_aligned(beta_k,BOX_ALIGN_JSTRIDE*sizeof(double));
    //__assume_aligned(x     ,BOX_ALIGN_JSTRIDE*sizeof(double));
    //__assume_aligned(res   ,BOX_ALIGN_JSTRIDE*sizeof(double));
    __assume(       (  jStride) % BOX_ALIGN_JSTRIDE == 0); // e.g. jStride%4==0 or jStride%8==0, hence x+jStride is aligned
    __assume(       (  kStride) % BOX_ALIGN_JSTRIDE == 0);
    __assume(               jStride >= BOX_ALIGN_JSTRIDE);
    __assume(               kStride >= BOX_ALIGN_JSTRIDE);
    __assume(                                   idim > 0);
    __assume(                                   jdim > 0);
    __assume(                                   kdim > 0);
    #elif __xlC__
    __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), rhs   );
    __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), alpha );
    __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_i);
    __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_j);
    __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_k);
    __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), x     );
    __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), res   );
    #endif
        

    int i,j,k,ij;

    #if 0
    for(k=0;k<kdim;k++){
    for(j=0;j<jdim;j++){
      #if (_OPENMP>=201307) // OpenMP 4.0
      #pragma omp simd aligned(alpha,beta_i,beta_j,beta_k,rhs,x,res:BOX_ALIGN_JSTRIDE*sizeof(double))
      #endif
      #pragma loop_count min=BOX_ALIGN_JSTRIDE, avg=128
      #pragma vector nontemporal // generally, we don't expect to reuse res
      for(i=0;i<idim;i++){
        int ijk = i + j*jStride + k*kStride;
        double Ax     = apply_op_ijk(x);
        res[ijk] = rhs[ijk]-Ax;
    }}}
    #else // fused IJ loop
    for(k=0;k<kdim;k++){
      #if (_OPENMP>=201307) // OpenMP 4.0
      #pragma omp simd aligned(alpha,beta_i,beta_j,beta_k,rhs,x,res:BOX_ALIGN_JSTRIDE*sizeof(double))
      #endif
      #pragma loop_count min=BOX_ALIGN_JSTRIDE, avg=512
      #pragma vector nontemporal // generally, we don't expect to reuse res
      for(ij=0;ij<jdim*jStride;ij++){  // This walks into the ghost zone but should be completely SIMDizable
        int ijk = ij + k*kStride;
        double Ax     = apply_op_ijk(x);
        res[ijk] = rhs[ijk]-Ax;
    }}
    #endif



  } // boxes
  level->timers.residual += (double)(getTime()-_timeStart);
}


//------------------------------------------------------------------------------------------------------------------------------
