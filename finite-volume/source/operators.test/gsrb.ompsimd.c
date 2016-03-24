//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#if (BLOCKCOPY_TILE_I != 10000)
#error gsrb.simd.c does not permit blocking in the unit stride (BLOCKCOPY_TILE_I != 10000)
#endif
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int x_id, int rhs_id, const double a, const double b){
  int block,s;
  for(s=0;s<2*NUM_SMOOTHS;s++){ // there are two sweeps per GSRB smooth

    // exchange the ghost zone...
    #ifdef GSRB_OOP // out-of-place GSRB ping pongs between x and VECTOR_TEMP
    if((s&1)==0){exchange_boundary(level,       x_id,stencil_get_shape());apply_BCs(level,       x_id,stencil_get_shape());}
            else{exchange_boundary(level,VECTOR_TEMP,stencil_get_shape());apply_BCs(level,VECTOR_TEMP,stencil_get_shape());}
    #else // in-place GSRB only operates on x
                 exchange_boundary(level,       x_id,stencil_get_shape());apply_BCs(level,        x_id,stencil_get_shape());
    #endif

    // apply the smoother...
    double _timeStart = getTime();

    // loop over all block/tiles this process owns...
    PRAGMA_THREAD_ACROSS_BLOCKS(level,block,level->num_my_blocks)
    for(block=0;block<level->num_my_blocks;block++){
      const int box  = level->my_blocks[block].read.box;
    //const int ilo  = level->my_blocks[block].read.i;
      const int jlo  = level->my_blocks[block].read.j;
      const int klo  = level->my_blocks[block].read.k;
      const int idim = level->box_dim;
    //const int idim = level->my_blocks[block].dim.i;
      const int jdim = level->my_blocks[block].dim.j;
      const int kdim = level->my_blocks[block].dim.k;

      const double h2inv = 1.0/(level->h*level->h);
      const int ghosts =  level->box_ghosts;
      const int jStride = level->box_jStride;
      const int kStride = level->box_kStride;
      //const int jStride = level->my_boxes[box].jStride;
      //const int kStride = level->my_boxes[box].kStride;

      const int offset = ghosts*(1+jStride+kStride) + (/*ilo +*/ jlo*jStride + klo*kStride); // offset from first ghost zone to first element this block operates on
      const double * __restrict__ rhs      = level->my_boxes[box].vectors[       rhs_id] + offset;
      const double * __restrict__ alpha    = level->my_boxes[box].vectors[VECTOR_ALPHA ] + offset;
      const double * __restrict__ beta_i   = level->my_boxes[box].vectors[VECTOR_BETA_I] + offset;
      const double * __restrict__ beta_j   = level->my_boxes[box].vectors[VECTOR_BETA_J] + offset;
      const double * __restrict__ beta_k   = level->my_boxes[box].vectors[VECTOR_BETA_K] + offset;
      const double * __restrict__ Dinv     = level->my_boxes[box].vectors[VECTOR_DINV  ] + offset;
      #ifdef GSRB_OOP
      const double * __restrict__ x_n;
            double * __restrict__ x_np1;
                     if((s&1)==0){x_n      = level->my_boxes[box].vectors[         x_id] + offset;
                                  x_np1    = level->my_boxes[box].vectors[VECTOR_TEMP  ] + offset;}
                             else{x_n      = level->my_boxes[box].vectors[VECTOR_TEMP  ] + offset;
                                  x_np1    = level->my_boxes[box].vectors[         x_id] + offset;}
      #else
      const double * __restrict__ x_n      = level->my_boxes[box].vectors[         x_id] + offset;
            double * __restrict__ x_np1    = level->my_boxes[box].vectors[         x_id] + offset;
      #endif

      #ifdef __INTEL_COMPILER
      //__assume_aligned(rhs   ,BOX_ALIGN_JSTRIDE*sizeof(double));
      //__assume_aligned(alpha ,BOX_ALIGN_JSTRIDE*sizeof(double));
      //__assume_aligned(beta_i,BOX_ALIGN_JSTRIDE*sizeof(double));
      //__assume_aligned(beta_j,BOX_ALIGN_JSTRIDE*sizeof(double));
      //__assume_aligned(beta_k,BOX_ALIGN_JSTRIDE*sizeof(double));
      //__assume_aligned(Dinv  ,BOX_ALIGN_JSTRIDE*sizeof(double));
      //__assume_aligned(x_n   ,BOX_ALIGN_JSTRIDE*sizeof(double));
      //__assume_aligned(x_np1 ,BOX_ALIGN_JSTRIDE*sizeof(double));
      __assume( (   jStride) % BOX_ALIGN_JSTRIDE == 0); // e.g. jStride%4==0 or jStride%8==0, hence x+jStride is aligned
      __assume( (   kStride) % BOX_ALIGN_KSTRIDE == 0);
      __assume( ( 2*jStride) % BOX_ALIGN_JSTRIDE == 0);
      __assume( ( 2*kStride) % BOX_ALIGN_KSTRIDE == 0);
      __assume( (      idim) % BOX_ALIGN_JSTRIDE == 0); // FIX... is this safe ?!!
      #elif __xlC__
      __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), rhs   );
      __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), alpha );
      __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_i);
      __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_j);
      __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), beta_k);
      __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), Dinv  );
      __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), x_n   );
      __alignx(BOX_ALIGN_JSTRIDE*sizeof(double), x_np1 );
      #endif
          

      int i,j,k;

      const int color000 = (level->my_boxes[box].low.i^level->my_boxes[box].low.j^level->my_boxes[box].low.k/*^ilo*/^jlo^klo^s)&1;  // is element 000 red or black on *THIS* sweep
      for(k=0;k<kdim;k++){
      for(j=0;j<jdim;j++){
        const double * __restrict__ RedBlack = level->RedBlack_FP + ghosts*(1+jStride) + jStride*((j^k^color000)&1);// exploit pencil symmetry (only need two pencils, not two planes == low cache pressure)
        #if (_OPENMP>=201307) // OpenMP 4.0
        #pragma omp simd aligned(alpha,beta_i,beta_j,beta_k,rhs,Dinv,x_n,x_np1,RedBlack:BOX_ALIGN_JSTRIDE*sizeof(double))
        #endif
        for(i=0;i<idim;i++){ 
          int ijk = i + j*jStride + k*kStride;
          double Ax     = apply_op_ijk(x_n);
          x_np1[ijk] = x_n[ijk] + RedBlack[i]*Dinv[ijk]*(rhs[ijk]-Ax);
      }}}



    } // boxes
    level->timers.smooth += (double)(getTime()-_timeStart);
  } // s-loop
}


//------------------------------------------------------------------------------------------------------------------------------
