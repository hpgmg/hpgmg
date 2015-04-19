//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int x_id, int rhs_id, double a, double b){
  int block,s;
  for(s=0;s<2*NUM_SMOOTHS;s++){ // there are two sweeps per GSRB smooth

    // exchange the ghost zone...
    #ifdef GSRB_OOP // out-of-place GSRB ping pongs between x and VECTOR_TEMP
    if((s&1)==0){exchange_boundary(level,       x_id,stencil_is_star_shaped());apply_BCs(level,       x_id,stencil_is_star_shaped());}
            else{exchange_boundary(level,VECTOR_TEMP,stencil_is_star_shaped());apply_BCs(level,VECTOR_TEMP,stencil_is_star_shaped());}
    #else // in-place GSRB only operates on x
                 exchange_boundary(level,       x_id,stencil_is_star_shaped());apply_BCs(level,        x_id,stencil_is_star_shaped());
    #endif

    // apply the smoother...
    uint64_t _timeStart = CycleTime();

    // loop over all block/tiles this process owns...
    PRAGMA_THREAD_ACROSS_BLOCKS(level,block,level->num_my_blocks)
    for(block=0;block<level->num_my_blocks;block++){
      const int box = level->my_blocks[block].read.box;
      const int ilo = level->my_blocks[block].read.i;
      const int jlo = level->my_blocks[block].read.j;
      const int klo = level->my_blocks[block].read.k;
      const int ihi = level->my_blocks[block].dim.i + ilo;
      const int jhi = level->my_blocks[block].dim.j + jlo;
      const int khi = level->my_blocks[block].dim.k + klo;
      int i,j,k;
      const int ghosts = level->box_ghosts;
      const int color000 = (level->my_boxes[box].low.i^level->my_boxes[box].low.j^level->my_boxes[box].low.k)&1;  // is element 000 red or black ???  (should only be an issue if box dimension is odd)
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const double h2inv = 1.0/(level->h*level->h);
      const double * __restrict__ rhs      = level->my_boxes[box].vectors[       rhs_id] + ghosts*(1+jStride+kStride);
      const double * __restrict__ alpha    = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_i   = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_j   = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_k   = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);
      const double * __restrict__ Dinv     = level->my_boxes[box].vectors[VECTOR_DINV  ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ valid    = level->my_boxes[box].vectors[VECTOR_VALID ] + ghosts*(1+jStride+kStride); // cell is inside the domain
      #ifdef GSRB_FP
      const double * __restrict__ RedBlack[2] = {level->RedBlack_FP[0] + ghosts*(1+jStride), 
                                                 level->RedBlack_FP[1] + ghosts*(1+jStride)};
      #endif

      #ifdef GSRB_OOP
      const double * __restrict__ x_n;
            double * __restrict__ x_np1;
                     if((s&1)==0){x_n      = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride);
                                  x_np1    = level->my_boxes[box].vectors[VECTOR_TEMP  ] + ghosts*(1+jStride+kStride);}
                             else{x_n      = level->my_boxes[box].vectors[VECTOR_TEMP  ] + ghosts*(1+jStride+kStride);
                                  x_np1    = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride);}
      #else
      const double * __restrict__ x_n      = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
            double * __restrict__ x_np1    = level->my_boxes[box].vectors[         x_id] + ghosts*(1+jStride+kStride); // i.e. [0] = first non ghost zone point
      #endif
          

      #if defined(GSRB_FP)
      #warning GSRB using pre-computed 1.0/0.0 FP array for Red-Black to facilitate vectorization...
      for(k=klo;k<khi;k++){
      for(j=jlo;j<jhi;j++){
      for(i=ilo;i<ihi;i++){
            int EvenOdd = (k^s^color000)&1;
            int ij  = i + j*jStride;
            int ijk = i + j*jStride + k*kStride;
            double Ax     = apply_op_ijk(x_n);
            double lambda =     Dinv_ijk();
            x_np1[ijk] = x_n[ijk] + RedBlack[EvenOdd][ij]*lambda*(rhs[ijk]-Ax); // compiler seems to get confused unless there are disjoint read/write pointers
      }}}
      #elif defined(GSRB_STRIDE2)
      for(k=klo;k<khi;k++){
      for(j=jlo;j<jhi;j++){
        #ifdef GSRB_OOP
        #warning GSRB using out-of-place and stride-2 accesses to minimie the number of flops
        // out-of-place must copy old value...
        for(i=ilo;i<ihi;i++){
          int ijk = i + j*jStride + k*kStride; 
          x_np1[ijk] = x_n[ijk];
        }
        #else
        #warning GSRB using stride-2 accesses to minimie the number of flops
        #endif
        for(i=ilo+((ilo^j^k^s^color000)&1);i<ihi;i+=2){ // stride-2 GSRB
          int ijk = i + j*jStride + k*kStride; 
          double Ax     = apply_op_ijk(x_n);
          double lambda =     Dinv_ijk();
          x_np1[ijk] = x_n[ijk] + lambda*(rhs[ijk]-Ax);
        }
      }}
      #elif defined(GSRB_OOP)
      #warning GSRB using out-of-place implementation with an if-then-else on loop indices...
      for(k=klo;k<khi;k++){
      for(j=jlo;j<jhi;j++){
      for(i=ilo;i<ihi;i++){
        int ijk = i + j*jStride + k*kStride;
        if((i^j^k^s^color000^1)&1){ // looks very clean when [0] is i,j,k=0,0,0 
          double Ax     = apply_op_ijk(x_n);
          double lambda =     Dinv_ijk();
          x_np1[ijk] = x_n[ijk] + lambda*(rhs[ijk]-Ax);
      }else{
          x_np1[ijk] = x_n[ijk]; // copy old value when sweep color != cell color
      }}}}
      #else
      #warning GSRB using if-then-else on loop indices...
      for(k=klo;k<khi;k++){
      for(j=jlo;j<jhi;j++){
      for(i=ilo;i<ihi;i++){
      if((i^j^k^s^color000^1)&1){ // looks very clean when [0] is i,j,k=0,0,0 
            int ijk = i + j*jStride + k*kStride;
            double Ax     = apply_op_ijk(x_n);
            double lambda =     Dinv_ijk();
            x_np1[ijk] = x_n[ijk] + lambda*(rhs[ijk]-Ax);
      }}}}
      #endif
    } // boxes
    level->cycles.smooth += (uint64_t)(CycleTime()-_timeStart);
  } // s-loop
}


//------------------------------------------------------------------------------------------------------------------------------
