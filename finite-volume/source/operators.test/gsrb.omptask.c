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
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int x_id, int rhs_id, double a, double b){
  int s;
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

    #pragma omp parallel if( (level->num_my_boxes>1) || (level->box_dim > BLOCKCOPY_TILE_J) || (level->box_dim > BLOCKCOPY_TILE_K) )
    {
    int box;
    #pragma omp for private(box) nowait schedule(static,1) // omp for to maximize task injection
    for(box=0;box<level->num_my_boxes;box++){ // loop over all boxes this process owns...
      int kk,jj;
      for(kk=0;kk<level->box_dim;kk+=BLOCKCOPY_TILE_K){
      for(jj=0;jj<level->box_dim;jj+=BLOCKCOPY_TILE_J){
      #pragma omp task firstprivate(kk,jj,box) if( (level->box_dim > BLOCKCOPY_TILE_J) || (level->box_dim > BLOCKCOPY_TILE_K) )
      {
      // either define these variables here, or make them all firstprivate
      const double h2inv = 1.0/(level->h*level->h);
      const int ghosts =  level->box_ghosts;
      const int jStride = level->box_jStride;
      const int kStride = level->box_kStride;
      const int color000 = (level->my_boxes[box].low.i^level->my_boxes[box].low.j^level->my_boxes[box].low.k^s)&1;  // is element 000 red or black on *THIS* sweep
      const double * __restrict__ rhs      = level->my_boxes[box].vectors[       rhs_id] + ghosts*(1+jStride+kStride);
      const double * __restrict__ alpha    = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_i   = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_j   = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_k   = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);
      const double * __restrict__ Dinv     = level->my_boxes[box].vectors[VECTOR_DINV  ] + ghosts*(1+jStride+kStride);
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
      // - - - - 
      const int khi = MIN(kk+BLOCKCOPY_TILE_K,level->box_dim);
      const int jhi = MIN(jj+BLOCKCOPY_TILE_J,level->box_dim);
      int i,j,k;

      #if defined(GSRB_FP)
      for(k=kk;k<khi;k++){
      for(j=jj;j<jhi;j++){
        const double * __restrict__ RedBlack = level->RedBlack_FP + ghosts*(1+jStride) + kStride*((k^color000)&0x1);
        for(i=0;i<level->box_dim;i++){
          int ij  = i + j*jStride;
          int ijk = i + j*jStride + k*kStride;
          double Ax     = apply_op_ijk(x_n);
          double lambda =     Dinv_ijk();
          x_np1[ijk] = x_n[ijk] + RedBlack[ij]*lambda*(rhs[ijk]-Ax);
          //x_np1[ijk] = ((i^j^k^color000)&1) ? x_n[ijk] : x_n[ijk] + lambda*(rhs[ijk]-Ax);
      }}} // i,j,k

      #elif defined(GSRB_STRIDE2)
      for(k=kk;k<khi;k++){
      for(j=jj;j<jhi;j++){
        #ifdef GSRB_OOP
        // out-of-place must copy old value...
        for(i=0;i<level->box_dim;i++){
          int ijk = i + j*jStride + k*kStride; 
          x_np1[ijk] = x_n[ijk];
        } // i copy
        #endif
        for(i=((j^k^color000)&1);i<level->box_dim;i+=2){ // stride-2 GSRB
          int ijk = i + j*jStride + k*kStride; 
          double Ax     = apply_op_ijk(x_n);
          double lambda =     Dinv_ijk();
          x_np1[ijk] = x_n[ijk] + lambda*(rhs[ijk]-Ax);
        } // i stencil
      }}  // j,k

      #elif defined(GSRB_BRANCH)
      for(k=kk;k<khi;k++){
      for(j=jj;j<jhi;j++){
        for(i=0;i<level->box_dim;i++){
          int ijk = i + j*jStride + k*kStride;
          if((i^j^k^color000^1)&1){ // looks very clean when [0] is i,j,k=0,0,0 
            double Ax     = apply_op_ijk(x_n);
            double lambda =     Dinv_ijk();
            x_np1[ijk] = x_n[ijk] + lambda*(rhs[ijk]-Ax);
          #ifdef GSRB_OOP
          }else{
            x_np1[ijk] = x_n[ijk]; // copy old value when sweep color != cell color
          #endif
          }
      }}} // i,j,k

      #else
      #error no GSRB implementation was specified
      #endif

      }}} // JJ,KK,task


    } // boxes
    } // parallel
    level->timers.smooth += (double)(getTime()-_timeStart);
  } // s-loop
}


//------------------------------------------------------------------------------------------------------------------------------
