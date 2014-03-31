//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// Based on Yousef Saad's Iterative Methods for Sparse Linear Algebra, Algorithm 12.1, page 399
//------------------------------------------------------------------------------------------------------------------------------
void smooth(level_type * level, int x_id, int rhs_id, double a, double b){
  if( (level->dominant_eigenvalue_of_DinvA<=0.0) && (level->my_rank==0) )printf("dominant_eigenvalue_of_DinvA <= 0.0 !\n");
  if((__CHEBYSHEV_DEGREE*__NUM_SMOOTHS)&1){
    printf("error... __CHEBYSHEV_DEGREE*__NUM_SMOOTHS must be even for the chebyshev smoother...\n");
    exit(0);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int box,s;
  int ghosts = level->box_ghosts;
  int radius     = __STENCIL_RADIUS;
  int starShaped = __STENCIL_STAR_SHAPED;
  int communicationAvoiding = ghosts > radius; 


  // compute the Chebyshev coefficients...
  double beta     = 1.000*level->dominant_eigenvalue_of_DinvA;
//double alpha    = 0.300000*beta;
//double alpha    = 0.250000*beta;
//double alpha    = 0.166666*beta;
  double alpha    = 0.125000*beta;
  double theta    = 0.5*(beta+alpha);		// center of the spectral ellipse
  double delta    = 0.5*(beta-alpha);		// major axis?
  double sigma = theta/delta;
  double rho_n = 1/sigma;			// rho_0
  double chebyshev_c1[__CHEBYSHEV_DEGREE];	// + c1*(x_n-x_nm1) == rho_n*rho_nm1
  double chebyshev_c2[__CHEBYSHEV_DEGREE];	// + c2*(b-Ax_n)
  chebyshev_c1[0] = 0.0;
  chebyshev_c2[0] = 1/theta;
  for(s=1;s<__CHEBYSHEV_DEGREE;s++){
    double rho_nm1 = rho_n;
    rho_n = 1.0/(2.0*sigma - rho_nm1);
    chebyshev_c1[s] = rho_n*rho_nm1;
    chebyshev_c2[s] = rho_n*2.0/delta;
  }


  // if communication-avoiding, need updated RHS for stencils in ghost zones
  if(communicationAvoiding)exchange_boundary(level,rhs_id,0); 

  for(s=0;s<__CHEBYSHEV_DEGREE*__NUM_SMOOTHS;s+=ghosts){
    // Chebyshev ping pongs between x_id and __temp
    if((s&1)==0){exchange_boundary(level,  x_id,starShaped && !communicationAvoiding);apply_BCs(level,  x_id);}
            else{exchange_boundary(level,__temp,starShaped && !communicationAvoiding);apply_BCs(level,__temp);}
    
    // now do ghosts communication-avoiding smooths on each box...
    uint64_t _timeStart = CycleTime();
    #pragma omp parallel for private(box) num_threads(level->concurrent_boxes)
    for(box=0;box<level->num_my_boxes;box++){
      int i,j,k,ss;
      const int jStride = level->my_boxes[box].jStride;
      const int kStride = level->my_boxes[box].kStride;
      const int     dim = level->my_boxes[box].dim;
      const double h2inv = 1.0/(level->h*level->h);
      const double * __restrict__ rhs      = level->my_boxes[box].components[  rhs_id] + ghosts*(1+jStride+kStride);
      const double * __restrict__ alpha    = level->my_boxes[box].components[__alpha ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_i   = level->my_boxes[box].components[__beta_i] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_j   = level->my_boxes[box].components[__beta_j] + ghosts*(1+jStride+kStride);
      const double * __restrict__ beta_k   = level->my_boxes[box].components[__beta_k] + ghosts*(1+jStride+kStride);
      const double * __restrict__ lambda   = level->my_boxes[box].components[__Dinv  ] + ghosts*(1+jStride+kStride);
      const double * __restrict__ valid    = level->my_boxes[box].components[__valid ] + ghosts*(1+jStride+kStride); // cell is inside the domain

      int ghostsToOperateOn=ghosts-1;
      for(ss=s;ss<s+ghosts;ss++,ghostsToOperateOn--){
              double * __restrict__ x_np1;
        const double * __restrict__ x_n;
        const double * __restrict__ x_nm1;
              if((ss&1)==0){x_n    = level->my_boxes[box].components[    x_id] + ghosts*(1+jStride+kStride);
                            x_nm1  = level->my_boxes[box].components[  __temp] + ghosts*(1+jStride+kStride); 
                            x_np1  = level->my_boxes[box].components[  __temp] + ghosts*(1+jStride+kStride);}
                       else{x_n    = level->my_boxes[box].components[  __temp] + ghosts*(1+jStride+kStride);
                            x_nm1  = level->my_boxes[box].components[    x_id] + ghosts*(1+jStride+kStride); 
                            x_np1  = level->my_boxes[box].components[    x_id] + ghosts*(1+jStride+kStride);}
        const double c1 = chebyshev_c1[ss%__CHEBYSHEV_DEGREE]; // limit polynomial to degree __CHEBYSHEV_DEGREE.
        const double c2 = chebyshev_c2[ss%__CHEBYSHEV_DEGREE]; // limit polynomial to degree __CHEBYSHEV_DEGREE.
        #pragma omp parallel for private(k,j,i) num_threads(level->threads_per_box) __OMP_COLLAPSE
        for(k=0-ghostsToOperateOn;k<dim+ghostsToOperateOn;k++){
        for(j=0-ghostsToOperateOn;j<dim+ghostsToOperateOn;j++){
        for(i=0-ghostsToOperateOn;i<dim+ghostsToOperateOn;i++){
          int ijk = i + j*jStride + k*kStride;
          // According to Saad... but his was missing a lambda[ijk] == D^{-1} !!!
          //  x_{n+1} = x_{n} + rho_{n} [ rho_{n-1}(x_{n} - x_{n-1}) + (2/delta)(b-Ax_{n}) ]
          //  x_temp[ijk] = x_n[ijk] + c1*(x_n[ijk]-x_temp[ijk]) + c2*lambda[ijk]*(rhs[ijk]-Ax_n);
          double Ax_n = __apply_op(x_n);
          x_np1[ijk] = x_n[ijk] + c1*(x_n[ijk]-x_nm1[ijk]) + c2*lambda[ijk]*(rhs[ijk]-Ax_n);
        }}}
      } // ss-loop
    } // box-loop
    level->cycles.smooth += (uint64_t)(CycleTime()-_timeStart);
  } // s-loop
}
