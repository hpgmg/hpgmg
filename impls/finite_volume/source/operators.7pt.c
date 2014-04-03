//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//#include<instrinsics.h>
//------------------------------------------------------------------------------------------------------------------------------
#include "timer.h"
#include "defines.h"
#include "level.h"
#include "operators.h"
//------------------------------------------------------------------------------------------------------------------------------
#define __STENCIL_STAR_SHAPED 1
#define __STENCIL_RADIUS      1
#define __STENCIL_VARIABLE_COEFFICIENT
#define __STENCIL_FUSE_BC
//------------------------------------------------------------------------------------------------------------------------------
#ifndef __OMP_COLLAPSE
#define __OMP_COLLAPSE collapse(2)
#endif
//------------------------------------------------------------------------------------------------------------------------------
// FIX... make #define
void apply_BCs(level_type * level, int x_id){
  #ifndef __STENCIL_FUSE_BCs 
  // This is a failure mode if (trying to do communication-avoiding) && (BC!=__BC_PERIODIC)
  apply_BCs_linear(level,x_id);
  #endif
}
//------------------------------------------------------------------------------------------------------------------------------
#ifdef __STENCIL_VARIABLE_COEFFICIENT
  #ifdef __STENCIL_FUSE_BCs                        
    #define __apply_op(x)                                                                        \
    (                                                                                            \
      a*alpha[ijk]*x[ijk] -b*h2inv*(                                                             \
       +beta_i[ijk        ]*( valid[ijk-1      ]*( (x)[ijk] + (x)[ijk-1      ]) - 2.0*(x)[ijk] ) \
       +beta_j[ijk        ]*( valid[ijk-jStride]*( (x)[ijk] + (x)[ijk-jStride]) - 2.0*(x)[ijk] ) \
       +beta_k[ijk        ]*( valid[ijk-kStride]*( (x)[ijk] + (x)[ijk-kStride]) - 2.0*(x)[ijk] ) \
       +beta_i[ijk+1      ]*( valid[ijk+1      ]*( (x)[ijk] + (x)[ijk+1      ]) - 2.0*(x)[ijk] ) \
       +beta_j[ijk+jStride]*( valid[ijk+jStride]*( (x)[ijk] + (x)[ijk+jStride]) - 2.0*(x)[ijk] ) \
       +beta_k[ijk+kStride]*( valid[ijk+kStride]*( (x)[ijk] + (x)[ijk+kStride]) - 2.0*(x)[ijk] ) \
      )                                                                                          \
    )
  #else
    #define __apply_op(x)                                          \
    (                                                              \
      a*alpha[ijk]*x[ijk] -b*h2inv*(                               \
         beta_i[ijk+1      ]*( (x)[ijk+1      ]-(x)[ijk        ] ) \
        -beta_i[ijk        ]*( (x)[ijk        ]-(x)[ijk-1      ] ) \
        +beta_j[ijk+jStride]*( (x)[ijk+jStride]-(x)[ijk        ] ) \
        -beta_j[ijk        ]*( (x)[ijk        ]-(x)[ijk-jStride] ) \
        +beta_k[ijk+kStride]*( (x)[ijk+kStride]-(x)[ijk        ] ) \
        -beta_k[ijk        ]*( (x)[ijk        ]-(x)[ijk-kStride] ) \
      )                                                            \
    )
  #endif
#else
  #error constant coefficient not yet implemented !!!
#endif
//------------------------------------------------------------------------------------------------------------------------------
void rebuild_operator(level_type * level, level_type *fromLevel, double a, double b){
  if(level->my_rank==0){printf("  rebuilding operator for level h=%e: ",level->h);fflush(stdout);}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // form restriction of alpha[], beta_*[] coefficients from fromLevel
  if(fromLevel != NULL){
    restriction(level,__alpha ,fromLevel,__alpha ,__RESTRICT_CELL  );
    restriction(level,__beta_i,fromLevel,__beta_i,__RESTRICT_FACE_I);
    restriction(level,__beta_j,fromLevel,__beta_j,__RESTRICT_FACE_J);
    restriction(level,__beta_k,fromLevel,__beta_k,__RESTRICT_FACE_K);
  } // else case assumes alpha/beta have been set


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // exchange alpha/beta/...  (must be done before calculating Dinv)
  exchange_boundary(level,__alpha ,0); // must be 0(faces,edges,corners) for CA version or 27pt
  exchange_boundary(level,__beta_i,0);
  exchange_boundary(level,__beta_j,0);
  exchange_boundary(level,__beta_k,0);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // calculate Dinv, L1inv, and estimate the dominant Eigenvalue
  uint64_t _timeStart = CycleTime();
  int printedError=0;
  int box;

  double dominant_eigenvalue = -1e9;
  #pragma omp parallel for private(box) num_threads(level->concurrent_boxes) reduction(max:dominant_eigenvalue) schedule(static)
  for(box=0;box<level->num_my_boxes;box++){
    int i,j,k;
    int lowi    = level->my_boxes[box].low.i;
    int lowj    = level->my_boxes[box].low.j;
    int lowk    = level->my_boxes[box].low.k;
    int jStride = level->my_boxes[box].jStride;
    int kStride = level->my_boxes[box].kStride;
    int  ghosts = level->my_boxes[box].ghosts;
    int     dim = level->my_boxes[box].dim;
    double h2inv = 1.0/(level->h*level->h);
    double * __restrict__ alpha  = level->my_boxes[box].components[__alpha ] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_i = level->my_boxes[box].components[__beta_i] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_j = level->my_boxes[box].components[__beta_j] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_k = level->my_boxes[box].components[__beta_k] + ghosts*(1+jStride+kStride);
    double * __restrict__   Dinv = level->my_boxes[box].components[  __Dinv] + ghosts*(1+jStride+kStride);
    double * __restrict__  L1inv = level->my_boxes[box].components[ __L1inv] + ghosts*(1+jStride+kStride);
    double * __restrict__  valid = level->my_boxes[box].components[ __valid] + ghosts*(1+jStride+kStride);
    double box_eigenvalue = -1e9;
    #pragma omp parallel for private(k,j,i) num_threads(level->threads_per_box) reduction(max:box_eigenvalue) schedule(static) __OMP_COLLAPSE
    for(k=0;k<dim;k++){
    for(j=0;j<dim;j++){
    for(i=0;i<dim;i++){
      int ijk = i + j*jStride + k*kStride;
      #if 0
      // FIX This looks wrong, but is faster... theory is because its doing something akin to SOR
      // radius of Gershgorin disc is the sum of the absolute values of the off-diagonal elements...
      double sumAbsAij = fabs(b*h2inv*beta_i[ijk]) + fabs(b*h2inv*beta_i[ijk+      1]) +
                         fabs(b*h2inv*beta_j[ijk]) + fabs(b*h2inv*beta_j[ijk+jStride]) +
                         fabs(b*h2inv*beta_k[ijk]) + fabs(b*h2inv*beta_k[ijk+kStride]);
      // centr of Gershgorin disc is the diagonal element...
      double    Aii = a*alpha[ijk] - b*h2inv*( 
                                       -beta_i[ijk]-beta_i[ijk+      1] 
                                       -beta_j[ijk]-beta_j[ijk+jStride] 
                                       -beta_k[ijk]-beta_k[ijk+kStride] 
                                     );
      #endif
      #if 1
      // radius of Gershgorin disc is the sum of the absolute values of the off-diagonal elements...
      double sumAbsAij = fabs(b*h2inv) * (
                      fabs( beta_i[ijk        ]*valid[ijk-1      ] )+
                      fabs( beta_j[ijk        ]*valid[ijk-jStride] )+
                      fabs( beta_k[ijk        ]*valid[ijk-kStride] )+
                      fabs( beta_i[ijk+1      ]*valid[ijk+1      ] )+
                      fabs( beta_j[ijk+jStride]*valid[ijk+jStride] )+
                      fabs( beta_k[ijk+kStride]*valid[ijk+kStride] )
                      );

      // centr of Gershgorin disc is the diagonal element...
      double    Aii = a*alpha[ijk] - b*h2inv*(
                                       beta_i[ijk        ]*( valid[ijk-1      ]-2.0 )+
                                       beta_j[ijk        ]*( valid[ijk-jStride]-2.0 )+
                                       beta_k[ijk        ]*( valid[ijk-kStride]-2.0 )+
                                       beta_i[ijk+1      ]*( valid[ijk+1      ]-2.0 )+
                                       beta_j[ijk+jStride]*( valid[ijk+jStride]-2.0 )+
                                       beta_k[ijk+kStride]*( valid[ijk+kStride]-2.0 ) 
                                     );

      #endif
      //if( (printedError==0) && (fabs(Aii) <  fabs(sumAbsAij)) ){printf("   not diagonally dominant at (%5d,%5d,%5d) %e vs %e\n",lowi+i,lowj+j,lowk+k,fabs(Aii),fabs(sumAbsAij));printedError=1;}
      //if( (printedError==0) && (fabs(Aii) == fabs(sumAbsAij)) ){printf("weakly diagonally dominant at (%5d,%5d,%5d) %e vs %e\n",lowi+i,lowj+j,lowk+k,fabs(Aii),fabs(sumAbsAij));printedError=1;}
      #if 0
      if( (Aii==0)       ){printf("Error... zero on the diagonal\n");}
      if( isinf(1.0/Aii) ){printf("Error... D^{-1} == inf\n");}
      #endif
       Dinv[ijk] = 1.0/Aii;							// inverse of the diagonal Aii
      L1inv[ijk] = 1.0/(fabs(Aii)+sumAbsAij);					// inverse of the L1 row norm
      double Di = (Aii + sumAbsAij)/Aii;if(Di>box_eigenvalue)box_eigenvalue=Di;	// upper limit to Gershgorin disc == bound on dominant eigenvalue
    }}}
    if(box_eigenvalue>dominant_eigenvalue){dominant_eigenvalue = box_eigenvalue;}
  }
  level->cycles.blas1 += (uint64_t)(CycleTime()-_timeStart);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Reduce the local estimates dominant eigenvalue to a global estimate
  #ifdef __MPI
  uint64_t _timeStartAllReduce = CycleTime();
  double send = dominant_eigenvalue;
  MPI_Allreduce(&send,&dominant_eigenvalue,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  uint64_t _timeEndAllReduce = CycleTime();
  level->cycles.collectives   += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce);
  #endif
  if(level->my_rank==0){printf("  eigenvalue_max <= %e\n",dominant_eigenvalue);fflush(stdout);}
  level->dominant_eigenvalue_of_DinvA = dominant_eigenvalue;


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // exchange Dinv/L1inv/...
  exchange_boundary(level,__Dinv ,0); // must be 0(faces,edges,corners) for CA version
  exchange_boundary(level,__L1inv,0);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}


//------------------------------------------------------------------------------------------------------------------------------
#ifdef  __USE_GSRB
#define __NUM_SMOOTHS      2 // RBRB
#include "operators/gsrb.c"
#elif   __USE_CHEBY
#define __NUM_SMOOTHS      1
#define __CHEBYSHEV_DEGREE 4 // i.e. one degree-4 polynomial smoother
#include "operators/chebyshev.c"
#elif   __USE_JACOBI
#define __NUM_SMOOTHS      6
#include "operators/jacobi.c"
#elif   __USE_L1JACOBI
#define __NUM_SMOOTHS      8
#include "operators/jacobi.c"
#elif   __USE_SYMGS
#define __NUM_SMOOTHS      2
#include "operators/symgs.c"
#else
#error You must compile with either -D__USE_GSRB, -D__USE_CHEBY, or -D__USE_JACOBI
#endif
#include "operators/residual.c"
#include "operators/apply_op.c"
//------------------------------------------------------------------------------------------------------------------------------
#include "operators/blockCopy.c"
#include "operators/misc.c"
#include "operators/exchange_boundary.c"
#include "operators/boundary_conditions.c"
#include "operators/matmul.c"
#include "operators/restriction.c"
#include "operators/interpolation_pc.c"
#include "operators/interpolation_pl.c"
//------------------------------------------------------------------------------------------------------------------------------
void interpolation_vcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_pc(level_f,id_f,prescale_f,level_c,id_c);}
void interpolation_fcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_pl(level_f,id_f,prescale_f,level_c,id_c);}
//------------------------------------------------------------------------------------------------------------------------------
#include "operators/problem.p4.c"
//------------------------------------------------------------------------------------------------------------------------------
