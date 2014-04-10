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
//#include <hpgmgconf.h>
#include "timer.h"
#include "defines.h"
#include "level.h"
#include "operators.h"
//------------------------------------------------------------------------------------------------------------------------------
#define STENCIL_STAR_SHAPED 1
#define STENCIL_RADIUS      1
//------------------------------------------------------------------------------------------------------------------------------
#define STENCIL_VARIABLE_COEFFICIENT
//#define STENCIL_FUSE_BC
//#define STENCIL_FUSE_DINV
//------------------------------------------------------------------------------------------------------------------------------
#define OMP_THREAD_ACROSS_BOXES(thread_teams    ) if(thread_teams    >1) num_threads(thread_teams    )
#define OMP_THREAD_WITHIN_A_BOX(threads_per_team) if(threads_per_team>1) num_threads(threads_per_team) collapse(2)
//#define OMP_THREAD_ACROSS_BOXES(thread_teams    ) if(0)
//#define OMP_THREAD_WITHIN_A_BOX(threads_per_team) collapse(2)
//#define OMP_THREAD_ACROSS_BOXES(thread_teams    )
//#define OMP_THREAD_WITHIN_A_BOX(threads_per_team) if(0)
//------------------------------------------------------------------------------------------------------------------------------
// fix... make #define...
void apply_BCs(level_type * level, int x_id){
  #ifndef STENCIL_FUSE_BC
  // This is a failure mode if (trying to do communication-avoiding) && (BC!=BC_PERIODIC)
  apply_BCs_linear(level,x_id);
  #endif
}
//------------------------------------------------------------------------------------------------------------------------------
// calculate Dinv?
#ifdef STENCIL_VARIABLE_COEFFICIENT
  #define calculate_Dinv()                                      \
  (                                                             \
    1.0 / (a*alpha[ijk] - b*h2inv*(                             \
             + beta_i[ijk        ]*( valid[ijk-1      ] - 2.0 ) \
             + beta_j[ijk        ]*( valid[ijk-jStride] - 2.0 ) \
             + beta_k[ijk        ]*( valid[ijk-kStride] - 2.0 ) \
             + beta_i[ijk+1      ]*( valid[ijk+1      ] - 2.0 ) \
             + beta_j[ijk+jStride]*( valid[ijk+jStride] - 2.0 ) \
             + beta_k[ijk+kStride]*( valid[ijk+kStride] - 2.0 ) \
          ))                                                    \
  )
#else // constant coefficient case... 
  #define calculate_Dinv()          \
  (                                 \
    1.0 / (a - b*h2inv*(            \
             + valid[ijk-1      ]   \
             + valid[ijk-jStride]   \
             + valid[ijk-kStride]   \
             + valid[ijk+1      ]   \
             + valid[ijk+jStride]   \
             + valid[ijk+kStride]   \
             - 12.0                 \
          ))                        \
  )
#endif

#if defined(STENCIL_FUSE_DINV) && defined(STENCIL_FUSE_BC)
#define Dinv_ijk() calculate_Dinv() // recalculate it
#else
#define Dinv_ijk() Dinv[ijk]        // simply retriev it rather than recalculating it
#endif
//------------------------------------------------------------------------------------------------------------------------------
#ifdef STENCIL_FUSE_BC

  #ifdef STENCIL_VARIABLE_COEFFICIENT
    #define apply_op_ijk(x)                                                                     \
    (                                                                                         \
      a*alpha[ijk]*x[ijk] - b*h2inv*(                                                         \
        + beta_i[ijk        ]*( valid[ijk-1      ]*( x[ijk] + x[ijk-1      ] ) - 2.0*x[ijk] ) \
        + beta_j[ijk        ]*( valid[ijk-jStride]*( x[ijk] + x[ijk-jStride] ) - 2.0*x[ijk] ) \
        + beta_k[ijk        ]*( valid[ijk-kStride]*( x[ijk] + x[ijk-kStride] ) - 2.0*x[ijk] ) \
        + beta_i[ijk+1      ]*( valid[ijk+1      ]*( x[ijk] + x[ijk+1      ] ) - 2.0*x[ijk] ) \
        + beta_j[ijk+jStride]*( valid[ijk+jStride]*( x[ijk] + x[ijk+jStride] ) - 2.0*x[ijk] ) \
        + beta_k[ijk+kStride]*( valid[ijk+kStride]*( x[ijk] + x[ijk+kStride] ) - 2.0*x[ijk] ) \
      )                                                                                       \
    )
  #else  // constant coefficient case...  
    #define apply_op_ijk(x)                                \
    (                                                    \
      a*x[ijk] - b*h2inv*(                               \
        + valid[ijk-1      ]*( x[ijk] + x[ijk-1      ] ) \
        + valid[ijk-jStride]*( x[ijk] + x[ijk-jStride] ) \
        + valid[ijk-kStride]*( x[ijk] + x[ijk-kStride] ) \
        + valid[ijk+1      ]*( x[ijk] + x[ijk+1      ] ) \
        + valid[ijk+jStride]*( x[ijk] + x[ijk+jStride] ) \
        + valid[ijk+kStride]*( x[ijk] + x[ijk+kStride] ) \
                       -12.0*( x[ijk]                  ) \
      )                                                  \
    )
  #endif // variable/constant coefficient

#endif


//------------------------------------------------------------------------------------------------------------------------------
#ifndef STENCIL_FUSE_BC

  #ifdef STENCIL_VARIABLE_COEFFICIENT
    #define apply_op_ijk(x)                                 \
    (                                                     \
      a*alpha[ijk]*x[ijk] - b*h2inv*(                     \
        + beta_i[ijk+1      ]*( x[ijk+1      ] - x[ijk] ) \
        + beta_i[ijk        ]*( x[ijk-1      ] - x[ijk] ) \
        + beta_j[ijk+jStride]*( x[ijk+jStride] - x[ijk] ) \
        + beta_j[ijk        ]*( x[ijk-jStride] - x[ijk] ) \
        + beta_k[ijk+kStride]*( x[ijk+kStride] - x[ijk] ) \
        + beta_k[ijk        ]*( x[ijk-kStride] - x[ijk] ) \
      )                                                   \
    )
  #else  // constant coefficient case...  
    #define apply_op_ijk(x)            \
    (                                \
      a*x[ijk] - b*h2inv*(           \
        + x[ijk+1      ]             \
        + x[ijk-1      ]             \
        + x[ijk+jStride]             \
        + x[ijk-jStride]             \
        + x[ijk+kStride]             \
        + x[ijk-kStride]             \
        - x[ijk        ]*6.0         \
      )                              \
    )
  #endif // variable/constant coefficient

#endif // BCs


//------------------------------------------------------------------------------------------------------------------------------
void rebuild_operator(level_type * level, level_type *fromLevel, double a, double b){
  if(level->my_rank==0){printf("  rebuilding operator for level...  h=%e  ",level->h);fflush(stdout);}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // form restriction of alpha[], beta_*[] coefficients from fromLevel
  if(fromLevel != NULL){
    restriction(level,STENCIL_ALPHA ,fromLevel,STENCIL_ALPHA ,RESTRICT_CELL  );
    restriction(level,STENCIL_BETA_I,fromLevel,STENCIL_BETA_I,RESTRICT_FACE_I);
    restriction(level,STENCIL_BETA_J,fromLevel,STENCIL_BETA_J,RESTRICT_FACE_J);
    restriction(level,STENCIL_BETA_K,fromLevel,STENCIL_BETA_K,RESTRICT_FACE_K);
  } // else case assumes alpha/beta have been set


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // exchange alpha/beta/...  (must be done before calculating Dinv)
  exchange_boundary(level,STENCIL_ALPHA ,0); // must be 0(faces,edges,corners) for CA version or 27pt
  exchange_boundary(level,STENCIL_BETA_I,0);
  exchange_boundary(level,STENCIL_BETA_J,0);
  exchange_boundary(level,STENCIL_BETA_K,0);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // calculate Dinv, L1inv, and estimate the dominant Eigenvalue
  uint64_t _timeStart = CycleTime();
  int printedError=0;
  int box;

  double dominant_eigenvalue = -1e9;
  #pragma omp parallel for private(box) OMP_THREAD_ACROSS_BOXES(level->concurrent_boxes) reduction(max:dominant_eigenvalue) schedule(static)
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
    double * __restrict__ alpha  = level->my_boxes[box].components[STENCIL_ALPHA ] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_i = level->my_boxes[box].components[STENCIL_BETA_I] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_j = level->my_boxes[box].components[STENCIL_BETA_J] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_k = level->my_boxes[box].components[STENCIL_BETA_K] + ghosts*(1+jStride+kStride);
    double * __restrict__   Dinv = level->my_boxes[box].components[STENCIL_DINV  ] + ghosts*(1+jStride+kStride);
    double * __restrict__  L1inv = level->my_boxes[box].components[STENCIL_L1INV ] + ghosts*(1+jStride+kStride);
    double * __restrict__  valid = level->my_boxes[box].components[STENCIL_VALID ] + ghosts*(1+jStride+kStride);
    double box_eigenvalue = -1e9;
    #pragma omp parallel for private(k,j,i) OMP_THREAD_WITHIN_A_BOX(level->threads_per_box) reduction(max:box_eigenvalue) schedule(static)
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
      #if 0
      if( (printedError==0) && (Aii==0)       ){printf("Error(%5d,%5d,%5d)... zero on the diagonal\n",lowi+i,lowj+j,lowk+k);printedError=1}
      if( (printedError==0) && isinf(1.0/Aii) ){printf("Error(%5d,%5d,%5d)... D^{-1} == inf       \n",lowi+i,lowj+j,lowk+k);printedError=1}
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
  #ifdef USE_MPI
  uint64_t _timeStartAllReduce = CycleTime();
  double send = dominant_eigenvalue;
  MPI_Allreduce(&send,&dominant_eigenvalue,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  uint64_t _timeEndAllReduce = CycleTime();
  level->cycles.collectives   += (uint64_t)(_timeEndAllReduce-_timeStartAllReduce);
  #endif
  if(level->my_rank==0){printf("eigenvalue_max<%e\n",dominant_eigenvalue);fflush(stdout);}
  level->dominant_eigenvalue_of_DinvA = dominant_eigenvalue;


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // exchange Dinv/L1inv/...
  exchange_boundary(level,STENCIL_DINV ,0); // must be 0(faces,edges,corners) for CA version
  exchange_boundary(level,STENCIL_L1INV,0);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}


//------------------------------------------------------------------------------------------------------------------------------
#ifdef  USE_GSRB
#define NUM_SMOOTHS      2 // RBRB
#include "operators/gsrb.c"
#elif   USE_CHEBY
#define NUM_SMOOTHS      1
#define CHEBYSHEV_DEGREE 4 // i.e. one degree-4 polynomial smoother
#include "operators/chebyshev.c"
#elif   USE_JACOBI
#define NUM_SMOOTHS      6
#include "operators/jacobi.c"
#elif   USE_L1JACOBI
#define NUM_SMOOTHS      8
#include "operators/jacobi.c"
#elif   USE_SYMGS
#define NUM_SMOOTHS      2
#include "operators/symgs.c"
#else
#error You must compile with either -DUSE_GSRB, -DUSE_CHEBY, or -DUSE_JACOBI
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
