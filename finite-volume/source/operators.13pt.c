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
//------------------------------------------------------------------------------------------------------------------------------
#include "timers.h"
#include "defines.h"
#include "level.h"
#include "operators.h"
//------------------------------------------------------------------------------------------------------------------------------
#define MyPragma(a) _Pragma(#a)
//------------------------------------------------------------------------------------------------------------------------------
#if (_OPENMP>=201107) // OpenMP 3.1 supports max reductions...
  // KNC does not like the num_threads() clause...
  #ifdef __xlC__ // XL C/C++ 12.1.09 sets _OPENMP to 201107, but does not support the max clause
  #define PRAGMA_THREAD_ACROSS_BLOCKS(    level,b,nb     )    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1)                     )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,b,nb,bsum)    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1) reduction(  +:bsum) )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,b,nb,bmax)    
  #warning not threading norm()
  #else
  #define PRAGMA_THREAD_ACROSS_BLOCKS(    level,b,nb     )    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1)                     )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,b,nb,bsum)    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1) reduction(  +:bsum) )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,b,nb,bmax)    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1) reduction(max:bmax) )
  #endif
#elif _OPENMP // older OpenMP versions don't support the max reduction clause
  #warning Threading max reductions requires OpenMP 3.1 (July 2011).  Please upgrade your compiler.                                                           
  #define PRAGMA_THREAD_ACROSS_BLOCKS(    level,b,nb     )    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1)                     )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,b,nb,bsum)    MyPragma(omp parallel for private(b) if(nb>1) schedule(static,1) reduction(  +:bsum) )
  #define PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,b,nb,bmax)    
#else // flat MPI should not define any threading...
  #define PRAGMA_THREAD_ACROSS_BLOCKS(    level,b,nb     )    
  #define PRAGMA_THREAD_ACROSS_BLOCKS_SUM(level,b,nb,bsum)    
  #define PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,b,nb,bmax)    
#endif
//------------------------------------------------------------------------------------------------------------------------------
void apply_BCs(level_type * level, int x_id, int justFaces){
  #warning linear is not sufficient !!!
  apply_BCs_linear(level,x_id,justFaces);
}
//------------------------------------------------------------------------------------------------------------------------------
#define STENCIL_COEF0 (-7.5000000000000000000)  // -90.0 / 12.0
#define STENCIL_COEF1 ( 1.3333333333333333333)	//  16.0 / 12.0
#define STENCIL_COEF2 (-0.0833333333333333333)	//  -1.0 / 12.0
//------------------------------------------------------------------------------------------------------------------------------
#ifdef STENCIL_VARIABLE_COEFFICIENT
  #error This implementation does not support variable-coefficient operators
#endif
#ifdef STENCIL_FUSE_BC
  #error This implementation does not support fusion of the boundary conditions with the operator
#endif
#ifndef USE_PERIODIC_BC
  #error 13pt MUST use -DUSE_PERIODIC_BC as D^{-1} is currently calculated assuming periodic BCs and a high-order BC has not been implemented
#endif
#ifdef USE_HELMHOLTZ
  //#error 13pt MAY not use -DUSE_HELMHOLTZ
#endif
//------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------
#define Dinv_ijk() Dinv[ijk]        // simply retriev it rather than recalculating it
//------------------------------------------------------------------------------------------------------------------------------
#define apply_op_ijk(x)				\
(						\
  a*x[ijk] - b*h2inv*(				\
      STENCIL_COEF2*(x[ijk-2*kStride] +		\
                     x[ijk-2*jStride] +		\
                     x[ijk-2        ] +		\
                     x[ijk+2        ] +		\
                     x[ijk+2*jStride] +		\
                     x[ijk+2*kStride] ) +	\
      STENCIL_COEF1*(x[ijk  -kStride] +		\
                     x[ijk  -jStride] +		\
                     x[ijk  -1      ] +		\
                     x[ijk  +1      ] +		\
                     x[ijk  +jStride] +		\
                     x[ijk  +kStride] ) +	\
      STENCIL_COEF0*(x[ijk          ] )		\
  )						\
)
//------------------------------------------------------------------------------------------------------------------------------
int stencil_get_radius()    {return(2);}
int stencil_is_star_shaped(){return(1);}
//------------------------------------------------------------------------------------------------------------------------------
void rebuild_operator(level_type * level, level_type *fromLevel, double a, double b){
  if(level->my_rank==0){fprintf(stdout,"  rebuilding 13pt CC operator for level...  h=%e  ",level->h);}

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // form restriction of alpha[], beta_*[] coefficients from fromLevel
  if(fromLevel != NULL){
    restriction(level,VECTOR_ALPHA ,fromLevel,VECTOR_ALPHA ,RESTRICT_CELL  );
    restriction(level,VECTOR_BETA_I,fromLevel,VECTOR_BETA_I,RESTRICT_FACE_I);
    restriction(level,VECTOR_BETA_J,fromLevel,VECTOR_BETA_J,RESTRICT_FACE_J);
    restriction(level,VECTOR_BETA_K,fromLevel,VECTOR_BETA_K,RESTRICT_FACE_K);
  } // else case assumes alpha/beta have been set


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // exchange alpha/beta/...  (must be done before calculating Dinv)
  exchange_boundary(level,VECTOR_ALPHA ,0); // must be 0(faces,edges,corners) for CA version or 27pt
  exchange_boundary(level,VECTOR_BETA_I,0);
  exchange_boundary(level,VECTOR_BETA_J,0);
  exchange_boundary(level,VECTOR_BETA_K,0);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // calculate Dinv, L1inv, and estimate the dominant Eigenvalue
  uint64_t _timeStart = CycleTime();
  int block;

  double dominant_eigenvalue = -1e9;

  PRAGMA_THREAD_ACROSS_BLOCKS_MAX(level,block,level->num_my_blocks,dominant_eigenvalue)
  for(block=0;block<level->num_my_blocks;block++){
    const int box = level->my_blocks[block].read.box;
    const int ilo = level->my_blocks[block].read.i;
    const int jlo = level->my_blocks[block].read.j;
    const int klo = level->my_blocks[block].read.k;
    const int ihi = level->my_blocks[block].dim.i + ilo;
    const int jhi = level->my_blocks[block].dim.j + jlo;
    const int khi = level->my_blocks[block].dim.k + klo;
    int i,j,k;
    const int jStride = level->my_boxes[box].jStride;
    const int kStride = level->my_boxes[box].kStride;
    const int  ghosts = level->my_boxes[box].ghosts;
    double h2inv = 1.0/(level->h*level->h);
    double * __restrict__ alpha  = level->my_boxes[box].vectors[VECTOR_ALPHA ] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_i = level->my_boxes[box].vectors[VECTOR_BETA_I] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_j = level->my_boxes[box].vectors[VECTOR_BETA_J] + ghosts*(1+jStride+kStride);
    double * __restrict__ beta_k = level->my_boxes[box].vectors[VECTOR_BETA_K] + ghosts*(1+jStride+kStride);
    double * __restrict__   Dinv = level->my_boxes[box].vectors[VECTOR_DINV  ] + ghosts*(1+jStride+kStride);
    double * __restrict__  L1inv = level->my_boxes[box].vectors[VECTOR_L1INV ] + ghosts*(1+jStride+kStride);
    double * __restrict__  valid = level->my_boxes[box].vectors[VECTOR_VALID ] + ghosts*(1+jStride+kStride);
    double block_eigenvalue = -1e9;

    for(k=klo;k<khi;k++){
    for(j=jlo;j<jhi;j++){
    for(i=ilo;i<ihi;i++){ 
      int ijk = i + j*jStride + k*kStride;
      // radius of Gershgorin disc is the sum of the absolute values of the off-diagonal elements...
                      double sumAbsAij = fabs(b*h2inv*6.0*STENCIL_COEF1) + fabs(b*h2inv*6.0*STENCIL_COEF2);
      // center of Gershgorin disc is the diagonal element...
                            double Aii = a - b*h2inv*( STENCIL_COEF0 );
                             Dinv[ijk] = 1.0/Aii;					// inverse of the diagonal Aii
                          //L1inv[ijk] = 1.0/(Aii+sumAbsAij);				// inverse of the L1 row norm... L1inv = ( D+D^{L1} )^{-1}
      // as suggested by eq 6.5 in Baker et al, "Multigrid smoothers for ultra-parallel computing: additional theory and discussion"...
      if(Aii>=1.5*sumAbsAij)L1inv[ijk] = 1.0/(Aii              ); 			//
                       else L1inv[ijk] = 1.0/(Aii+0.5*sumAbsAij);			// 
      double Di = (Aii + sumAbsAij)/Aii;if(Di>block_eigenvalue)block_eigenvalue=Di;	// upper limit to Gershgorin disc == bound on dominant eigenvalue
    }}}
    if(block_eigenvalue>dominant_eigenvalue){dominant_eigenvalue = block_eigenvalue;}
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
  if(level->my_rank==0){fprintf(stdout,"eigenvalue_max<%e\n",dominant_eigenvalue);}
  level->dominant_eigenvalue_of_DinvA = dominant_eigenvalue;


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // exchange Dinv/L1inv/...
  exchange_boundary(level,VECTOR_DINV ,0); // must be 0(faces,edges,corners) for CA version
  exchange_boundary(level,VECTOR_L1INV,0);
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
}


//------------------------------------------------------------------------------------------------------------------------------
#ifdef  USE_GSRB
#error 13-point operator cannot use GSRB
#elif   USE_CHEBY
#define NUM_SMOOTHS      1
#define CHEBYSHEV_DEGREE 4 // i.e. one degree-4 polynomial smoother
#include "operators/chebyshev.c"
#elif   USE_JACOBI
#define NUM_SMOOTHS      6
#include "operators/jacobi.c"
#elif   USE_L1JACOBI
#define NUM_SMOOTHS      6
#include "operators/jacobi.c"
#elif   USE_SYMGS
#define NUM_SMOOTHS      2
#include "operators/symgs.c"
#else
#error You must compile with either -DUSE_GSRB, -DUSE_CHEBY, -DUSE_JACOBI, -DUSE_L1JACOBI, or -DUSE_SYMGS
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
//#include "operators/interpolation_pq.c"
//------------------------------------------------------------------------------------------------------------------------------
void interpolation_vcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_pc(level_f,id_f,prescale_f,level_c,id_c);}
void interpolation_fcycle(level_type * level_f, int id_f, double prescale_f, level_type *level_c, int id_c){interpolation_pl(level_f,id_f,prescale_f,level_c,id_c);}
//------------------------------------------------------------------------------------------------------------------------------
#include "operators/problem.p6.c"
//------------------------------------------------------------------------------------------------------------------------------
