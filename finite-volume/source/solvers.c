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
#include "defines.h"
#include "level.h"
#include "operators.h"
//------------------------------------------------------------------------------------------------------------------------------
#ifdef USE_BICGSTAB
#include "solvers/bicgstab.c"
#elif  USE_CG
#include "solvers/cg.c"
#elif  USE_CABICGSTAB
#include "solvers/cabicgstab.c"
#elif  USE_CACG
#include "solvers/cacg.c"
#endif
//------------------------------------------------------------------------------------------------------------------------------
void IterativeSolver(level_type * level, int u_id, int f_id, double a, double b, double desired_reduction_in_norm){ 
  if(!level->active)return;
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(level->alpha_is_zero==-1)level->alpha_is_zero = (dot(level,VECTOR_ALPHA,VECTOR_ALPHA) == 0.0);  // haven't determined if alpha[] == 0
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef USE_BICGSTAB
    BiCGStab(level,u_id,f_id,a,b,desired_reduction_in_norm);
  #elif  USE_CG
    CG(level,u_id,f_id,a,b,desired_reduction_in_norm);
  #elif  USE_CABICGSTAB
    CABiCGStab(level,u_id,f_id,a,b,desired_reduction_in_norm);
  #elif  USE_CACG
    CACG(level,u_id,f_id,a,b,desired_reduction_in_norm);
  #else // just point relaxation via multiple smooth()'s
    #if 1 
                     residual(level,VECTOR_TEMP,u_id,f_id,a,b);
                  mul_vectors(level,VECTOR_TEMP,1.0,VECTOR_TEMP,VECTOR_DINV); //  Using ||D^{-1}(b-Ax)||_{inf} as convergence criteria...
     double norm_of_r0 = norm(level,VECTOR_TEMP);
    int s=0,maxSmoothsBottom=10,converged=0;
    while( (s<maxSmoothsBottom) && !converged){
      s++;
      level->Krylov_iterations++;
                       smooth(level,u_id,f_id,a,b);
                     residual(level,VECTOR_TEMP,u_id,f_id,a,b);
                  mul_vectors(level,VECTOR_TEMP,1.0,VECTOR_TEMP,VECTOR_DINV); //  Using ||D^{-1}(b-Ax)||_{inf} as convergence criteria...
      double norm_of_r = norm(level,VECTOR_TEMP);
      if(norm_of_r == 0.0){converged=1;break;}
      if(norm_of_r < desired_reduction_in_norm*norm_of_r0){converged=1;break;}
    }
    #else
    int s=0;int maxSmoothsBottom=10;
    while( (s<maxSmoothsBottom) ){
      smooth(level,u_id,f_id,a,b);
      s++;
    }
    #endif
  #endif
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
}


//------------------------------------------------------------------------------------------------------------------------------
int IterativeSolver_NumVectors(){
  // additionally number of grids required by an iterative solver...
  #ifdef USE_BICGSTAB
  return(8);                  // BiCGStab requires additional grids r0,r,p,s,Ap,As
  #elif  USE_CG
  return(5);                  // CG requires extra grids r0,r,p,Ap,z
  #elif  USE_CABICGSTAB
  return(4+4*CA_KRYLOV_S);    // CABiCGStab requires additional grids rt,p,r,P[2s+1],R[2s].
  #elif  USE_CACG
  return(4+2*CA_KRYLOV_S);    // CACG requires additional grids r0,p,r,P[s+1],R[s].
  #endif
  return(0);                  // simply doing multiple smooths requires no extra grids
}
//------------------------------------------------------------------------------------------------------------------------------
