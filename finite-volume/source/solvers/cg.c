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
void CG(level_type * level, int x_id, int R_id, double a, double b, double desired_reduction_in_norm){
  // Algorithm 6.18 in Iterative Methods for Sparse Linear Systems(Yousef Saad)
  int  r0_id = VECTORS_RESERVED+0;
  int   r_id = VECTORS_RESERVED+1;
  int   p_id = VECTORS_RESERVED+2;
  int  Ap_id = VECTORS_RESERVED+3;

  int jMax=200;
  int j=0;
  int CGFailed    = 0;
  int CGConverged = 0;
  residual(level,r0_id,x_id,R_id,a,b);                                          // r0[] = R_id[] - A(x_id)
  scale_vector(level,r_id,1.0,r0_id);                                              // r[] = r0[]
  scale_vector(level,p_id,1.0,r0_id);                                              // p[] = r0[]
  double norm_of_r0 = norm(level,r_id);                                         // the norm of the initial residual...
  if(norm_of_r0 == 0.0){CGConverged=1;}                                        // entered CG with exact solution
  double r_dot_r = dot(level,r_id,r_id);                                         // r_dot_r = dot(r,r)
  while( (j<jMax) && (!CGFailed) && (!CGConverged) ){                          // while(not done){
    j++;level->Krylov_iterations++;                                            //
    apply_op(level,Ap_id,p_id,a,b);                                              //   Ap = A(p)
    double Ap_dot_p = dot(level,Ap_id,p_id);                                     //   Ap_dot_p = dot(Ap,p)
    if(Ap_dot_p == 0.0){CGFailed=1;break;}                                     //   pivot breakdown ???
    double alpha = r_dot_r / Ap_dot_p;                                         //   alpha = r_dot_r / Ap_dot_p
    if(isinf(alpha)){CGFailed=1;break;}                                        //   ???
    add_vectors(level,x_id,1.0,x_id, alpha,p_id );                                //   x_id[] = x_id[] + alpha*p[]
    add_vectors(level,r_id,1.0,r_id,-alpha,Ap_id);                                //   r[]    = r[]    - alpha*Ap[]   (intermediate residual?)
    double norm_of_r = norm(level,r_id);                                        //   norm of intermediate residual
    if(norm_of_r == 0.0){CGConverged=1;break;}                                 //
    if(norm_of_r < desired_reduction_in_norm*norm_of_r0){CGConverged=1;break;} //
    double r_dot_r_new = dot(level,r_id,r_id);                                   //   r_dot_r_new = dot(r_{j+1},r_{j+1})
    if(r_dot_r_new == 0.0){CGFailed=1;break;}                                  //   Lanczos breakdown ???
    double beta = (r_dot_r_new/r_dot_r);                                       //   beta = (r_dot_r_new/r_dot_r)
    if(isinf(beta)){CGFailed=1;break;}                                         //   ???
    add_vectors(level,p_id,1.0,r_id,beta,p_id );                                    //   p[] = r[] + beta*p[]
    r_dot_r = r_dot_r_new;                                                     //   r_dot_r = r_dot_r_new   (save old r_dot_r)
  }                                                                            // }
}
