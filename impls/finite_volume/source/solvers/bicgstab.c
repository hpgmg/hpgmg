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
#define DiagonallyPrecondition
//------------------------------------------------------------------------------------------------------------------------------
void BiCGStab(level_type * level, int x_id, int R_id, double a, double b, double desired_reduction_in_norm){
  // Algorithm 7.7 in Iterative Methods for Sparse Linear Systems(Yousef Saad)
  // uses "right" preconditioning...  AD^{-1}(Dx) = b ... AD^{-1}y = b ... solve for y, then solve for x = D^{-1}y
  int  __r0 = __Components+0;
  int  __r  = __Components+1;
  int  __p  = __Components+2;
  int  __s  = __Components+3;
  int  __Ap = __Components+4;
  int  __As = __Components+5;

  int jMax=200;
  int j=0;
  int BiCGStabFailed    = 0;
  int BiCGStabConverged = 0;
  residual(level,__r0,x_id,R_id,a,b);                                           // r0[] = R_id[] - A(x_id)
  scale_grid(level,__r,1.0,__r0);                                               // r[] = r0[]
  scale_grid(level,__p,1.0,__r0);                                               // p[] = r0[]
  double r_dot_r0 = dot(level,__r,__r0);                                        // r_dot_r0 = dot(r,r0)
  double norm_of_r0 = norm(level,__r);                                          // the norm of the initial residual...
  if(r_dot_r0   == 0.0){BiCGStabConverged=1;}                                   // entered BiCGStab with exact solution
  if(norm_of_r0 == 0.0){BiCGStabConverged=1;}                                   // entered BiCGStab with exact solution
  while( (j<jMax) && (!BiCGStabFailed) && (!BiCGStabConverged) ){               // while(not done){
    j++;level->Krylov_iterations++;                                             //
    #ifdef DiagonallyPrecondition                                               //
    mul_grids(level,__temp,1.0,__Dinv,__p);                                     //   temp[] = Dinv[]*p[]
    apply_op(level,__Ap,__temp,a,b);                                            //   Ap = AD^{-1}(p)
    #else                                                                       //
    apply_op(level,__Ap,__p,a,b);                                               //   Ap = A(p)
    #endif                                                                      //
    double Ap_dot_r0 = dot(level,__Ap,__r0);                                    //   Ap_dot_r0 = dot(Ap,r0)
    if(Ap_dot_r0 == 0.0){BiCGStabFailed=1;break;}                               //   pivot breakdown ???
    double alpha = r_dot_r0 / Ap_dot_r0;                                        //   alpha = r_dot_r0 / Ap_dot_r0
    if(isinf(alpha)){BiCGStabFailed=2;break;}                                   //   pivot breakdown ???
    add_grids(level,x_id,1.0,x_id, alpha,__p );                                 //   x_id[] = x_id[] + alpha*p[]
    add_grids(level,__s ,1.0,__r ,-alpha,__Ap);                                 //   s[]    = r[]    - alpha*Ap[]   (intermediate residual?)
    double norm_of_s = norm(level,__s);                                         //   FIX - redundant??  norm of intermediate residual
  //if(level->my_rank==0)printf("norm(s)/norm(r0) = %e\n",norm_of_s/norm_of_r0);
    if(norm_of_s == 0.0){BiCGStabConverged=1;break;}                            //   FIX - redundant??  if As_dot_As==0, then As must be 0 which implies s==0
    if(norm_of_s < desired_reduction_in_norm*norm_of_r0){BiCGStabConverged=1;break;}
    #ifdef DiagonallyPrecondition                                               //
    mul_grids(level,__temp,1.0,__Dinv,__s);                                     //   temp[] = Dinv[]*s[]
    apply_op(level,__As,__temp,a,b);                                            //   As = AD^{-1}(s)
    #else                                                                       //
    apply_op(level,__As,__s,a,b);                                               //   As = A(s)
    #endif                                                                      //
    double As_dot_As = dot(level,__As,__As);                                    //   As_dot_As = dot(As,As)
    double As_dot_s  = dot(level,__As, __s);                                    //   As_dot_s  = dot(As, s)
    if(As_dot_As == 0.0){BiCGStabConverged=1;break;}                            //   converged ?
    double omega = As_dot_s / As_dot_As;                                        //   omega = As_dot_s / As_dot_As
    if(omega == 0.0){BiCGStabFailed=3;break;}                                   //   stabilization breakdown ???
    if(isinf(omega)){BiCGStabFailed=4;break;}                                   //   stabilization breakdown ???
    add_grids(level,  x_id,  1.0,x_id, omega,__s   );                           //   x_id[] = x_id[] + omega*s[]
    add_grids(level,__r   ,  1.0,__s ,-omega,__As  );                           //   r[]    = s[]    - omega*As[]  (recursively computed / updated residual)
    double norm_of_r = norm(level,__r);                                         //   norm of recursively computed residual (good enough??)
  //if(level->my_rank==0)printf("norm(r)/norm(r0) = %e\n",norm_of_r/norm_of_r0);
    if(norm_of_r == 0.0){BiCGStabConverged=1;break;}                            //
    if(norm_of_r < desired_reduction_in_norm*norm_of_r0){BiCGStabConverged=1;break;}
    #ifdef __DEBUG                                                              //
    residual(level,__temp,x_id,R_id,a,b);                                       //
    double norm_of_residual = norm(level,__temp);                               //
    if(level->my_rank==0)printf("j=%8d, norm=%12.6e, norm_inital=%12.6e, reduction=%e\n",j,norm_of_residual,norm_of_r0,norm_of_residual/norm_of_r0);   //
    #endif                                                                      //
    double r_dot_r0_new = dot(level,__r,__r0);                                  //   r_dot_r0_new = dot(r,r0)
    if(r_dot_r0_new == 0.0){BiCGStabFailed=5;break;}                            //   Lanczos breakdown ???
    double beta = (r_dot_r0_new/r_dot_r0) * (alpha/omega);                      //   beta = (r_dot_r0_new/r_dot_r0) * (alpha/omega)
    if(isinf(beta)){BiCGStabFailed=6;break;}                                    //   ???
    add_grids(level,__temp,1.0,__p,-omega,__Ap  );                              //   __temp =         (p[]-omega*Ap[])
    add_grids(level,__p   ,1.0,__r,  beta,__temp);                              //   p[] = r[] + beta*(p[]-omega*Ap[])
    r_dot_r0 = r_dot_r0_new;                                                    //   r_dot_r0 = r_dot_r0_new   (save old r_dot_r0)
  }                                                                             // }
    #ifdef DiagonallyPrecondition                                               //
    mul_grids(level,x_id,1.0,__Dinv,x_id);                                      //   x_id[] = Dinv[]*x_id[] // i.e. x = D^{-1}x'
    #endif                                                                      //
  #ifdef __DEBUG
  if(BiCGStabFailed)if(level->my_rank==0)printf("BiCGStab Failed... error = %d\n",BiCGStabFailed);
  #endif
}
