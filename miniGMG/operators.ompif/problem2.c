//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include "../timer.h"
//------------------------------------------------------------------------------------------------------------------------------
void initialize_problem(domain_type *domain, int level, double hLevel, double a, double b){
  double NPi = 2.0*M_PI;
  double Bmin =  1.0;
  double Bmax = 10.0;
  double c2 = (Bmax-Bmin)/2;
  double c1 = (Bmax+Bmin)/2;
  double c3=10.0; // how sharply (B)eta transitions
  double c4 = -5.0/0.25;

  int box;
  for(box=0;box<domain->subdomains_per_rank;box++){
    memset(domain->subdomains[box].levels[level].grids[__u_exact],0,domain->subdomains[box].levels[level].volume*sizeof(double));
    memset(domain->subdomains[box].levels[level].grids[__f      ],0,domain->subdomains[box].levels[level].volume*sizeof(double));
    int i,j,k;
    #pragma omp parallel for private(k,j,i) collapse(2)
    for(k=0;k<domain->subdomains[box].levels[level].dim.k;k++){
    for(j=0;j<domain->subdomains[box].levels[level].dim.j;j++){
    for(i=0;i<domain->subdomains[box].levels[level].dim.i;i++){
      double x = hLevel*((double)(i+domain->subdomains[box].levels[level].low.i)+0.5);
      double y = hLevel*((double)(j+domain->subdomains[box].levels[level].low.j)+0.5);
      double z = hLevel*((double)(k+domain->subdomains[box].levels[level].low.k)+0.5);
      int ijk =                                              (i+domain->subdomains[box].levels[level].ghosts)+
                domain->subdomains[box].levels[level].pencil*(j+domain->subdomains[box].levels[level].ghosts)+
                domain->subdomains[box].levels[level].plane *(k+domain->subdomains[box].levels[level].ghosts);
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      double r2   = pow((x-0.50),2) +  pow((y-0.50),2) +  pow((z-0.50),2); // distance from center squared
      double r2x  = 2.0*(x-0.50);
      double r2y  = 2.0*(y-0.50);
      double r2z  = 2.0*(z-0.50);
      double r2xx = 2.0;
      double r2yy = 2.0;
      double r2zz = 2.0;
      double r    = pow(r2,0.5);
      double rx   = 0.5*r2x*pow(r2,-0.5);
      double ry   = 0.5*r2y*pow(r2,-0.5);
      double rz   = 0.5*r2z*pow(r2,-0.5);
      double rxx  = 0.5*r2xx*pow(r2,-0.5) - 0.25*r2x*r2x*pow(r2,-1.5);
      double ryy  = 0.5*r2yy*pow(r2,-0.5) - 0.25*r2y*r2y*pow(r2,-1.5);
      double rzz  = 0.5*r2zz*pow(r2,-0.5) - 0.25*r2z*r2z*pow(r2,-1.5);
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // beta (B) is a bubble...
      double A  = 1.0;
      double B  =           c1+c2*tanh( c3*(r-0.25) );
      double Bx = c2*c3*rx*(1-pow(tanh( c3*(r-0.25) ),2));
      double By = c2*c3*ry*(1-pow(tanh( c3*(r-0.25) ),2));
      double Bz = c2*c3*rz*(1-pow(tanh( c3*(r-0.25) ),2)); 
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      // should be continuous in u, u', and u''
      // v(w) = w^4 - 2w^3 + w^2
      // u(x,y,z) = v(x)v(y)v(z)
      double X   =  1.0*pow(x,4) -  2.0*pow(x,3) + 1.0*pow(x,2) - 1.0/30.0;
      double Y   =  1.0*pow(y,4) -  2.0*pow(y,3) + 1.0*pow(y,2) - 1.0/30.0;
      double Z   =  1.0*pow(z,4) -  2.0*pow(z,3) + 1.0*pow(z,2) - 1.0/30.0;
      double Xx  =  4.0*pow(x,3) -  6.0*pow(x,2) + 2.0*x;
      double Yy  =  4.0*pow(y,3) -  6.0*pow(y,2) + 2.0*y;
      double Zz  =  4.0*pow(z,3) -  6.0*pow(z,2) + 2.0*z;
      double Xxx = 12.0*pow(x,2) - 12.0*x        + 2.0;
      double Yyy = 12.0*pow(y,2) - 12.0*y        + 2.0;
      double Zzz = 12.0*pow(z,2) - 12.0*z        + 2.0;
      double u   = X*Y*Z;
      double ux  = Xx*Y*Z;
      double uy  = X*Yy*Z;
      double uz  = X*Y*Zz;
      double uxx = Xxx*Y*Z;
      double uyy = X*Yyy*Z;
      double uzz = X*Y*Zzz;
      double f = a*A*u - b*( (Bx*ux + By*uy + Bz*uz)  +  B*(uxx + uyy + uzz) );
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      domain->subdomains[box].levels[level].grids[__alpha  ][ijk] = A;
      domain->subdomains[box].levels[level].grids[__beta   ][ijk] = B;
      domain->subdomains[box].levels[level].grids[__u_exact][ijk] = u;
      domain->subdomains[box].levels[level].grids[__f      ][ijk] = f;
    }}}
  }

  double average_value_of_f = mean(domain,level,__f);
  if(domain->rank==0){printf("\n  average value of f = %20.12e\n",average_value_of_f);fflush(stdout);}
  if(a!=0){
  shift_grid(domain,level,__f      ,__f      ,-average_value_of_f);
  shift_grid(domain,level,__u_exact,__u_exact,-average_value_of_f/a);
  }
  // what if a==0 and average_value_of_f != 0 ???
}
//------------------------------------------------------------------------------------------------------------------------------
