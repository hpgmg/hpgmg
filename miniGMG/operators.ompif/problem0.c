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
      // constant coefficient
      double A  = 1.0;
      double B  = 1.0;
      double Bx = 0.0;
      double By = 0.0;
      double Bz = 0.0;
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      //double u   =          sin(NPi*x)*sin(NPi*y)*sin(NPi*z);
      //double ux  =      NPi*cos(NPi*x)*sin(NPi*y)*sin(NPi*z);
      //double uy  =      NPi*sin(NPi*x)*cos(NPi*y)*sin(NPi*z);
      //double uz  =      NPi*sin(NPi*x)*sin(NPi*y)*cos(NPi*z);
      //double uxx = -NPi*NPi*sin(NPi*x)*sin(NPi*y)*sin(NPi*z);
      //double uyy = -NPi*NPi*sin(NPi*x)*sin(NPi*y)*sin(NPi*z);
      //double uzz = -NPi*NPi*sin(NPi*x)*sin(NPi*y)*sin(NPi*z);
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
      //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
