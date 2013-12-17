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
#include "cuda_runtime.h"
#include "cuda.h"
//#define cudaMemcpyHostToDevice 1
//------------------------------------------------------------------------------------------------------------------------------
#include <omp.h>
#ifdef __MPI
#include <mpi.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include "defines.h"
#include "box.h"
#include "mg.h"
#include "kernels.h"

//==============================================================================
void initialize_rhs(double * grid_id, int lowi, int lowj, int lowk, int dim, int volume,int ghosts, int pencil, int plane, double h){
  int i,j,k;
  double twoPi = 2.0 * 3.1415926535;
  memset(grid_id,0,volume*sizeof(double)); // zero out the grid and ghost zones
  for(k=0;k<dim;k++){
  for(j=0;j<dim;j++){
  for(i=0;i<dim;i++){
    double x = h*(double)((i+lowi)+0.5);
    double y = h*(double)((j+lowj)+0.5);
    double z = h*(double)((k+lowk)+0.5);
    int ijk = (i+ghosts) + (j+ghosts)*pencil + (k+ghosts)*plane;
    double value = sin(twoPi*x)*sin(twoPi*y)*sin(twoPi*z);
    grid_id[ijk] = value;
  }}}
}

//==============================================================================
void initialize_exact(domain_type *domain, int level, double hLevel, double a, double b){
  zero_grid(domain,level,__alpha  ); 
  zero_grid(domain,level,__beta   ); 
  zero_grid(domain,level,__beta_i ); 
  zero_grid(domain,level,__beta_j ); 
  zero_grid(domain,level,__beta_k ); 
  zero_grid(domain,level,__lambda ); 
  zero_grid(domain,level,__u_exact); 
  zero_grid(domain,level,__f      ); 
  initialize_exact_on_gpu(domain,level,hLevel,a,b);
  double average_value_of_f = mean(domain,level,__f);
  if(a!=0){
  shift_grid(domain,level,__f      ,__f      ,-average_value_of_f);
  shift_grid(domain,level,__u_exact,__u_exact,-average_value_of_f/a);
  }
}


//==============================================================================
int main(int argc, char **argv){
  int MPI_Rank=0;
  int MPI_Tasks=1;
  int OMP_Threads = 1;
/*
  #pragma omp parallel 
  {
    #pragma omp master
    {
      OMP_Threads = omp_get_num_threads();
    }
  }
*/    

  #ifdef __MPI
  #warning Compiling for MPI...
  int MPI_threadingModel          = -1;
//int MPI_threadingModelRequested = MPI_THREAD_SINGLE;
  int MPI_threadingModelRequested = MPI_THREAD_FUNNELED;
//int MPI_threadingModelRequested = MPI_THREAD_MULTIPLE;
  MPI_Init_thread(&argc, &argv, MPI_threadingModelRequested, &MPI_threadingModel);
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_Tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_Rank);

  if(MPI_threadingModel>MPI_threadingModelRequested)MPI_threadingModel=MPI_threadingModelRequested;
  if(MPI_Rank==0){
       if(MPI_threadingModelRequested == MPI_THREAD_MULTIPLE  )printf("Requested MPI_THREAD_MULTIPLE, ");
  else if(MPI_threadingModelRequested == MPI_THREAD_SINGLE    )printf("Requested MPI_THREAD_SINGLE, ");
  else if(MPI_threadingModelRequested == MPI_THREAD_FUNNELED  )printf("Requested MPI_THREAD_FUNNELED, ");
  else if(MPI_threadingModelRequested == MPI_THREAD_SERIALIZED)printf("Requested MPI_THREAD_SERIALIZED, ");
  else if(MPI_threadingModelRequested == MPI_THREAD_MULTIPLE  )printf("Requested MPI_THREAD_MULTIPLE, ");
  else                                                printf("got Unknown MPI_threadingModel (%d)\n",MPI_threadingModel);
       if(MPI_threadingModel == MPI_THREAD_MULTIPLE  )printf("got MPI_THREAD_MULTIPLE\n");
  else if(MPI_threadingModel == MPI_THREAD_SINGLE    )printf("got MPI_THREAD_SINGLE\n");
  else if(MPI_threadingModel == MPI_THREAD_FUNNELED  )printf("got MPI_THREAD_FUNNELED\n");
  else if(MPI_threadingModel == MPI_THREAD_SERIALIZED)printf("got MPI_THREAD_SERIALIZED\n");
  else if(MPI_threadingModel == MPI_THREAD_MULTIPLE  )printf("got MPI_THREAD_MULTIPLE\n");
  else                                                printf("got Unknown MPI_threadingModel (%d)\n",MPI_threadingModel);
  fflush(stdout);  }
  #endif

  if(MPI_Rank==0){
  printf("configuration: %2dx%dx%d, ",TBDIMX,TBDIMY,1);
  #ifdef __LOCALITY_VIA_SHARED
  printf("phi in shared, ");
  #else
  printf("phi in L1,     ");
  #endif
  #ifdef __POINTERS_IN_SHARED
  printf("pointers in shared,  ");
  #else
  printf("pointers in private, ");
  #endif
  #ifdef __PREFER_SHARED
  printf("16K L1 + 48K shared, ");
  #else
  printf("48K L1 + 16K shared, ");
  #endif
  #ifdef VL
  printf("using 1D thread blocks with VL=%d ",VL);
  #endif
  printf("\n");
  }

  int log2_subdomain_dim = 6;
//    log2_subdomain_dim = 7;
//    log2_subdomain_dim = 5;
//    log2_subdomain_dim = 2;
  int subdomains_per_rank_in_i=256 / (1<<log2_subdomain_dim);
  int subdomains_per_rank_in_j=256 / (1<<log2_subdomain_dim);
  int subdomains_per_rank_in_k=256 / (1<<log2_subdomain_dim);
  int ranks_in_i=1;
  int ranks_in_j=1;
  int ranks_in_k=1;

  if(argc==2){
          log2_subdomain_dim=atoi(argv[1]);
          subdomains_per_rank_in_i=256 / (1<<log2_subdomain_dim);
          subdomains_per_rank_in_j=256 / (1<<log2_subdomain_dim);
          subdomains_per_rank_in_k=256 / (1<<log2_subdomain_dim);
  }else if(argc==5){
          log2_subdomain_dim=atoi(argv[1]);
    subdomains_per_rank_in_i=atoi(argv[2]);
    subdomains_per_rank_in_j=atoi(argv[3]);
    subdomains_per_rank_in_k=atoi(argv[4]);
  }else if(argc==8){
          log2_subdomain_dim=atoi(argv[1]);
    subdomains_per_rank_in_i=atoi(argv[2]);
    subdomains_per_rank_in_j=atoi(argv[3]);
    subdomains_per_rank_in_k=atoi(argv[4]);
                  ranks_in_i=atoi(argv[5]);
                  ranks_in_j=atoi(argv[6]);
                  ranks_in_k=atoi(argv[7]);
  }else if(argc!=1){
    if(MPI_Rank==0){printf("usage: ./a.out [log2_subdomain_dim]   [subdomains per rank in i,j,k]  [ranks in i,j,k]\n");}
    #ifdef __MPI
    MPI_Finalize();
    #endif
    exit(0);
  }

/*
  if(log2_subdomain_dim!=6){
    if(MPI_Rank==0){printf("error, log2_subdomain_dim(%d)!=6\n",log2_subdomain_dim);}
    #ifdef __MPI
    MPI_Finalize();
    #endif
    exit(0);
  }
*/

  if(ranks_in_i*ranks_in_j*ranks_in_k != MPI_Tasks){
    if(MPI_Rank==0){printf("error, ranks_in_i*ranks_in_j*ranks_in_k(%d*%d*%d=%d) != MPI_Tasks(%d)\n",ranks_in_i,ranks_in_j,ranks_in_k,ranks_in_i*ranks_in_j*ranks_in_k,MPI_Tasks);}
    #ifdef __MPI
    MPI_Finalize();
    #endif
    exit(0);
  }

  if(MPI_Rank==0)printf("%d MPI Tasks of %d threads\n",MPI_Tasks,OMP_Threads);

  int subdomain_dim_i=1<<log2_subdomain_dim;
  int subdomain_dim_j=1<<log2_subdomain_dim;
  int subdomain_dim_k=1<<log2_subdomain_dim;
  //    dim = 128 64 32 16 8 4
  // levels =   6  5  4  3 2 1
  int log2_coarse_dim = 1; // i.e. coarsen to 2^3
  int levels_in_vcycle=1+log2_subdomain_dim-log2_coarse_dim; // ie 1+log2(fine grid size)-log2(bottom grid size)
//int levels_in_vcycle=(log2_subdomain_dim+1)-2; // ie -log2(bottom size)

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int box,grid,level,i,j,k,n;
  int ghosts = 1; // MUST BE 1 for Sam's version of the GPU code...

  domain_type domain;
  create_domain(&domain,subdomain_dim_i,subdomain_dim_j,subdomain_dim_k,
                        subdomains_per_rank_in_i,subdomains_per_rank_in_j,subdomains_per_rank_in_k,
                        ranks_in_i,ranks_in_j,ranks_in_k,
                        MPI_Rank,__NumGrids,ghosts,levels_in_vcycle);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  double h0=1.0/((double)(domain.dim.i));
  double ConvergenceCriteria = 1e-15;
  if(MPI_Rank==0){printf("calculating convergence criteria... %e\n",ConvergenceCriteria);fflush(stdout);}
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  double  a=0.9;
  double  b=0.9;

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #if 0
  // define __alpha, __beta*, etc...
  if(MPI_Rank==0){printf("initializing coefficients...");fflush(stdout);}
  initialize_grid_to_scalar(&domain,0,__alpha ,h0,1.0);
  initialize_grid_to_scalar(&domain,0,__beta  ,h0,1.0);
  if(MPI_Rank==0){printf("done\n");fflush(stdout);}
 
  // define RHS
  if(MPI_Rank==0){printf("initializing RHS...");fflush(stdout);}
  double * tempo = (double*)malloc( domain.subdomains[0].levels[0].volume*sizeof(double));  
  for(box=0;box<domain.numsubdomains;box++){ // initialize RHS locally, then copy to device...
    initialize_rhs(&tempo[0],  domain.subdomains[box].levels[0].low.i, domain.subdomains[box].levels[0].low.j, domain.subdomains[box].levels[0].low.k, 
                    domain.subdomains[box].levels[0].dim.i,  domain.subdomains[box].levels[0].volume,ghosts, 
                    domain.subdomains[box].levels[0].pencil, domain.subdomains[box].levels[0].plane,  h0);
    cudaMemcpy(domain.subdomains[box].levels[0].grids[__f],tempo, domain.subdomains[box].levels[0].volume*sizeof(double),cudaMemcpyHostToDevice);
  }free(tempo);
  if(MPI_Rank==0){printf("done\n");fflush(stdout);}
  #else
  if(MPI_Rank==0){printf("initializing problem...");fflush(stdout);}
  initialize_exact(&domain,0,h0,a,b);
  if(MPI_Rank==0){printf("done\n");fflush(stdout);}
  #endif
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // MGBuild
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  uint64_t _timeStartMGBuild;
  _timeStartMGBuild = CycleTime();

  // initialize timers...
  for(level=0;level<10;level++){
  domain.cycles.communication[level] = 0;
  domain.cycles.recv[level]          = 0;
  domain.cycles.send[level]          = 0;
  domain.cycles.wait[level]          = 0;
  domain.cycles.collectives[level]   = 0;
  domain.cycles.PCIe[level]          = 0;
  }
  domain.cycles.build                = 0;
  domain.cycles.vcycles              = 0;
  domain.cycles.MGSolve              = 0;

  // form all restrictions of alpha[] for all boxes...
  if(MPI_Rank==0){printf("restricting alpha...");fflush(stdout);}
  for(level=0;level<domain.numLevels-1;level++){      restriction(&domain,level,__alpha,__alpha);}
  for(level=0;level<domain.numLevels  ;level++){exchange_boundary(&domain,level,__alpha ,1,1,1);} // FIX, only necessary if CA version
  if(MPI_Rank==0){printf("done\n");fflush(stdout);}

  // form all restrictions of beta_*[] for all boxes...
  if(MPI_Rank==0){printf("restricting betas...");fflush(stdout);}
     exchange_boundary(&domain,0,__beta,1,1,1);
  project_cell_to_face(&domain,0,__beta,__beta_i,0);
  project_cell_to_face(&domain,0,__beta,__beta_j,1);
  project_cell_to_face(&domain,0,__beta,__beta_k,2);
  for(level=0;level<domain.numLevels;level++){
                         exchange_boundary(&domain,level,__beta_i,1,1,1); // lambda needs high betas
                         exchange_boundary(&domain,level,__beta_j,1,1,1);
                         exchange_boundary(&domain,level,__beta_k,1,1,1);
    if(level<domain.numLevels-1)restriction_betas(&domain,level);
  }
  if(MPI_Rank==0){printf("done\n");fflush(stdout);}

  // form all restrictions of lambda[] for all boxes...
  if(MPI_Rank==0){printf("building lambda...");fflush(stdout);}
              for(level=0;level<domain.numLevels;level++){   rebuild_lambda(&domain,level,a,b,h0 * (double)(1<<level));}
  if(ghosts>1)for(level=0;level<domain.numLevels;level++){exchange_boundary(&domain,level,__lambda,1,1,1);} // FIX, only necessary if CA version
  if(MPI_Rank==0){printf("done\n");fflush(stdout);}

  domain.cycles.build += (uint64_t)(CycleTime()-_timeStartMGBuild);

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #warning Write CUDA generate problem
  #warning Write CUDA project cell centered to face centered
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int s,sMax=4;
  #ifdef __MPI
  sMax=8;
  #endif
  for(s=0;s<sMax;s++){MGSolve(&domain,__u,__f,a,b,h0,ConvergenceCriteria);}
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  print_timing(&domain);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // calculate error...
  double h3 = h0*h0*h0;
                 add_grids(&domain,0,__temp,1.0,__u_exact,-1.0,__u);       // __temp = __u_exact - __u
  double   max =      norm(&domain,0,__temp);                              // max norm of error function
  double error = sqrt( dot(&domain,0,__temp,__temp)*h3);                   // normalized L2 error ?
  if(MPI_Rank==0){printf("Error test: h = %e, max = %e\n",h0,max);}
  if(MPI_Rank==0){printf("Error test: h = %e, L2  = %e\n",h0,error);}
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  destroy_domain(&domain);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef __MPI
  MPI_Finalize();
  #endif
  return(0);
}
