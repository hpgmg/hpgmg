//------------------------------------------------------------------------------------------------------------------------------
// Copyright Notice 
//------------------------------------------------------------------------------------------------------------------------------
// HPGMG, Copyright (c) 2014, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to
// receipt of any required approvals from the U.S. Dept. of Energy).  All
// rights reserved.
// 
// If you have questions about your rights to use or distribute this
// software, please contact Berkeley Lab's Technology Transfer Department
// at  TTD@lbl.gov.
// 
// NOTICE.  This software is owned by the U.S. Department of Energy.  As
// such, the U.S. Government has been granted for itself and others
// acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
// license in the Software to reproduce, prepare derivative works, and
// perform publicly and display publicly.  Beginning five (5) years after
// the date permission to assert copyright is obtained from the U.S.
// Department of Energy, and subject to any subsequent five (5) year
// renewals, the U.S. Government is granted for itself and others acting
// on its behalf a paid-up, nonexclusive, irrevocable, worldwide license
// in the Software to reproduce, prepare derivative works, distribute
// copies to the public, perform publicly and display publicly, and to
// permit others to do so.
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
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include "defines.h"
#include "level.h"
#include "mg.h"
#include "operators.h"
#include "solvers.h"
//------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv){
  int my_rank=0;
  int num_tasks=1;
  int OMP_Threads = 1;
  int OMP_Nested = 0;

  #ifdef _OPENMP
  #pragma omp parallel 
  {
    #pragma omp master
    {
      OMP_Threads = omp_get_num_threads();
      OMP_Nested  = omp_get_nested();
    }
  }
  #endif
    

  #ifdef USE_MPI
  int    actual_threading_model = -1;
  int requested_threading_model = -1;
      requested_threading_model = MPI_THREAD_SINGLE;
    //requested_threading_model = MPI_THREAD_FUNNELED;
    //requested_threading_model = MPI_THREAD_SERIALIZED;
    //requested_threading_model = MPI_THREAD_MULTIPLE;
  //MPI_Init(&argc, &argv);
  #ifdef _OPENMP
      requested_threading_model = MPI_THREAD_FUNNELED;
    //requested_threading_model = MPI_THREAD_SERIALIZED;
    //requested_threading_model = MPI_THREAD_MULTIPLE;
  //MPI_Init_thread(&argc, &argv, requested_threading_model, &actual_threading_model);
  #endif
  MPI_Init_thread(&argc, &argv, requested_threading_model, &actual_threading_model);
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//if(actual_threading_model>requested_threading_model)actual_threading_model=requested_threading_model;
  if(my_rank==0){
       if(requested_threading_model == MPI_THREAD_MULTIPLE  )printf("Requested MPI_THREAD_MULTIPLE, ");
  else if(requested_threading_model == MPI_THREAD_SINGLE    )printf("Requested MPI_THREAD_SINGLE, ");
  else if(requested_threading_model == MPI_THREAD_FUNNELED  )printf("Requested MPI_THREAD_FUNNELED, ");
  else if(requested_threading_model == MPI_THREAD_SERIALIZED)printf("Requested MPI_THREAD_SERIALIZED, ");
  else if(requested_threading_model == MPI_THREAD_MULTIPLE  )printf("Requested MPI_THREAD_MULTIPLE, ");
  else                                                       printf("Requested Unknown MPI Threading Model (%d), ",requested_threading_model);
       if(actual_threading_model    == MPI_THREAD_MULTIPLE  )printf("got MPI_THREAD_MULTIPLE\n");
  else if(actual_threading_model    == MPI_THREAD_SINGLE    )printf("got MPI_THREAD_SINGLE\n");
  else if(actual_threading_model    == MPI_THREAD_FUNNELED  )printf("got MPI_THREAD_FUNNELED\n");
  else if(actual_threading_model    == MPI_THREAD_SERIALIZED)printf("got MPI_THREAD_SERIALIZED\n");
  else if(actual_threading_model    == MPI_THREAD_MULTIPLE  )printf("got MPI_THREAD_MULTIPLE\n");
  else                                                       printf("got Unknown MPI Threading Model (%d)\n",actual_threading_model);
  }
  #ifdef USE_HPM // IBM HPM counters for BGQ...
  HPM_Init();
  #endif
  #endif // USE_MPI


  int log2_box_dim = 6;
  int target_boxes_per_rank = 1;

  if(argc==3){
           log2_box_dim=atoi(argv[1]);
     target_boxes_per_rank=atoi(argv[2]);
  }else{
    if(my_rank==0){printf("usage: ./a.out  [log2_box_dim]  [target_boxes_per_rank]\n");}
    #ifdef USE_MPI
    MPI_Finalize();
    #endif
    exit(0);
  }

  if(log2_box_dim<4){
    if(my_rank==0){printf("log2_box_dim must be at least 4\n");}
    #ifdef USE_MPI
    MPI_Finalize();
    #endif
    exit(0);
  }

  if(target_boxes_per_rank<1){
    if(my_rank==0){printf("target_boxes_per_rank must be at least 1\n");}
    #ifdef USE_MPI
    MPI_Finalize();
    #endif
    exit(0);
  }

  if(my_rank==0){
    if(OMP_Nested)printf("%d MPI Tasks of %d threads (OMP_NESTED=TRUE)\n\n" ,num_tasks,OMP_Threads);
             else printf("%d MPI Tasks of %d threads (OMP_NESTED=FALSE)\n\n",num_tasks,OMP_Threads);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  // calculate the problem size...
  int box_dim=1<<log2_box_dim;
  int target_boxes = target_boxes_per_rank*num_tasks;
  int boxes_in_i = 1000; // FIX, int64_t?  could we really have >2e9 boxes?
  int total_boxes = boxes_in_i*boxes_in_i*boxes_in_i;
  while(total_boxes>target_boxes){
    boxes_in_i--;
    total_boxes = boxes_in_i*boxes_in_i*boxes_in_i;
  }

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // create the fine level...
  int ghosts=1;
  level_type fine_grid;
  //create_level(&fine_grid,boxes_in_i,box_dim,ghosts,VECTORS_RESERVED,BC_PERIODIC ,my_rank,num_tasks);double h0=1.0/( (double)boxes_in_i*(double)box_dim );double a=2.0;double b=1.0; // Helmholtz w/Periodic
  //create_level(&fine_grid,boxes_in_i,box_dim,ghosts,VECTORS_RESERVED,BC_PERIODIC ,my_rank,num_tasks);double h0=1.0/( (double)boxes_in_i*(double)box_dim );double a=0.0;double b=1.0; //   Poisson w/Periodic
  //create_level(&fine_grid,boxes_in_i,box_dim,ghosts,VECTORS_RESERVED,BC_DIRICHLET,my_rank,num_tasks);double h0=1.0/( (double)boxes_in_i*(double)box_dim );double a=2.0;double b=1.0; // Helmholtz w/Dirichlet
    create_level(&fine_grid,boxes_in_i,box_dim,ghosts,VECTORS_RESERVED,BC_DIRICHLET,my_rank,num_tasks);double h0=1.0/( (double)boxes_in_i*(double)box_dim );double a=0.0;double b=1.0; //   Poisson w/Dirichlet
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  initialize_problem(&fine_grid,h0,a,b);
  rebuild_operator(&fine_grid,NULL,a,b); // i.e. calculate Dinv and lambda_max
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  mg_type all_grids;
  int minCoarseDim = 1;
  MGBuild(&all_grids,&fine_grid,a,b,minCoarseDim);
  fflush(stdout);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  int doTiming;
  for(doTiming=0;doTiming<=1;doTiming++){ // first pass warms up, second pass times
    MGResetTimers(&all_grids);
    #ifdef USE_HPM // IBM performance counters for BGQ...
    if(doTiming)HPM_Start("FMGSolve()");
    #endif
       int numSolves = 0; // solves completed
       int minSolves = 5; // do at least minSolves MGSolves
    #ifdef USE_MPI
    double minTime   = 10.0; // minimum time in seconds that the benchmark should run
    double startTime = MPI_Wtime();
    while( (numSolves<minSolves) || ((MPI_Wtime()-startTime)<minTime) ){
    #else
    while( (numSolves<minSolves)                                      ){
    #endif
      zero_vector(all_grids.levels[0],VECTOR_U);
      #ifdef USE_FCYCLES
      FMGSolve(&all_grids,VECTOR_U,VECTOR_F,a,b,1e-15);
      #else
       MGSolve(&all_grids,VECTOR_U,VECTOR_F,a,b,1e-15);
      #endif
      numSolves++;
    }
    #ifdef USE_HPM // IBM performance counters for BGQ...
    if(doTiming)HPM_Stop("FMGSolve()");
    #endif
  }
  MGPrintTiming(&all_grids); // don't include the error check in the timing results
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(my_rank==0){printf("calculating error...\n");}
  double fine_error = error(&fine_grid,VECTOR_U,VECTOR_UTRUE);if(my_rank==0){printf(" h = %22.15e  ||error|| = %22.15e\n\n",h0,fine_error);fflush(stdout);}
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // MGDestroy()
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  #ifdef USE_MPI
  #ifdef USE_HPM // IBM performance counters for BGQ...
  HPM_Print();
  #endif
  MPI_Finalize();
  #endif
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  return(0);
}
