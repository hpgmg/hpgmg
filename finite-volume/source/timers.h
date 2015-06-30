//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifndef TIMER_H
#define TIMER_H

  #include<stdint.h>

  #ifdef _OPENMP
    #include <omp.h>
    #define CycleTime() ((uint64_t)(1e9*omp_get_wtime()))

  #elif USE_MPI
    #include <mpi.h>
    #define CycleTime() ((uint64_t)(1e9*MPI_Wtime()))

  #else
    // user must provide a function CycleTime and include it in timers.c
    // if calibration is necesary, then the user must #define CALIBRATE_TIMER
    uint64_t CycleTime();
  #endif

#endif
