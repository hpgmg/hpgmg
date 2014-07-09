//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include <mpi.h>
uint64_t CycleTime(){
  return((uint64_t)(1e9*MPI_Wtime())); // convert DP time in seconds to 64b integer nanosecond counter...
}
