//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include <omp.h>
uint64_t CycleTime(){
  return((uint64_t)(1e9*omp_get_wtime())); // convert DP time in seconds to 64b integer nanosecond counter...
}
