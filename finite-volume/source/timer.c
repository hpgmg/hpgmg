//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifdef _OPENMP
#include "timer.omp.c"
#elif USE_MPI
#include "timer.mpi.c"
#else
#error You need to include a custom timer routine
#endif
