//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifdef _OPENMP
// CycleTime in OpenMP is now defined as a preprocessor macro
//#include "./timers/omp.c"
#elif USE_MPI
// CycleTime in MPI is now defined as a preprocessor macro
//#include "./timers/mpi.c"
#else
#error You need to include a custom timer routine
#endif
