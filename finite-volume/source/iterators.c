//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include "iterators.h"
#include "level.h"
//------------------------------------------------------------------------------------------------------------------------------
#ifdef USE_OMP_WORKSHARE
#include "./iterators/omp.c"
#elif  USE_OMP_TASK
#include "./iterators/omptask.c"
#else // default case...
void defaultThreadingForLevel(level_type *level){level->threads_per_box=1;level->concurrent_boxes=1;}
void tuneThreadingForLevel(level_type *level){}
#endif
//------------------------------------------------------------------------------------------------------------------------------
