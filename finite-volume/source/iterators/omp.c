//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <omp.h>
//------------------------------------------------------------------------------------------------------------------------------
void defaultThreadingForLevel(level_type *level){
  int omp_threads = 1;
  int omp_nested  = 0;
  #pragma omp parallel 
  {
    #pragma omp master
    {
      omp_threads = omp_get_num_threads();
      omp_nested  = omp_get_nested();
    }
  }
  // set default parameters for threading...
  level->threads_per_box  = omp_threads;
  level->concurrent_boxes = 1;
}
//------------------------------------------------------------------------------------------------------------------------------
void tuneThreadingForLevel(level_type *level){
  int omp_threads = 1;
  int omp_nested  = 0;
  #pragma omp parallel 
  {
    #pragma omp master
    {
      omp_threads = omp_get_num_threads();
      omp_nested  = omp_get_nested();
    }
  }
  // inspect omp_nested, omp_num_threads, the number of boxes, and the box size, and choose the optimal varlues for
  // threads_per_box
  // concurrent_boxes
}
//------------------------------------------------------------------------------------------------------------------------------

