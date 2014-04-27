#include "fefas.h"

#ifdef __bgq__
#  include <spi/include/kernel/memory.h>
#endif

PetscErrorCode MemoryGetUsage(double *heapused,double *heapavail) {
  PetscFunctionBegin;
  *heapused = -1;
  *heapavail = -1;
#ifdef __bgq__
  {
    uint64_t heap,avail;
    Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP,&heap);
    Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL,&avail);
    *heapused = (double)heap;
    *heapavail = (double)avail;
  }
#else
  {
    PetscErrorCode ierr;
    ierr = PetscMemoryGetCurrentUsage(heapused);CHKERRQ(ierr);
  }
#endif
  PetscFunctionReturn(0);
}
