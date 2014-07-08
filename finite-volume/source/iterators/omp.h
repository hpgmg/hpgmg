//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#define FOR_ALL_BOXES(level,box)			\
  _Pragma("omp parallel for private(box) if(level->concurrent_boxes>1) num_threads(level->concurrent_boxes)")\
  for(box=0;box<level->num_my_boxes;box++){		\

#define FOR_ALL_BOXES_MAX(level,box,rv)			\
  _Pragma("omp parallel for private(box) if(level->concurrent_boxes>1) num_threads(level->concurrent_boxes) reduction(max:rv) schedule(static)")\
  for(box=0;box<level->num_my_boxes;box++){		\

#define FOR_ALL_BOXES_SUM(level,box,rv)			\
  _Pragma("omp parallel for private(box) if(level->concurrent_boxes>1) num_threads(level->concurrent_boxes) reduction(+:rv) schedule(static)")\
  for(box=0;box<level->num_my_boxes;box++){		\

#define END_FOR_ALL_BOXES }


#define FOR_ALL_ELEMENTS(level,box,i,j,k,ghostsToOperateOn) 		\
  int dim = level->box_dim;						\
  _Pragma("omp parallel for private(i,j,k) if(level->threads_per_box>1) num_threads(level->threads_per_box) collapse(2)")\
  for(k=-ghostsToOperateOn;k<dim+ghostsToOperateOn;k++){		\
  for(j=-ghostsToOperateOn;j<dim+ghostsToOperateOn;j++){		\
  for(i=-ghostsToOperateOn;i<dim+ghostsToOperateOn;i++){		\

#define FOR_ALL_ELEMENTS_MAX(level,box,i,j,k,ghostsToOperateOn,rv) 	\
  int dim = level->box_dim;						\
  _Pragma("omp parallel for private(i,j,k) if(level->threads_per_box>1) num_threads(level->threads_per_box) collapse(2) reduction(max:rv) schedule(static)")\
  for(k=-ghostsToOperateOn;k<dim+ghostsToOperateOn;k++){		\
  for(j=-ghostsToOperateOn;j<dim+ghostsToOperateOn;j++){		\
  for(i=-ghostsToOperateOn;i<dim+ghostsToOperateOn;i++){		\

#define FOR_ALL_ELEMENTS_SUM(level,box,i,j,k,ghostsToOperateOn,rv) 	\
  int dim = level->box_dim;						\
  _Pragma("omp parallel for private(i,j,k) if(level->threads_per_box>1) num_threads(level->threads_per_box) collapse(2) reduction(+:rv) schedule(static)")\
  for(k=-ghostsToOperateOn;k<dim+ghostsToOperateOn;k++){		\
  for(j=-ghostsToOperateOn;j<dim+ghostsToOperateOn;j++){		\
  for(i=-ghostsToOperateOn;i<dim+ghostsToOperateOn;i++){		\

#define END_FOR_ALL_ELEMENTS }}}

/*
FOR_REDBLACK_ELEMENTS(level,box,i,j,k,ghostsToOperateOn,color)
*/
