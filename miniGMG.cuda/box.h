//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  struct {int i, j, k;}low;               // global coordinates of the first (non-ghost) element of subdomain
  struct {int i, j, k;}dim;               // dimensions of this box's core (not counting ghost zone)
  struct {int i, j, k;}dim_with_ghosts;   // dimensions of this box's core (not counting ghost zone)
  int ghosts;                             // ghost zone depth
  int pencil,plane,volume;                // useful for offsets
  int                       bufsizes[27]; // = sizes of extracted surfaces and ghost zones (pointer to array of 27 elements)
  double * __restrict__ surface_bufs[27]; // = extracted surface (rhs on the way down, correction on the way up)  (GPU pointers)
  double * __restrict__   ghost_bufs[27]; // = incoming ghost zone (rhs on the way down, correction on the way up)  (GPU pointers)

  int numGrids;
  double ** __restrict__ grids; // CPU pointers to an array of GPU pointers to grid data
} box_type;
//------------------------------------------------------------------------------------------------------------------------------
void destroy_box(box_type *box);
 int create_box(box_type *box, int numGrids, int low_i, int low_j, int low_k, int dim_i, int dim_j, int dim_k, int ghosts);
