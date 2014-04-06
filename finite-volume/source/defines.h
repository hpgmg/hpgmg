//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#define __DEFAULT_BOTTOM_NORM 1e-3
//------------------------------------------------------------------------------------------------------------------------------
#define __RESTRICT_CELL   0
#define __RESTRICT_FACE_I 1
#define __RESTRICT_FACE_J 2
#define __RESTRICT_FACE_K 3
#define __RESTRICT_NODAL  4
//------------------------------------------------------------------------------------------------------------------------------
// Lu = a*alpha[]*u[] - b*divergence beta[]*gradient u[]
//------------------------------------------------------------------------------------------------------------------------------
#define  __temp        0 // 
#define  __u_exact     1 // exact solution used to generate f
#define  __f_minus_Av  2 // cell centered residual (f-Av)
//------------------------------------------------------------------------------------------------------------------------------
#define  __f           3 // original right-hand side (Au=f), cell centered
#define  __u           4 // numerical solution
#define  __alpha       5 // cell centered coefficient
#define  __beta        6 // cell centered coefficient
//------------------------------------------------------------------------------------------------------------------
#define  __beta_i      7 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define  __beta_j      8 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define  __beta_k      9 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
#define  __Dinv       10 // cell centered relaxation parameter (e.g. inverse of the diagonal)
#define  __L1inv      11 // cell centered relaxation parameter (e.g. inverse of the L1 norm of each row)
#define  __valid      12
//------------------------------------------------------------------------------------------------------------------
#define  __Components 13 // total number of grids and the starting location for any auxillary bottom solver grids
