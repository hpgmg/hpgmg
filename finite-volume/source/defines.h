//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// Lu = a*alpha[]*u[] - b*divergence beta[]*gradient u[]
//------------------------------------------------------------------------------------------------------------------------------
#define STENCIL_TEMP         0 // 
#define STENCIL_UTRUE        1 // exact solution used to generate f
#define STENCIL_F_MINUS_AV   2 // cell centered residual (f-Av)
//------------------------------------------------------------------------------------------------------------------------------
#define STENCIL_F            3 // original right-hand side (Au=f), cell centered
#define STENCIL_U            4 // numerical solution
#define STENCIL_ALPHA        5 // cell centered coefficient
#define STENCIL_BETA_I       6 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define STENCIL_BETA_J       7 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define STENCIL_BETA_K       8 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
//------------------------------------------------------------------------------------------------------------------
#define STENCIL_DINV         9 // cell centered relaxation parameter (e.g. inverse of the diagonal)
#define STENCIL_L1INV       10 // cell centered relaxation parameter (e.g. inverse of the L1 norm of each row)
#define STENCIL_VALID       11 // cell centered array noting which cells are actually present
//------------------------------------------------------------------------------------------------------------------
#define COMPONENTS_RESERVED 12 // total number of grids and the starting location for any auxillary bottom solver grids
