//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// Lu = a*alpha[]*u[] - b*divergence beta[]*gradient u[]
//------------------------------------------------------------------------------------------------------------------------------
#ifndef DEFINES_H
#define DEFINES_H
//------------------------------------------------------------------------------------------------------------------------------
#define  VECTOR_TEMP         0 // 
#define  VECTOR_UTRUE        1 // exact solution used to generate f
#define  VECTOR_F_MINUS_AV   2 // cell centered residual (f-Av)
//------------------------------------------------------------------------------------------------------------------------------
#define  VECTOR_F            3 // original right-hand side (Au=f), cell centered
#define  VECTOR_U            4 // numerical solution
#define  VECTOR_ALPHA        5 // cell centered coefficient
#define  VECTOR_BETA_I       6 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define  VECTOR_BETA_J       7 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define  VECTOR_BETA_K       8 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
//------------------------------------------------------------------------------------------------------------------
#define  VECTOR_DINV         9 // cell centered relaxation parameter (e.g. inverse of the diagonal)
#define  VECTOR_L1INV       10 // cell centered relaxation parameter (e.g. inverse of the L1 norm of each row)
#define  VECTOR_VALID       11 // cell centered array noting which cells are actually present
//------------------------------------------------------------------------------------------------------------------
#define VECTORS_RESERVED    12 // total number of grids and the starting location for any auxillary bottom solver grids
//------------------------------------------------------------------------------------------------------------------------------
#endif
