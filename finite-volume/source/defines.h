//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
// Lu = a*alpha[]*u[] - b*divergence( beta[]*gradient(u[]) )
//------------------------------------------------------------------------------------------------------------------------------
#ifndef DEFINES_H
#define DEFINES_H
//------------------------------------------------------------------------------------------------------------------------------

#if defined(USE_HELMHOLTZ) || defined(USE_L1JACOBI)
// for those using HPGMG for research...
#warning OVERRIDING BENCHMARK DEFAULTS AND CREATING ALPHA AND L1INV
#define  VECTOR_TEMP         0 // 
#define  VECTOR_U            1 // numerical solution
#define  VECTOR_F            2 // original right-hand side (Au=f), cell centered
#define  VECTOR_E            3 // error used in residual correction FMG
#define  VECTOR_R            4 // cell centered residual (f-Av)
#define  VECTOR_DINV         5 // cell centered relaxation parameter (e.g. inverse of the diagonal)
#define  VECTOR_BETA_I       6 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define  VECTOR_BETA_J       7 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define  VECTOR_BETA_K       8 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
#define  VECTOR_ALPHA        9 // cell centered coefficient used only in Helmholtz
#define  VECTOR_L1INV       10 // cell centered relaxation parameter (e.g. inverse of the L1 norm of each row)
#define VECTORS_RESERVED    11 // total number of vectors and the starting location for any auxillary bottom solver vectors
#else 
// Default HPGMG benchmark configuration (saves memory)
#define  VECTOR_TEMP         0 // 
#define  VECTOR_U            1 // numerical solution
#define  VECTOR_F            2 // original right-hand side (Au=f), cell centered
#define  VECTOR_E            3 // error used in residual correction FMG
#define  VECTOR_R            4 // cell centered residual (f-Av)
#define  VECTOR_DINV         5 // cell centered relaxation parameter (e.g. inverse of the diagonal)
#define  VECTOR_BETA_I       6 // face centered coefficient (n.b. element 0 is the left face of the ghost zone element)
#define  VECTOR_BETA_J       7 // face centered coefficient (n.b. element 0 is the back face of the ghost zone element)
#define  VECTOR_BETA_K       8 // face centered coefficient (n.b. element 0 is the bottom face of the ghost zone element)
#define VECTORS_RESERVED     9 // total number of vectors and the starting location for any auxillary bottom solver vectors
#endif

//------------------------------------------------------------------------------------------------------------------------------
#endif
