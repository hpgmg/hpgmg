#ifndef _pointwise_h
#define _pointwise_h

#include <petscsys.h>
#include "fefas-align.h"

static PetscErrorCode PointwiseJacobianInvert(PetscInt ne,PetscInt Q,const PetscReal w[Q],PetscScalar dx[3][3][Q][ne],PetscScalar wdxdet[Q][ne])
{
  PetscInt i,j,k,e;

  for (i=0; i<Q; i++) {
    PetscScalar a[3][3][ne]_align;
    for (e=0; e<ne; e++) {
      PetscScalar b0,b3,b6,det,idet;
      for (j=0; j<3; j++) {
        for (k=0; k<3; k++) {
          a[j][k][e] = dx[j][k][i][e];
        }
      }
      b0 =  (a[1][1][e]*a[2][2][e] - a[2][1][e]*a[1][2][e]);
      b3 = -(a[1][0][e]*a[2][2][e] - a[2][0][e]*a[1][2][e]);
      b6 =  (a[1][0][e]*a[2][1][e] - a[2][0][e]*a[1][1][e]);
      det = a[0][0][e]*b0 + a[0][1][e]*b3 + a[0][2][e]*b6;
      idet = 1.0 / det;
      dx[0][0][i][e] =  idet*b0;
      dx[0][1][i][e] = -idet*(a[0][1][e]*a[2][2][e] - a[2][1][e]*a[0][2][e]);
      dx[0][2][i][e] =  idet*(a[0][1][e]*a[1][2][e] - a[1][1][e]*a[0][2][e]);
      dx[1][0][i][e] =  idet*b3;
      dx[1][1][i][e] =  idet*(a[0][0][e]*a[2][2][e] - a[2][0][e]*a[0][2][e]);
      dx[1][2][i][e] = -idet*(a[0][0][e]*a[1][2][e] - a[1][0][e]*a[0][2][e]);
      dx[2][0][i][e] =  idet*b6;
      dx[2][1][i][e] = -idet*(a[0][0][e]*a[2][1][e] - a[2][0][e]*a[0][1][e]);
      dx[2][2][i][e] =  idet*(a[0][0][e]*a[1][1][e] - a[1][0][e]*a[0][1][e]);
      wdxdet[i][e] =  det*w[i];
    }
  }
  PetscLogFlops(Q*ne*(14 + 1/* division */ + 27 + 1));
  return 0;
}

#endif
