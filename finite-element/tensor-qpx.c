#include "tensorimpl.h"

#ifdef __bgq__

#define NE 4

#define Pragma_(a) _Pragma(#a)
#define Pragma(a) Pragma_(a)

// xlc refuses to inline this function, so we use an evil macro.
// static inline PetscErrorCode TensorContract_QPX(PetscInt dof,PetscInt P,PetscInt Q,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[])

#define TensorContract_QPX(dof_,P_,Q_,Rf,Sf,Tf,tmode,xx,yy) do {	\
    PetscInt dof = dof_,P = P_,Q = Q_;					\
    /* PetscFunctionBegin; */						\
    if (tmode == TENSOR_TRANSPOSE) {PetscInt tmp = Q; Q = P; P = tmp;}	\
    {									\
      PetscReal R[Q][P],S[Q][P],T[Q][P];				\
      PetscScalar (*x)[P*P*P][NE] = (PetscScalar(*)[P*P*P][NE])xx; /* should be const, but xlc vector intrinsics don't understand const */ \
      PetscScalar (*y)[P*P*P][NE] = (PetscScalar(*)[Q*Q*Q][NE])yy;	\
      PetscScalar u[dof][Q*P*P][NE]_align,v[dof][Q*Q*P][NE]_align;	\
									\
      for (PetscInt i=0; i<Q; i++) {					\
	for (PetscInt j=0; j<P; j++) {					\
	  R[i][j] = tmode == TENSOR_EVAL ? Rf[i*P+j] : Rf[j*Q+i];	\
	  S[i][j] = tmode == TENSOR_EVAL ? Sf[i*P+j] : Sf[j*Q+i];	\
	  T[i][j] = tmode == TENSOR_EVAL ? Tf[i*P+j] : Tf[j*Q+i];	\
	}								\
      }									\
									\
      /* u[l,a,j,k] = R[a,i] x[l,i,j,k] */				\
      for (PetscInt l=0; l<dof; l++) {					\
	for (PetscInt a=0; a<Q; a++) {					\
	  vector4double r[P];						\
	  for (PetscInt i=0; i<P; i++) r[i] = vec_lds(0,&R[a][i]);	\
	  Pragma(unrollandfuse(P_))	/* P==Q required for this */	\
	    for (PetscInt jk=0; jk<P*P; jk++) {				\
	      vector4double u_lajk = {0.};				\
	      for (PetscInt i=0; i<P; i++) {				\
		u_lajk = vec_madd(r[i],vec_lda(0,x[l][i*P*P+jk]),u_lajk); \
	      }								\
	      vec_sta(u_lajk,0,u[l][a*P*P+jk]);				\
	    }								\
	}								\
      }									\
									\
      /* v[l,a,b,k] = S[b,j] u[l,a,j,k] */				\
      for (PetscInt l=0; l<dof; l++) {					\
	for (PetscInt b=0; b<Q; b++) {					\
	  vector4double s[P];						\
	  for (int j=0; j<P; j++) s[j] = vec_lds(0,&S[b][j]);		\
	  for (PetscInt a=0; a<Q; a++) {				\
	    Pragma(unrollandfuse(P_))					\
	      for (PetscInt k=0; k<P; k++) {				\
		vector4double v_labk = {0.};				\
		for (PetscInt j=0; j<P; j++) {				\
		  v_labk = vec_madd(s[j],vec_lda(0,u[l][(a*P+j)*P+k]),v_labk); \
		}							\
		vec_sta(v_labk,0,v[l][(a*Q+b)*P+k]);			\
	      }								\
	  }								\
	}								\
      }									\
									\
      /* y[l,a,b,c] = T[c,k] v[l,a,b,k] */				\
      for (PetscInt l=0; l<dof; l++) {					\
	for (PetscInt c=0; c<Q; c++) {					\
	  vector4double t[P];						\
	  for (int k=0; k<P; k++) t[k] = vec_lds(0,&T[c][k]);		\
	  Pragma(unrollandfuse(Q_))					\
	    for (PetscInt ab=0; ab<Q*Q; ab++) {				\
	      vector4double y_labc = vec_lda(0,y[l][ab*Q+c]);		\
	      for (PetscInt k=0; k<P; k++) {				\
		y_labc = vec_madd(t[k],vec_lda(0,v[l][ab*P+k]),y_labc);	\
	      }								\
	      vec_sta(y_labc,0,y[l][ab*Q+c]);				\
	    }								\
	}								\
      }									\
      PetscLogFlops(dof*(Q*P*P*P+Q*Q*P*P+Q*Q*Q*P)*NE*2);		\
    }									\
    /* PetscFunctionReturn(0); */					\
  } while (0)

static PetscErrorCode TensorContract_QPX_4_1_2_2(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  TensorContract_QPX(1,2,2,Rf,Sf,Tf,tmode,xx,yy);
  return 0;
}
static PetscErrorCode TensorContract_QPX_4_3_2_2(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  TensorContract_QPX(3,2,2,Rf,Sf,Tf,tmode,xx,yy);
  return 0;
}

static PetscErrorCode TensorContract_QPX_4_1_3_3(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  TensorContract_QPX(1,3,3,Rf,Sf,Tf,tmode,xx,yy);
  return 0;
}
static PetscErrorCode TensorContract_QPX_4_3_3_3(Tensor ten,const PetscReal Rf[],const PetscReal Sf[],const PetscReal Tf[],TensorMode tmode,const PetscScalar xx[],PetscScalar yy[]) {
  TensorContract_QPX(3,3,3,Rf,Sf,Tf,tmode,xx,yy);
  return 0;
}

#endif

// Choose our optimized functions if available
PetscErrorCode TensorSelect_QPX(Tensor ten) {

  PetscFunctionBegin;
#ifdef __bgq__
  if (ten->ne == 4) {
    PetscInt P = ten->P,Q = ten->Q;
    switch (ten->dof) {
    case 1: // Scalar problems with Q1 or Q2 elements
      if (P == 2 && Q == 2)      ten->Contract = TensorContract_QPX_4_1_2_2;
      else if (P == 3 && Q == 3) ten->Contract = TensorContract_QPX_4_1_3_3;
      break;
    case 3: // Coordinates or elasticity
      if (P == 2 && Q == 2)      ten->Contract = TensorContract_QPX_4_3_2_2;
      else if (P == 3 && Q == 3) ten->Contract = TensorContract_QPX_4_3_3_3;
      break;
    }
  }
#endif
  PetscFunctionReturn(0);
}
