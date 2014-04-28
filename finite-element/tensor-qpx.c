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
      vector4double R[Q][P],S[Q][P],T[Q][P];				\
      PetscScalar (*restrict x)[P*P*P][NE] = (PetscScalar(*)[P*P*P][NE])xx; /* should be const, but xlc vector intrinsics don't understand const */ \
      PetscScalar (*restrict y)[P*P*P][NE] = (PetscScalar(*)[Q*Q*Q][NE])yy; \
      PetscScalar u[dof][Q*P*P][NE]_align,v[dof][Q*Q*P][NE]_align;	\
									\
      for (PetscInt i=0; i<Q; i++) {					\
	for (PetscInt j=0; j<P; j++) {					\
	  R[i][j] = vec_splats(tmode == TENSOR_EVAL ? Rf[i*P+j] : Rf[j*Q+i]); \
	  S[i][j] = vec_splats(tmode == TENSOR_EVAL ? Sf[i*P+j] : Sf[j*Q+i]); \
	  T[i][j] = vec_splats(tmode == TENSOR_EVAL ? Tf[i*P+j] : Tf[j*Q+i]); \
	}								\
      }									\
									\
      /* u[l,a,j,k] = R[a,i] x[l,i,j,k] */				\
      for (PetscInt l=0; l<dof; l++) {					\
	for (PetscInt jk=0; jk<P*P; jk++) {				\
	  vector4double u_lAjk[Q],x_lIjk[P];				\
	  Pragma(loopid(ld_x_lIjk)) for (PetscInt i=0; i<P; i++) x_lIjk[i] = vec_lda(0,x[l][i*P*P+jk]); \
	  Pragma(loopid(z_u_lAjk)) for (PetscInt a=0; a<Q; a++) u_lAjk[a] = vec_splats(0.); \
	  Pragma(unroll(P_)) Pragma(loopid(R_aI)) for (PetscInt i=0; i<P; i++) { \
	    Pragma(unroll(Q_)) Pragma(loopid(R_Ai)) for (PetscInt a=0; a<Q; a++) { \
	      u_lAjk[a] = vec_madd(R[a][i],x_lIjk[i],u_lAjk[a]);	\
	    }								\
	  }								\
	  Pragma(loopid(st_u_lAjk)) for (PetscInt a=0; a<Q; a++) vec_sta(u_lAjk[a],0,u[l][a*P*P+jk]); \
	}								\
      }									\
									\
      /* v[l,a,b,k] = S[b,j] u[l,a,j,k] */				\
      for (PetscInt l=0; l<dof; l++) {					\
	for (PetscInt a=0; a<Q; a++) {					\
	  for (PetscInt k=0; k<P; k++) {				\
	    vector4double v_laBk[Q],u_laJk[P];				\
	    Pragma(loopid(ld_u_laJk)) for (PetscInt j=0; j<P; j++) u_laJk[j] = vec_lda(0,u[l][(a*P+j)*P+k]); \
	    Pragma(loopid(z_v_laBk)) for (PetscInt b=0; b<Q; b++) v_laBk[b] = vec_splats(0.); \
	    Pragma(unroll(P_)) Pragma(loopid(S_bJ)) for (PetscInt j=0; j<P; j++) { \
	      Pragma(unroll(Q_)) Pragma(loopid(S_Bj)) for (PetscInt b=0; b<Q; b++) { \
		v_laBk[b] = vec_madd(S[b][j],u_laJk[j],v_laBk[b]);	\
	      }								\
	    }								\
	    Pragma(loopid(st_v_laBk)) for (PetscInt b=0; b<Q; b++) vec_sta(v_laBk[b],0,v[l][(a*Q+b)*P+k]); \
	  }								\
	}								\
      }									\
									\
      /* y[l,a,b,c] = T[c,k] v[l,a,b,k] */				\
      for (PetscInt l=0; l<dof; l++) {					\
        for (PetscInt ab=0; ab<Q*Q; ab++) {				\
	  vector4double y_labC[Q],v_labK[P];				\
	  Pragma(loopid(ld_v_labK)) for (PetscInt k=0; k<P; k++) v_labK[k] = vec_lda(0,v[l][ab*P+k]); \
	  Pragma(loopid(ld_y_labC)) for (PetscInt c=0; c<Q; c++) y_labC[c] = vec_lda(0,y[l][ab*Q+c]); \
	  Pragma(unroll(P_)) Pragma(loopid(T_cK)) for (PetscInt k=0; k<P; k++) { \
	    Pragma(unroll(Q_)) Pragma(loopid(T_Ck)) for (PetscInt c=0; c<Q; c++) { \
	      y_labC[c] = vec_madd(T[c][k],v_labK[k],y_labC[c]);	\
	    }								\
	  }								\
	  Pragma(loopid(st_y_labC)) for (PetscInt c=0; c<Q; c++) vec_sta(y_labC[c],0,y[l][ab*Q+c]); \
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
