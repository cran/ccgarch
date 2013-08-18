#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemm, a matrix multiplication */
#include <R_ext/Lapack.h> /* for dpotrf, the Cholesky decomposition */

SEXP eccc_sim(SEXP n, SEXP a0, SEXP Arch, SEXP Garch, SEXP V, SEXP inih, SEXP nu)
{
  int i, j, nobs = asInteger(n), ndim = LENGTH(a0), info, ione = 1;
  double one = 1.0, zero = 0.0, df = asReal(nu),
         *rz, *ra, *rA, *rB,  *rR, *rhini, *reps, *reps_row, *rh_row, *rh, *rhl, *rel2, *rz_, *rz_row;
  SEXP a, A, B, R, hini, eps, eps_row, z, h, hl, el2, z_, z_row, h_row, output;
  
  PROTECT(z = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(z_ = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(eps = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(h = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(eps_row = allocVector(REALSXP, ndim));
  PROTECT(h_row = allocVector(REALSXP, ndim));
  PROTECT(z_row = allocVector(REALSXP, ndim));
  PROTECT(el2 = allocVector(REALSXP, ndim));
  PROTECT(hl = duplicate(inih));
  PROTECT(output = allocVector(VECSXP, 2));
  PROTECT(a = duplicate(a0));
  PROTECT(A = duplicate(Arch));
  PROTECT(B = duplicate(Garch));
  PROTECT(R = duplicate(V));
  PROTECT(hini = duplicate(inih));
  ra = REAL(a);
  rA = REAL(A);
  rB = REAL(B);
  rR = REAL(R);
  rhini = REAL(hini);
  rz = REAL(z);
  rz_ = REAL(z_);
  rz_row = REAL(z_row);
  reps = REAL(eps);
  reps_row = REAL(eps_row);
  rh_row = REAL(h_row);
  rh = REAL(h);
  rhl = REAL(hl);
  rel2 = REAL(el2);

  for(j=0; j<ndim; j++){
    rel2[j] = 0.0;
    rhl[j] = 0.0;
    rh_row[j] = 0.0;
  }
  /* copying initial value of the conditional variance */
    F77_CALL(dcopy)(&ndim, rhini, &ione, rel2, &ione);
    F77_CALL(dcopy)(&ndim, rhini, &ione, rhl, &ione);

  /* carry out the Cholesky decomposition of R */
 	for (j = 0; j < ndim; j++) {   /* setting lower left block = 0 */
	  for (i = j+1; i < ndim; i++) {
		  rR[i + j * ndim] = 0.0;
	  }
	}
	F77_CALL(dpotrf)("Upper", &ndim, rR, &ndim, &info); /* rR is a upper triangular matrix such that t(rR)%*%rR = R */
                                                        /* "info" returns an index for the result of the Chlesky decomposition */
  /* generating a vector from N(0, I) or from t_df(0, I) */
    GetRNGstate();
      if(!R_FINITE(df)){
          for(i=0; i<(nobs*ndim); i++){
                rz[i] = rnorm(zero, one);
          }
      } else {
          for(i=0; i<(nobs*ndim); i++){
                   rz[i] = rt(df);    
          }
      }
    PutRNGstate();

  /* transforming z - N(0, I) into N(0, R) through rz%*%rV*/
    F77_CALL(dgemm)("N", "N", &nobs, &ndim, &ndim, &one, rz, &nobs, rR, &ndim, &zero, rz_, &nobs);  /* rz_ ~ N(0, R). */            

  /* the main loop: creating eps and conditional volatilities */
  for(i=0; i<nobs; i++){
    /* creating h */
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rA, &ndim, rel2, &ione, &zero, rh_row, &ione);   
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rB, &ndim, rhl, &ione, &one, rh_row, &ione);     
    F77_CALL(daxpy)(&ndim, &one, ra, &ione, rh_row, &ione);

    for(j=0; j<ndim; j++){
      rz_row[j] = rz_[i+j*nobs];
    }
    
    /* saving and initialising elements */
    for(j=0; j<ndim; j++){
      reps_row[j] = sqrt(rh_row[j])*rz_row[j];  /* reps_ = D%*%z */
      reps[i+j*nobs] = reps_row[j];             /* saving simulated eps */
      rh[i+j*nobs] = rh_row[j];                 /* saving simulated volatilities */
      rel2[j] = R_pow_di(reps_row[j],2);        /* eps_{-1}^2: used in the next step of loop */
      rhl[j] = rh_row[j];                       /* h_{-1}: used in the next step of loop */
    }
  }
  
  SET_VECTOR_ELT(output, 0, h);
  SET_VECTOR_ELT(output, 1, eps);

  UNPROTECT(15);
  return(output);
}
