#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemm, a matrix multiplication */
#include <R_ext/Lapack.h> /* for dpotrf, the Cholesky decomposition */

SEXP stcc_sim(SEXP n, SEXP a0, SEXP Arch, SEXP Garch, SEXP R, SEXP inih, SEXP nu){ 
  int i, j, k, nobs = asInteger(n), ndim = LENGTH(a0), cordim = ndim*ndim, info, ione = 1;
  double one = 1.0, zero = 0.0, df = asReal(nu), 
      *ra, *rA, *rB, *rhini, *rel2, *rh, *rhl, *rh_row,                           /* for conditional variance equation */
      *rvR, *tmprR,                                                               /* for correlation matrices */
      *reps, *rz, *rz_row, *rz_tmp, *reps_row;                                    /* for residuals or standardised residuals */
  
  SEXP a, A, B, hini, el2, h, hl, h_row, vR, tmpR,
       eps, z, z_row, z_tmp, eps_row,
       output;

  PROTECT(eps = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(z = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(z_row = allocVector(REALSXP, ndim));
  PROTECT(z_tmp = allocVector(REALSXP, ndim));
  PROTECT(eps_row = allocVector(REALSXP, ndim));
    
  /* a, A, B, hini, el2, h, hl, h_row */
  PROTECT(a = duplicate(a0));
  PROTECT(A = duplicate(Arch));
  PROTECT(B = duplicate(Garch));
  PROTECT(vR = duplicate(R));
  PROTECT(tmpR = allocMatrix(REALSXP, ndim, ndim));       /* a temporal matrix for R_{t}*/
  PROTECT(hini = duplicate(inih));
  PROTECT(el2 = allocVector(REALSXP, ndim));
  PROTECT(h = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(hl = allocVector(REALSXP, ndim));
  PROTECT(h_row = allocVector(REALSXP, ndim));

  PROTECT(output = allocVector(VECSXP, 2));
  
  reps = REAL(eps);
  rz = REAL(z);
  rz_row = REAL(z_row);
  rz_tmp = REAL(z_tmp);
  reps_row = REAL(eps_row);

  rvR = REAL(vR);
  tmprR = REAL(tmpR);

  ra = REAL(a);
  rA = REAL(A);
  rB = REAL(B);
  rhini = REAL(hini);
  rel2 = REAL(el2);
  rh = REAL(h);
  rhl = REAL(hl);
  rh_row = REAL(h_row);
/*****************************************************************************************/
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

  /* copying initial value of the conditional variance */
    F77_CALL(dcopy)(&ndim, rhini, &ione, rel2, &ione);
    F77_CALL(dcopy)(&ndim, rhini, &ione, rhl, &ione);
/*****************************************************************************************/
for(i=0; i<nobs; i++){
    for(j=0; j<cordim; j++){
      tmprR[j] = rvR[i+j*nobs]; 
    }
   	for (k = 0; k < ndim; k++) {   /* setting upper right block = 0 */
  	  for (j = k+1; j < ndim; j++) {
  		  tmprR[k + j * ndim] = 0.0;
  	  }
  	}
  /* carry out the Cholesky decomposition of R */
  	F77_CALL(dpotrf)("Lower", &ndim, tmprR, &ndim, &info); /* tmprR is a lower triangular matrix */
    for(j=0; j<ndim; j++){
      rz_tmp[j] = rz[i+j*nobs];     /* i is determined in the other loop */
    }
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, tmprR, &ndim, rz_tmp, &ione, &zero, rz_row, &ione);   /* rz_row is the standardised residual */
  /* creating h */
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rA, &ndim, rel2, &ione, &zero, rh_row, &ione);   
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rB, &ndim, rhl, &ione, &one, rh_row, &ione);     
    F77_CALL(daxpy)(&ndim, &one, ra, &ione, rh_row, &ione);
    /* saving elements */
    for(j=0; j<ndim; j++){
      reps_row[j] = sqrt(rh_row[j])*rz_row[j];
      reps[i+j*nobs] = reps_row[j];           /* saving simulated eps */
      rh[i+j*nobs] = rh_row[j];               /* saving simulated volatilities */
      rel2[j] = R_pow_di(reps_row[j],2);      /* eps_{-1}^2: used in the next step of loop */
      rhl[j] = rh_row[j];                     /* h_{-1}: used in the next step of loop */
    }
}
/*****************************************************************************************/
  SET_VECTOR_ELT(output, 0, h);
  SET_VECTOR_ELT(output, 1, eps);
      
  UNPROTECT(16);
  return(output);
}
