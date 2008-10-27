#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemm, a matrix multiplication */
#include <R_ext/Lapack.h> /* for dpotrf, the Cholesky decomposition */

SEXP eccc_sim(SEXP n, SEXP a0, SEXP Arch, SEXP Garch, SEXP V, SEXP inih, SEXP nu)
{
  int i, j, nobs = asInteger(n), ndim = Rf_nrows(V), info, ione = 1;
  double one = 1.0, zero = 0.0, df = asReal(nu),
         *rz, *ra, *rA, *rB,  *rR, *rhini, *rD, *reps, *reps_row, *rh_row, *rh, *rhl, *rel2, *rz_, *rz_row, *ry;
  SEXP a, A, B, R, hini, D, eps, eps_row, z, h, hl, el2, z_, z_row, h_row, y, output;
  
  PROTECT(z = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(z_ = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(eps = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(h = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(D = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(eps_row = allocVector(REALSXP, ndim));
  PROTECT(h_row = allocVector(REALSXP, ndim));
  PROTECT(z_row = allocVector(REALSXP, ndim));
  PROTECT(el2 = allocVector(REALSXP, ndim));
  PROTECT(hl = duplicate(inih));
  PROTECT(output = allocVector(VECSXP, 2));
  PROTECT(y = allocVector(REALSXP,1 ));
  a = duplicate(a0);
  A = duplicate(Arch);
  B = duplicate(Garch);
  R = duplicate(V);
  hini = duplicate(inih);
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
  rD = REAL(D);
  rhl = REAL(hl);
  rel2 = REAL(el2);
  ry = REAL(y);

  for(j=0; j<ndim; j++){
    rel2[j] = rhini[j];
    rhl[j] = rhini[j];
    rh_row[j] = 0.0;
    for(i=0; i<ndim; i++){
      rD[i+j*ndim] = 0.0;
    }
  }

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
                rz[i] = norm_rand();
          }
      } else {
          for(i=0; i<nobs; i++){
              ry[0] = sqrt(df/rchisq(df));
                 for(j=0; j<ndim; j++){
                   rz[i + j*nobs] = norm_rand()*ry[0];    
                 }
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
      if(ISNAN(rh_row[j]))
        error("'rtmp' contains 'NaN'");
      rD[j*(ndim+1)] = sqrt(rh_row[j]);
      rz_row[j] = rz_[i+j*nobs];
    }
    /* creating eps */
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rD, &ndim, rz_row, &ione, &zero, reps_row, &ione);   /* reps_ = D%*%z */
    
    /* saving and initialising elements */
    for(j=0; j<ndim; j++){
      reps[i+j*nobs] = reps_row[j];           /* saving simulated eps */
      rh[i+j*nobs] = rh_row[j];               /* saving simulated volatilities */
      rel2[j] = R_pow_di(reps_row[j],2);      /* eps_{-1}^2: used in the next step of loop */
      rhl[j] = rh_row[j];                     /* h_{-1}: used in the next step of loop */
      reps_row[j] = 0.0;
      rh_row[j] = 0.0;
    }
  }
  
  SET_VECTOR_ELT(output, 0, h);
  SET_VECTOR_ELT(output, 1, eps);
      
  UNPROTECT(12);
  return(output);
}

