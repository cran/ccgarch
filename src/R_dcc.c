#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemm, a matrix multiplication */
/* #include <R_ext/Lapack.h> for dpotrf, the Cholesky decomposition */

SEXP dcc_est(SEXP z, SEXP R, SEXP dcca, SEXP dccb){
  int i, j, nobs = Rf_nrows(z), ndim= Rf_ncols(z), cordim = ndim*ndim, ione = 1;
  double one = 1.0, zero = 0.0, 
      *rQ, *rQbar, *rQlag, *rdiagQ, *tmprQ, *tmprR, *rDCC, *routQ, *rdab, *rd_a, *rd_b,  /* for DCC equation */
      *rz_, *rz_lag;                                    /* for residuals or standardised residuals */
  
  SEXP Q, Qbar, Qlag, diagQ, tmpQ, tmpR, DCC, outQ, dab, d_a, d_b, 
       z_, z_lag,
       output;

  PROTECT(Q = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(Qbar = allocMatrix(REALSXP, ndim, ndim));                      /* the unconditional correlation matrix */
  PROTECT(Qlag = duplicate(R));
  PROTECT(diagQ = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(tmpQ = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(tmpR = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(DCC = allocMatrix(REALSXP, nobs, cordim));
  PROTECT(outQ = allocMatrix(REALSXP, nobs, cordim));
  PROTECT(dab = allocVector(REALSXP, 1));
  PROTECT(d_a = duplicate(dcca));                       /*  */
  PROTECT(d_b = duplicate(dccb));                       /*  */
  
  PROTECT(z_ = duplicate(z));
  PROTECT(z_lag = allocVector(REALSXP, ndim));

  PROTECT(output = allocVector(VECSXP, 2));

  rQ = REAL(Q);
  rQbar = REAL(Qbar);
  rQlag = REAL(Qlag);
  rdiagQ = REAL(diagQ);
  tmprQ = REAL(tmpQ);
  tmprR = REAL(tmpR);
  rDCC = REAL(DCC);
  routQ = REAL(outQ);
  rdab = REAL(dab);
  rd_a = REAL(d_a);
  rd_b = REAL(d_b);
  
  rz_ = REAL(z_);
  rz_lag = REAL(z_lag);

  /* a coefficient on rQbar */
    rdab[0] = 1.0-rd_a[0]-rd_b[0];
  /* initial values = the average of eps^2 */
  for(j=0; j<ndim; j++){
      rz_lag[j] = 0.0;
    for(i=0; i<nobs; i++){
      rz_lag[j] += rz_[i + j*nobs]/nobs ;   
    }
  }

  /* assigning parameters to Qbar matrix and initial values */
    for(j=0; j<cordim; j++){
      rdiagQ[j] = 0.0;
      rQbar[j] = rdab[0]*rQlag[j]; 
    }

for(i=0; i<nobs; i++){
  /* outer product of z_(t-1) */
    F77_CALL(dgemm)("N", "T", &ndim, &ndim, &ione, rd_a, rz_lag, &ndim, rz_lag, &ndim, rd_b, rQlag, &ndim);      /* a*ee'+b*Q */

  /* creating Q matrix for observation t */
    for(j=0; j<cordim; j++){
        rQ[j] = rQbar[j]+rQlag[j];        /* Qbar+(a*ee'+b*Q) */
        rQlag[j] = rQ[j];                            /* for the next round */
    }
  /* creating a dcc matrix for observation t */
    for(j=0; j<ndim; j++){
        rdiagQ[(ndim+1)*j] = 1.0/sqrt(rQ[(ndim+1)*j]);
    }
    F77_CALL(dgemm)("N", "N", &ndim, &ndim, &ndim, &one, rdiagQ, &ndim, rQ, &ndim, &zero, tmprQ, &ndim);      /* */
    F77_CALL(dgemm)("N", "N", &ndim, &ndim, &ndim, &one, tmprQ, &ndim, rdiagQ, &ndim, &zero, tmprR, &ndim);   /* */
  /* tmprR must be saved */
    for(j=0; j<cordim; j++){
      rDCC[i+j*nobs] = tmprR[j];
      routQ[i+j*nobs] = rQ[j];
    }
    for(j=0; j<ndim; j++){
       rz_lag[j] = rz_[i+j*nobs];                  /* for the next round */
    }
}
  SET_VECTOR_ELT(output, 0, DCC);
  SET_VECTOR_ELT(output, 1, outQ);

  UNPROTECT(14);
  return(output); 
}

