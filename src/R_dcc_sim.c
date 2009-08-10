#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemm, a matrix multiplication */
#include <R_ext/Lapack.h> /* for dpotrf, the Cholesky decomposition */
SEXP dcc_sim(SEXP n, SEXP a0, SEXP Arch, SEXP Garch, SEXP inih, SEXP uncR, SEXP dcca, SEXP dccb, SEXP nu){       /* nobs, uncR, alpha, beta, df */
  int i, j, k, nobs = asInteger(n), ndim = LENGTH(a0), cordim = ndim*ndim, info, ione = 1;
  double one = 1.0, zero = 0.0, df = asReal(nu),
      *rQ, *rQbar, *rQlag, *rR, *rdiagQ, *tmprQ, *tmprR, *rDCC, *rdab, *rd_a, *rd_b,  /* for DCC equation */
      *ra, *rA, *rB, *rhini, *rel2, *rh, *rhl, *rh_row,                      /* for conditional variance equation */
      *rz, *rz_row, *rz_tmp, *rstd_z, *reps, *reps_row;                                    /* for residuals or standardised residuals */
  
  SEXP Q, Qbar, Qlag, diagQ, tmpQ, tmpR, DCC, R, dab, d_a, d_b, 
       a, A, B, hini, el2, h, hl, h_row,
       z, z_row, z_tmp, std_z, eps, eps_row, 
       output;
  
  PROTECT(Q = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(Qbar = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(Qlag = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(diagQ = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(tmpQ = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(tmpR = allocMatrix(REALSXP, ndim, ndim));
  PROTECT(DCC = allocMatrix(REALSXP, nobs, cordim));
  PROTECT(R = duplicate(uncR));                       /* unconditional correlation */
  PROTECT(dab = allocVector(REALSXP, 1));
  PROTECT(d_a = duplicate(dcca));                       /*  */
  PROTECT(d_b = duplicate(dccb));                       /*  */
  PROTECT(a = duplicate(a0));
  PROTECT(A = duplicate(Arch));
  PROTECT(B = duplicate(Garch));
  PROTECT(hini = duplicate(inih));
  PROTECT(el2 = allocVector(REALSXP, ndim));
  PROTECT(h = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(hl = allocVector(REALSXP, ndim));
  PROTECT(h_row = allocVector(REALSXP, ndim));
  PROTECT(z = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(z_row = allocVector(REALSXP, ndim));
  PROTECT(z_tmp = allocVector(REALSXP, ndim));
  PROTECT(std_z = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(eps = allocMatrix(REALSXP, nobs, ndim));
  PROTECT(eps_row = allocVector(REALSXP, ndim));

  PROTECT(output = allocVector(VECSXP, 5));
  
  /* DCC part */
  rQ = REAL(Q);
  rQbar = REAL(Qbar);
  rQlag = REAL(Qlag);
  rdiagQ = REAL(diagQ);
  tmprQ = REAL(tmpQ);
  tmprR = REAL(tmpR);
  rDCC = REAL(DCC);
  rR = REAL(R);
  rdab = REAL(dab);
  rd_a = REAL(d_a);
  rd_b = REAL(d_b);
  
  /* standardised residuals and residuals */
  rz = REAL(z);
  rz_row = REAL(z_row);
  rz_tmp = REAL(z_tmp);
  rstd_z = REAL(std_z);
  reps = REAL(eps);
  reps_row = REAL(eps_row);

  /* volatility part */
  ra = REAL(a);
  rA = REAL(A);
  rB = REAL(B);
  rhini = REAL(hini);
  rel2 = REAL(el2);
  rh = REAL(h);
  rhl = REAL(hl);
  rh_row = REAL(h_row);

/*****************************************************************************************/
  /* a coefficient on rQbar */
    rdab[0] = 1.0-rd_a[0]-rd_b[0];
  /* generating a vector from N(0, I) or from t_df(0, I) */
    GetRNGstate();
      if(!R_FINITE(df)){
          for(i=0; i<(nobs*ndim); i++){
                rz[i] = rnorm(zero, one);
          }
          for(j=0; j<ndim; j++){
              rz_tmp[j] = rnorm(zero, one);      /* initial values */
          }
      } else {
          for(i=0; i<(nobs*ndim); i++){
                   rz[i] = rt(df);    
          }
          for(j=0; j<ndim; j++){
              rz_tmp[j] = rt(df);
          }
      }
    PutRNGstate();

  /* creating Qbar matrix */
    for(j=0; j<cordim; j++){
      rdiagQ[j] = 0.0;
      rQbar[j] = rdab[0]*rR[j];
      rQlag[j] = rQbar[j];
    }
  /* copying initial value of the conditional variance */
    F77_CALL(dcopy)(&ndim, rhini, &ione, rel2, &ione);
    F77_CALL(dcopy)(&ndim, rhini, &ione, rhl, &ione);
    /* F77_CALL(dcopy)(&ndim, rQbar, &ione, rQlag, &ione); */

  /* carry out the Cholesky decomposition of R */
   	for (k = 0; k < ndim; k++) {   
  	  for (j = k+1; j < ndim; j++) {
  		  rR[k + j * ndim] = 0.0;        /* setting upper right block = 0 */
  	  }
  	}
  	F77_CALL(dpotrf)("Lower", &ndim, rR, &ndim, &info); /* rR is a lower triangular matrix */
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rR, &ndim, rz_tmp, &ione, &zero, rz_row, &ione);   /* rz_row is the standardised residual */
    
/*****************************************************************************************/
for(i=0; i<nobs; i++){
  /* outer product of z_(t-1) */
    F77_CALL(dgemm)("N", "T", &ndim, &ndim, &ione, rd_a, rz_row, &ndim, rz_row, &ndim, rd_b, rQlag, &ndim);      /* a*ee'+b*Q */
  /* creating Q matrix for obserbation t */
    for(j=0; j<cordim; j++){
      rQ[j] = rQbar[j]+rQlag[j];        /* Qbar+(a*ee'+b*Q) */
      rQlag[j] = rQ[j]; /* saving for the next round */
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
    }
  /* carry out the Cholesky decomposition of R */
   	for (k = 0; k < ndim; k++) {   
  	  for (j = k+1; j < ndim; j++) {
  		  tmprR[k + j * ndim] = 0.0;        /* setting upper right block = 0 */
  	  }
  	}
  	F77_CALL(dpotrf)("Lower", &ndim, tmprR, &ndim, &info); /* tmprR is a lower triangular matrix */
    for(j=0; j<ndim; j++){
      rz_tmp[j] = rz[i+j*nobs];     /* i is determined in the other loop */
    }
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, tmprR, &ndim, rz_tmp, &ione, &zero, rz_row, &ione);   /* rz_row is the standardised residual */

    /* saving and initialising elements */
    for(j=0; j<ndim; j++){
      rstd_z[i+j*nobs] = rz_row[j];
    }
  /* creating h */
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rA, &ndim, rel2, &ione, &zero, rh_row, &ione);   
    F77_CALL(dgemv)("N", &ndim, &ndim, &one, rB, &ndim, rhl, &ione, &one, rh_row, &ione);     
    F77_CALL(daxpy)(&ndim, &one, ra, &ione, rh_row, &ione);
    for(j=0; j<ndim; j++){
      reps_row[j] = sqrt(rh_row[j])*rz_row[j];    /* reps_row = D%*%z */
    /*  */
      reps[i+j*nobs] = reps_row[j];           /* saving simulated eps */
      rh[i+j*nobs] = rh_row[j];               /* saving simulated volatilities */
      rhl[j] = rh_row[j];                     /* h_{t-1}: for the next loop */
      /* rel2[j] = R_pow_di(reps[i+j*nobs], 2);      eps_{-1}^2: used in the next step of loop */
      rel2[j] = R_pow_di(reps_row[j], 2);  /*     eps_{-1}^2: used in the next step of loop */
    }
}
/*****************************************************************************************/
  SET_VECTOR_ELT(output, 0, z);
  SET_VECTOR_ELT(output, 1, std_z);
  SET_VECTOR_ELT(output, 2, DCC);
  SET_VECTOR_ELT(output, 3, h);
  SET_VECTOR_ELT(output, 4, eps);

  UNPROTECT(26);
  return(output);
}
