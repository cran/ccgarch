#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Applic.h> /* for dgemv and daxpy */

SEXP vector_garch(SEXP res, SEXP a0, SEXP Arch, SEXP Garch)
{
  int i, j, nobs, ndim, ione = 1;
  double *rx, *ra, *rA, *rB, *rh, *rtmp, *rx1, *rx2,
          one = 1.0, zero = 0.0;  /* x = eps^2. */
  SEXP h, tmp, x1, x2, a, A, B, x;
  
  nobs = Rf_nrows(res); ndim = Rf_ncols(res);
  PROTECT(tmp = allocVector(REALSXP, ndim));
  PROTECT(x1 = allocVector(REALSXP, ndim));
  PROTECT(x2 = allocVector(REALSXP, ndim));
  PROTECT(a = duplicate(a0));
  PROTECT(A = duplicate(Arch));
  PROTECT(B = duplicate(Garch));
  PROTECT(x = duplicate(res));
  rtmp = REAL(tmp);
  ra = REAL(a);
  rA = REAL(A);
  rB = REAL(B);
  rx = REAL(x);
  rx1 = REAL(x1);   /* eps_{0}^2*/
  rx2 = REAL(x2);   /* h_{0} */
    
  PROTECT(h = allocMatrix(REALSXP, nobs, ndim));  
  rh = REAL(h);                                   

  for(j=0; j<ndim; j++){
    rx1[j] = 0.0;
    rx2[j] = 0.0;
  }

  /* initial values = the average of eps^2 */
  for(j=0; j<ndim; j++){
    for(i=0; i<nobs; i++){
      rx1[j] += rx[i + j*nobs]/nobs ;   /* eps_{0}^2*/
      rx2[j] += rx[i + j*nobs]/nobs ;   /* h_{0} */
    }
  }
  
  /* recursion from T=1 to the end */
  for(i=0; i<nobs; i++){
     F77_CALL(dgemv)("N", &ndim, &ndim, &one, rA, &ndim, rx1, &ione, &zero, rtmp, &ione); 
     F77_CALL(dgemv)("N", &ndim, &ndim, &one, rB, &ndim, rx2, &ione, &one, rtmp, &ione); 
     F77_CALL(daxpy)(&ndim, &one, ra, &ione, rtmp ,&ione); 
    
    for(j=0; j<ndim; j++){
      rh[i+j*nobs] = rtmp[j];   /* setting elements in h */
      rx1[j] = rx[i+j*nobs];    /* for the next loop */
      rx2[j] = rtmp[j];         /* for the next loop */
    }
  }
  
  UNPROTECT(8);
  return(h);
}
