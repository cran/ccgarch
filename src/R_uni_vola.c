#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

SEXP uni_vola(SEXP a, SEXP usq)
{
  int n = length(usq), i;
  double *rusq = REAL(usq), *para = REAL(a), *rh, m1 = 0.0, m2 = 0.0;
  SEXP ansh;
  
  PROTECT(ansh = allocVector(REALSXP, n));
  rh = REAL(ansh);
  
  for(i=0; i<n;i++){
    m1 += rusq[i]/n;  /* setting e_{0}^2 = average(usq)*/
    m2 += rusq[i]/n;  /* setting h_{0}^2 = average(usq)*/
  }
  
  for(i=0; i < n; i++){
    rh[i] = para[0] + para[1]*m1+ para[2]*m2;
    m1 = rusq[i];
    m2 = rh[i];
  }
  UNPROTECT(1);
  return(ansh);
}

