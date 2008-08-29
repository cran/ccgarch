/**************************************************
Simulating univariate GARCH(1,1) volatility:
dyn.load("R_uni_vola_t_sim.dll")
uni.vola.sim = function(a, nobs, df = Inf, cut = 1000){
   nobs = nobs + cut
   cond.h = .Call("uni_vola_sim", nobs, a, df)
   h = cond.h[[1]]; eps = cond.h[[2]]
   list(h = h[(cut+1):nobs], eps = eps[(cut+1):nobs])
}
nobs = 1000; cut = 1000; df = 8
a = c(0.1,0.2,0.7)
temp.t = uni.vola.sim(a, nobs, df = df, cut= cut)
temp = uni.vola.sim(a, nobs, cut = cut)
par(mfcol=c(2,2))
	plot(temp.t$h,type="l")
	plot(temp.t$eps,type="l")
	plot(temp$h,type="l")
	plot(temp$eps,type="l")


dyn.unload("R_uni_vola_t_sim.dll")
**************************************************/
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

SEXP uni_vola_sim(SEXP n, SEXP a, SEXP df)
{
  int nobs = asInteger(n), i;
  double *rh, *reps, *rz, *rpara = REAL(a), *rel2, *rhl, nu = asReal(df);
  SEXP h, eps, z, el2, hl, output;
  
  PROTECT(h = allocVector(REALSXP, nobs));
  PROTECT(eps = allocVector(REALSXP, nobs));
  PROTECT(z = allocVector(REALSXP, nobs));
  PROTECT(output = allocVector(VECSXP, 2));
  el2 = allocVector(REALSXP, 1);
  hl = allocVector(REALSXP,1);
  rh = REAL(h);
  reps = REAL(eps);
  rz = REAL(z);
  rel2 = REAL(el2);
  rhl = REAL(hl);
  
  /* generating a vector from t(0, 1). If df = Inf, from N(0, 1). */
    GetRNGstate();
      for(i=0; i<nobs; i++){
          rz[i] = rt(nu);
      }
    PutRNGstate();
  /* initial observation */
  rel2[0] = rpara[0]/(1.0 - rpara[1] - rpara[2]);
  rhl[0] = rpara[0]/(1.0 - rpara[1] - rpara[2]);
  
  for(i=0; i < nobs; i++){
    rh[i] = rpara[0] + rpara[1]*rel2[0]+ rpara[2]*rhl[0];
    reps[i] = sqrt(rh[i])*rz[i];
    rel2[0] = R_pow_di(reps[i],2);    /* used in the next loop */
    rhl[0] = rh[i];                   /* used in the next loop */
  }

  SET_VECTOR_ELT(output, 0, h);
  SET_VECTOR_ELT(output, 1, eps);

  UNPROTECT(4);
  return(output);
}

