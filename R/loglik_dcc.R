# the full log-likelihood function of the DCC-GARCH(1,1) model
loglik.dcc <- function (param, dvar, model) {
    nobs <- dim(dvar)[1]
    ndim <- dim(dvar)[2]
    Id <- diag(ndim)
    npar <- length(param)
    dcc.param <- param[(npar-1):npar]
    param <- c(param[1:(npar-2)], rep(0, ndim*(ndim - 1)/2))
    para.mat <- p.mat(param, model, ndim)
    a <- para.mat$a
    A <- para.mat$A
    B <- para.mat$B
    h <- vector.garch(dvar, a, A, B)
    z <- dvar/sqrt(h)
    DCC <- dcc.est(z, dcc.param)$DCC
    
    lf1 <- -0.5*ndim*log(2*pi) - rowSums(log(h)) -0.5*rowSums(z^2)
    lf2 <- numeric(nobs)
   for( i in 1:nobs){                        
      R <- matrix(DCC[i,], ndim, ndim)
      invR <- solve(R)
      lf2[i] <- -0.5*(log(det(R)) +sum(z[i,]*crossprod(invR,z[i,])))
   }
    lf1 + lf2
}
