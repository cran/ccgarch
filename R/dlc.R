# the gradient and Hessian of the DCC part of the log-likelihood function
dlc <- function(dcc.para, B, u, h, model){
   nobs <- dim(h)[1]
   ndim <- dim(h)[2]
   dccpar1 <- dcc.para[1]
   dccpar2 <- dcc.para[2]
   dhw <- vec.garch.derivative(u, B, h)
      inv.sq.h <- 1/sqrt(h)
      z <- u*inv.sq.h
   Q <- cov(z)
   dcc <- dcc.est(z, dcc.para)
   P <- dcc$DCC
   Qt <- dcc$Q
   
   In <- diag(ndim)
      ind <- as.vector(rbind(1, In, In))
      if(model=="diagonal"){
         npar.h <- 3*ndim
      } else {
         npar.h <- ndim*(2*ndim+1)
      }
#   Qtilde <- 1/sqrt(Qt[,as.vector(In)==1])
   Qtilde <- sqrt(Qt[,as.vector(In)==1])      # modified on 09.10.11
   dvecQ <- matrix(0, nobs, 2*ndim^2)
   dvecP <- matrix(0, nobs, 2*ndim^2)
   dlc <- matrix(0, nobs, 2)
   d2lc <- matrix(0, nobs, 4)
   dvecQt <- matrix(0, 2, ndim^2)
      dwdvecQt <- matrix(0, npar.h, ndim^2)
      dwdvecQ <- matrix(0, nobs, npar.h*ndim^2)
      dfdwd2lc <- matrix(0, nobs, 2*npar.h)
   const <- -rbind(as.vector(Q), as.vector(Q))
   for ( i in 2:nobs){
      zz <- outer(z[i,], z[i,])
      dvecQt <- const + rbind(as.vector(zz), Qt[i-1,]) + dccpar2*dvecQt
      dvecQ[i,] <- as.vector(dvecQt)
      invQ.tilde <- diag(1/Qtilde[i,])
      Pt <- matrix(P[i,], ndim, ndim)
      invPt <- solve(Pt)
      Z <- Pt%x%invQ.tilde
      QPQ <- invQ.tilde%x%invQ.tilde - 0.5*diag(as.vector(invQ.tilde))%*%( Z + t(Z) )
      dvecPt <- dvecQt%*%QPQ
      dvecPt[,as.vector(In)==1] <- 0
      dvecP[i,] <- as.vector(dvecPt)
      dlc[i,] <- -0.5*as.vector(dvecPt%*%as.vector(invPt - invPt%*%zz%*%invPt))
   
      d2lc[i,] <- as.vector(dvecPt%*%(invPt%x%invPt)%*%t(dvecPt))
         dhwt <- matrix(dhw[i,],ncol=ndim)
         if(model=="diagonal"){
            dhwt <- dhwt[ind==1,]
         }
         dwdvecVt <- matrix(0, npar.h, ndim^2)
         for(j in 1:ndim){
            dwdvecVt[,(j+(j-1)*ndim)] <- dhwt[,j]
         }
         invD <- diag(inv.sq.h[i,])
         dgD <- diag(as.vector(invD))
         ZD <- outer(z[i,], z[i,])%x%invD
         
         dwdvecQt <- -0.5*dccpar1*dwdvecVt%*%dgD%*%(ZD + t(ZD)) + dccpar2*dwdvecQt
         dwdvecQ[i,] <- as.vector(dwdvecQt)
         dwdvecPt <- t(dwdvecQt%*%QPQ)
         
         dfdwd2lc[i,] <- as.vector(0.5*dvecPt%*%(invPt%x%invPt)%*%dwdvecPt + 0.25*dvecPt%*%((invPt%*%invD)%x%In + In%x%(invPt%*%invD))%*%t(dwdvecVt))
   }
   list(dlc=dlc, dvecP = dvecP, dvecQ = dvecQ, d2lc=0.5*d2lc, dfdwd2lc=dfdwd2lc)
}
