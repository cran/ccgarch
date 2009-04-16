      dcc.estimation <- function(inia, iniA, iniB, ini.dcc, dvar, model, method="BFGS", gradient=1, message=1){
         dvar <- as.matrix(dvar)
         ndim <- dim(dvar)[2]
         In <- diag(ndim)
         
         if(!is.matrix(iniA)||!is.matrix(iniB)){
            stop("iniA or iniB or both must be matrices")
         }
         
         first.stage <- dcc.estimation1(dvar=dvar, a=inia, A=iniA, B=iniB, model=model, method=method)
       #  if(first.stage$convergence != 0){
       #     cat("********************************************\n")
       #     cat("* The first stage optimization is failed.  *\n")
       #     cat("* See a list variable "first" for details. *\n")
       #     cat("********************************************\n")
       #  }
            # cat("**************************************************\n")
         tmp.para <- c(first.stage$par, In[lower.tri(In)])
         estimates <- p.mat(tmp.para, model=model, ndim=ndim)
         esta <- estimates$a
         estA <- estimates$A
         estB <- estimates$B

         h <- vector.garch(dvar, esta, estA, estB)    # estimated conditional variances
         std.resid <- dvar/sqrt(h)                    # std. residuals

         second.stage <- dcc.estimation2(std.resid, ini.dcc, gradient=gradient)
         
         if(second.stage$convergence != 0){
            cat("*********************************************\n")
            cat("* The second stage optimization is failed.  *\n")
            cat("* See a list variable 'second' for details. *\n")
            cat("*********************************************\n")
         } else if(message != 0) {
            cat("***********************************************************\n")
            cat("*  Estimation has been completed.                         *\n")
            cat("*  The outputs are saved in a list with components:       *\n")
            cat("*    out    : the estimates and their standard errors     *\n")
            cat("*    h      : a matrix of estimated conditional variances *\n")
            cat("*    DCC    : a matrix of DCC estimates                   *\n")
            cat("*    std.resid : a matrix of the standardised residuals   *\n")
            cat("*    first  : the results of the first stage estimation   *\n")
            cat("*    second : the results of the second stage estimation  *\n")
            cat("***********************************************************\n")
         }

         dcc <- dcc.est(std.resid, second.stage$par)$DCC
         
         output <- dcc.results(dvar, first.stage$par, second.stage$par, h, model=model)
         list(out=output, h=h, DCC=dcc, std.resid=std.resid, first=first.stage, second=second.stage)
      }
