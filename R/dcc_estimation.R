      dcc.estimation <- function(inia, iniA, iniB, ini.dcc, dvar, model){
         ndim <- dim(dvar)[2]
         In <- diag(ndim)
            cat("***********************************************************\n")
            cat("*  Doing the first stage estimation..............         *\n")
            cat("*                                                         *\n")
         first.stage <- dcc.estimation1(dvar=dvar, a=inia, A=iniA, B=iniB, model=model)
         if(first.stage$convergence == 0){
            cat("*  The first stage completed!                             *\n")
            cat("*                                                         *\n")
            cat("*  Proceed to the second stage estimation.                *\n")
            cat("*                                                         *\n")
            cat("*  Please wait for a while.......................         *\n")
         } else {
            cat("* Optimization is FAILED.                                 *\n")
            cat("* Fine tuning is required. Unfortunately, this has to be  *\n")
            cat("* done manually. However, some functions like             *\n")
            cat("* loglik.dcc1 can be used. See the manual for details.    *\n")
            cat("***********************************************************\n")
         }
            # cat("**************************************************\n")
         tmp.para <- c(first.stage$par, In[lower.tri(In)])
         estimates <- p.mat(tmp.para, model=model, ndim=ndim)
         esta <- estimates$a
         estA <- estimates$A
         estB <- estimates$B

         h <- vector.garch(dvar, esta, estA, estB)    # estimated volatility
         std.resid <- dvar/sqrt(h)                          # std. residuals

         second.stage <- dcc.estimation2(std.resid, ini.dcc, gradient=1)
         
         if(second.stage$convergence != 0){
            cat("* Optimization is failed.                                  \n")
            cat("***********************************************************\n")
         } else {
            cat("*                                                         *\n")
            cat("*  All the estimation has been completed.                 *\n")
            cat("*  The outputs are saved in a list with components:       *\n")
            cat("*    out    : the estimates and their standard errors     *\n")
            cat("*    h      : a matrix of volatility estimates            *\n")
            cat("*    DCC    : a matrix of DCC estimates                   *\n")
            cat("*    first  : the results of the first stage estimation   *\n")
            cat("*    second : the results of the second stage estimation  *\n")
            cat("***********************************************************\n")
         }

         dcc <- dcc.est(std.resid, second.stage$par)$DCC
         
         output <- dcc.results(dvar, first.stage$par, second.stage$par, h, model=model)
         list(out=output, h=h, DCC=dcc, first=first.stage, second=second.stage)
      }
