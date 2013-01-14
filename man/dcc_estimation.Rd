\name{dcc.estimation}
\alias{dcc.estimation}
\title{ Estimating an (E)DCC-GARCH model }
\description{
    This function carries out the two step estimation of the (E)DCC-GARCH model and returns 
    estimates, standardised residuals, the estimated conditional variances, and the dynamic conditional 
    correlations.
}
\usage{
    dcc.estimation(inia, iniA, iniB, ini.dcc, dvar, model,
method="BFGS", gradient=1, message=1)
}
\arguments{
    \item{inia}{a vector of initial values for the constants in the GARCH equation 
            \code{length(inia)=N}}
    \item{iniA}{a matrix of initial values for the ARCH parameter matrix
      \eqn{(N \times N)}}
    \item{iniB}{a matrix of initial values for the GARCH parameter
      matrix \eqn{(N \times N)}}
    \item{ini.dcc}{a vector of initial values for the DCC parameters \eqn{(2 \times 1)}}
    \item{dvar}{a matrix of the data \eqn{(T \times N)}}
    \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
               and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
    \item{method}{a character string specifying the optimisation method in \code{optim}. 
               There are three choices, namely, \code{Nelder-Mead},
	       \code{BFGS} (default) and \code{CG}.}
   \item{gradient}{a switch variable that determines the optimisation
     algorithm in the second stage optimisation. If \code{gradient=0}
     \code{Nelder-Mead} is invokded. Otherwise \code{BFGS} is used (default).}
   \item{message}{a switch variable to turn off the display of the message when 
                  the estimation is completed. If \code{message=0}, the message 
                  is suppressed. Otherwise, the message is displayed (default)}
}
\value{
    a list with components:
    \item{out}{the parameter estimates and their standard errors}
    \item{loglik}{the value of the log-likelihood at the estimates}
    \item{h}{a matrix of the estimated conditional variances \eqn{(T
	\times N)}}
    \item{DCC}{a matrix of the estimated dynamic conditional correlations \eqn{(T \times N^{2})}}
    \item{std.resid}{a matrix of the standardised residuals \eqn{(T \times N}). See \emph{Note}.}
    \item{first}{the results of the first stage estimation}
    \item{second}{the results of the second stage estimation}
}




\note{
    The standardised residuals are calculated by dividing the original
    series \code{dvar} by the estimated conditional standard deviations \code{sqrt(h)}. 
    See Engle (2002), in particular the equations (2) and (14), for details.
    
    The details of the first and second stage estimation are also saved in \code{first} 
    and \code{second}, respectively.
    
    The switch variable \code{simulation} is useful when one uses \code{dcc.estimation} for simulation. 
    It supresses the display of the completion message. 
}

\references{
  Engle, R.F. and K. Sheppard (2001),
  \dQuote{Theoretical and Empirical Properties of Dynamic
  Conditional Correlation Multivariate GARCH.}
  \emph{Stern Finance Working Paper Series}
  FIN-01-027 (Revised in Dec. 2001),
  New York University Stern School of Business.
  
  Engle, R.F. (2002),
  \dQuote{Dynamic Conditional Correlation: A Simple Class of
  Multivariate Generalized Autoregressive Conditional
  Heteroskedasticity Models.}
  \emph{Journal of Business and Economic Statistics}
  \bold{20}, 339--350.
}

\seealso{
    \code{\link{dcc.sim}}
}

\examples{
# Simulating data from the original DCC-GARCH(1,1) process
  nobs <- 1000; cut <- 1000
  a <- c(0.003, 0.005, 0.001)
  A <- diag(c(0.2,0.3,0.15))
  B <- diag(c(0.75, 0.6, 0.8))
  uncR <- matrix(c(1.0, 0.4, 0.3, 0.4, 1.0, 0.12, 0.3, 0.12, 1.0),3,3)
  dcc.para <- c(0.01,0.98)
  dcc.data <- dcc.sim(nobs, a, A, B, uncR, dcc.para, model="diagonal")

\dontrun{
# Estimating a DCC-GARCH(1,1) model
  dcc.results <- dcc.estimation(inia=a, iniA=A, iniB=B, ini.dcc=dcc.para, 
        dvar=dcc.data$eps, model="diagonal")

# Parameter estimates and their robust standard errors
  dcc.results$out
}
}

\keyword{ts}
\keyword{multivariate}
\keyword{optimize}
\keyword{models}


