\name{loglik.dcc}
\alias{loglik.dcc}
\title{ The log-likelihood function for the (E)DCC GARCH model  }

\description{
    This function returns a log-likelihood of the (E)DCC-GARCH model.
}

\usage{
   loglik.dcc(param, dvar, model)
}

\arguments{
  \item{param}{a vector of all the parameters in the (E)DCC-GARCH model }
  \item{dvar}{a matrix of the data used for estimating the (E)DCC-GARCH model \eqn{(T \times N)}}
  \item{model}{a character string describing the model. "diagonal" for the diagonal model 
               and "extended" for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    the negative of the full log-likelihood of the (E)DCC-GARCH model
}

\note{
    \code{param} must be made by stacking all the parameter matrices.
}

\references{
  Robert F. Engle and Kevin Sheppard (2001),
  \dQuote{Theoretical and Empirical Properties of Dynamic
  Conditional Correlation Multivariate GARCH.}
  \emph{Stern Finance Working Paper Series}
  FIN-01-027 (Revised in Dec. 2001),
  New York University Stern School of Business.
  
  Robert F. Engle (2002),
  \dQuote{Dynamic Conditional Correlation: A Simple Class of
  Multivariate Generalised Autoregressive Conditional
  Heteroskedasticity Models.}
  \emph{Journal of Business and Economic Statistics}
  \bold{20}, 339--350.
}


\examples{
\dontrun{
# Simulating data from the original DCC-GARCH(1,1) process
  nobs <- 1000; cut <- 1000
  a <- c(0.003, 0.005, 0.001)
  A <- diag(c(0.2,0.3,0.15))
  B <- diag(c(0.75, 0.6, 0.8))
  uncR <- matrix(c(1.0, 0.4, 0.3, 0.4, 1.0, 0.12, 0.3, 0.12, 1.0),3,3)
  dcc.para <- c(0.01,0.98)
  dcc.data <- dcc.sim(nobs, a, A, B, uncR, dcc.para, model="diagonal")

# Estimating a DCC-GARCH(1,1) model
  dcc.results <- dcc.estimation(inia=a, iniA=A, iniB=B, ini.dcc=dcc.para, 
    dvar=dcc.data$eps, model="diagonal")

# Parameter estimates and their robust standard errors
  dcc.results$out

# Computing the value of the log-likelihood at the estimates
  loglik.dcc(dcc.results$out[1,], dcc.data$eps, model="diagonal")
}
}


\keyword{ts}
\keyword{models}
\keyword{multivariate}

