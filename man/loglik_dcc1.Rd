\name{loglik.dcc1}
\alias{loglik.dcc1}
\title{ The 1st stage log-likelihood function for the (E)DCC GARCH  }

\description{
    This function returns a log-likelihood of the (E)DCC-GARCH model in the first stage estimation.
}

\usage{
   loglik.dcc1(param, dvar, model)
}

\arguments{
  \item{param}{initial values for a vector of the parameters \eqn{(npar \times 1)} }
  \item{dvar}{a matrix of the data \eqn{(T \times N)}}
  \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
               and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    the negative of the first stage log-likelihood
}

\note{
    The function is used in \code{optim} in \code{\link{dcc.estimation1}}. 
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
    \code{\link{dcc.estimation}}, 
    \code{\link{dcc.estimation1}}
}


\keyword{ts}
\keyword{models}
\keyword{multivariate}

