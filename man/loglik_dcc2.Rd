\name{loglik.dcc2}
\alias{loglik.dcc2}
\title{ The 2nd stage log-likelihood function for the (E)DCC GARCH  }

\description{
    This function returns a log-likelihood of the (E)DCC-GARCH model in the 2nd step estimation.
}

\usage{
   loglik.dcc2(param, dvar)
}

\arguments{
  \item{param}{initial values for the DCC parameters \eqn{(2 \times 1)}}
  \item{dvar}{a matrix of the standardised residuals \eqn{(T \times N)}}
}

\value{
  the negative of the second stage log-likelihood
}

\note{
    The function is used in \code{constrOptim} in \code{\link{dcc.estimation2}}. 
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
  Multivariate GeneralizSed Autoregressive Conditional
  Heteroskedasticity Models.}
  \emph{Journal of Business and Economic Statistics}
  \bold{20}, 339--350.
}

\seealso{
    \code{\link{dcc.estimation}}, 
    \code{\link{dcc.estimation2}}
}

\keyword{ts}
\keyword{models}
\keyword{multivariate}

