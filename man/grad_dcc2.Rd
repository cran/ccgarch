\name{grad.dcc2}
\alias{grad.dcc2}
\title{ Numerical gradient of the DCC part of the log-likelihood function}

\description{
    This function computes numerical gradient of the second stage log-likelihood function 
    of the DCC-GARCH model w.r.t. its parameters.
}

\usage{
   grad.dcc2(param, dvar, d=1e-5)
}

\arguments{
    \item{param}{a vector of DCC parameters \eqn{(2 \times 1)}}
    \item{dvar}{a matrix of the standardised residuals \eqn{(T \times N)}}
    \item{d}{a step size for computing numerical gradient}
}

\value{
    a vector of partial derivatives \eqn{(2 \times 1)}
}

\note{
    The function is used only from \code{\link{dcc.estimation2}} when the switch variable 
    \code{gradient!=0}. \code{dlc$dlc} in \code{\link{dlc}} contains the analytical gradient of 
    the second stage log-likelihood function. 
}

\seealso{
    \code{\link{dlc}}, 
    \code{\link{dcc.estimation}}
}

\keyword{ts}
\keyword{models}
\keyword{multivariate}


