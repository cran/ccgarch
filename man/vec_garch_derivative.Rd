\name{vec.garch.derivative}
\alias{vec.garch.derivative}
\title{ Computing partial derivatives of a vector GARCH(1, 1) equation}

\description{
    This function computes partial derivatives of a vector GARCH(1, 1) equation with respect to 
    its parameters.
}

\usage{
   vec.garch.derivative(dvar, B, h)
}

\arguments{
    \item{dvar}{a matrix of the data used for estimating an ECCC or DCC GARCH model  \eqn{(T \times N)}}
    \item{B}{a GARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
    \item{h}{a matrix of conditional variances \eqn{(T \times N)}}
}

\value{
    a vector of partial derivatives \eqn{(T \times N*npar.h)}
}

\references{
  Nakatani, T. and T. Ter\"asvirta (2009),
  \dQuote{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model},
  \emph{Econometrics Journal}, \bold{12}, 147--163.

  Nakatani, T. and T. Ter\"asvirta (2008),
  \dQuote{Appendix to \emph{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model}}
  Department of Economic Statistics, Stockholm School of Economics, 
  available at \url{http://swopec.hhs.se/hastef/abs/hastef0649.htm}.
}

\keyword{ts}
\keyword{models}
\keyword{multivariate}


