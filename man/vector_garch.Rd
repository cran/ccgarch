\name{vector.garch}
\alias{vector.garch}
\title{ A vector GARCH(1,1) volatility }

\description{
This function computes a vector GARCH(1,1) volatility.
}

\usage{
    vector.garch(dvar, a, A, B)
}

\arguments{
  \item{dvar}{a matrix of the observed residuals \eqn{(T \times N)}}
  \item{a}{a vector of constants \eqn{(N \times 1)}}
  \item{A}{an ARCH parameter matrix \eqn{(N \times N)}}
  \item{B}{a GARCH parameter matrix \eqn{(N \times N)}}
}

\value{
    a matrix of volatilities \eqn{(T \times N)}
}

\references{
  Nakatani, T. and T. er\"{a}svirta (2008),
  \dQuote{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model},
  \emph{Econometrics Journal}, forthcoming.

  Nakatani, T. and T. Ter\"{a}svirta (2008),
  \dQuote{Appendix to \emph{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model}}
  Department of Economic Statistics, Stockholm School of Economics, 
  available at \url{http://swopec.hhs.se/hastef/abs/hastef0649.htm}.
}


\keyword{ts},
\keyword{models}

