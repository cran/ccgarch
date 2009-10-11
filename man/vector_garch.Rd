\name{vector.garch}
\alias{vector.garch}
\title{ A vector GARCH(1,1) conditional variances }

\description{
    This function computes a vector GARCH(1,1) conditional variances.
}

\usage{
    vector.garch(dvar, a, A, B)
}

\arguments{
  \item{dvar}{a matrix of the data, used as epsilon \eqn{(T \times N)}}
  \item{a}{initial values for constants in the vector GARCH equation \eqn{(N \times 1)}}
  \item{A}{initial values for an ARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
  \item{B}{initial values for a GARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
}

\value{
    a matrix of conditional variances \eqn{(T \times N)}
}

\references{
  Nakatani, T. and T. Ter\"asvirta (2009),
  \dQuote{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model},
  \emph{Econometrics Journal}, \bold{12}, 147-163.

  Nakatani, T. and T. Ter\"asvirta (2008),
  \dQuote{Appendix to \emph{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model}}
  Department of Economic Statistics, Stockholm School of Economics, 
  available at \url{http://swopec.hhs.se/hastef/abs/hastef0649.htm}.
}


\keyword{ts}
\keyword{models}

