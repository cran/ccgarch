\name{vdR}
\alias{vdR}
\title{ Computing partial derivatives of the CCC matrix }

\description{
    This function computes partial derivatives of the CCC 
    matrix with respect to its correlation coefficients.
}

\usage{
   vdR(n)
}

\arguments{
    \item{n}{ the number of dimension of the model }
}

\value{
    a matrix of zeros and ones \eqn{( (N(N-1))/2 \times N^{2} )}
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
\keyword{multivariate}
