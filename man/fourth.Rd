\name{fourth}
\alias{fourth}
\title{ Fourth-order moment condition for the vector GARCH equation }

\description{
    This function computes the fourth-order moment condition for the vector GARCH equation 
    in the (E)CCC-GARCH models.
}

\usage{
    fourth(A, B, R)
}
\arguments{
    \item{A}{an ARCH parameter matrix \eqn{(N \times N)}}
    \item{B}{a GARCH parameter matrix \eqn{(N \times N)}}
    \item{R}{a constant conditional correlation matrix \eqn{(N \times N)}}
}

\value{
    a scalar. If strictly less than unity, the condition is satisfied. 
}

\references{
   He, C. and T. Ter\"asvirta (2004):
   \dQuote{An Extended Constant Conditional Correlation GARCH model and its Fourth-moment Structure}, 
   \emph{Econometric Theory}, \bold{20}, 904--926.

  Nakatani, T. and T. Ter\"asvirta (2009),
  \dQuote{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model},
  \emph{Econometrics Journal}, \bold{12}, 147--163.

  Nakatani, T. and T. Ter\"asvirta (2008),
  \dQuote{Appendix to \emph{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model}}
  Department of Economic Statistics, Stockholm School of Economics, 
  available at \url{http://swopec.hhs.se/hastef/abs/hastef0649.htm}.
}

\seealso{
    \code{\link{stationarity}}
}

\keyword{algebra}
\keyword{multivariate}

