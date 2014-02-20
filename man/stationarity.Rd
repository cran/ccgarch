\name{stationarity}

\alias{stationarity}

\title{ The stationarity condition in Extended CC-GARCH models  }

\description{
    A utility function that checks if the two parameter matrices in a vector GARCH model 
    satisfy the stationarity condition. 
}
\usage{ 
   stationarity(A,B)
}

\arguments{
  \item{A}{an ARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
  \item{B}{a GARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
}

\value{
    a scalar. If strictly less than unity, the condition is satisfied. 
}

\references{
   He, C. and T. Ter\"asvirta (2004):
   \dQuote{An Extende Constant Conditional Correlation GARCH model and its Fourth-moment Structure}, 
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
  \code{\link{fourth}}
}


\keyword{ts}
\keyword{multivariate}
