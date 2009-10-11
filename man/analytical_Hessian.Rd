\name{analytical.Hessian}
\alias{analytical.Hessian}
\title{ Analytical Hessian of the (E)CCC-GARCH}

\description{
    This function computes the analytical Hessian of the log-likelihood function of the (E)CCC-GARCH model.
}

\usage{
   analytical.Hessian(a, A, B, R, u, model)
}

\arguments{
  \item{a}{a vector of constants in the vector GARCH equation \eqn{(N \times 1)}}
  \item{A}{an ARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
  \item{B}{a GARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
  \item{R}{a constant conditional correlation matrix \eqn{(N \times N)}}
  \item{u}{a matrix of the data data used for estimating the (E)CCC-GARCH(1,1) model \eqn{(T \times N)}}
  \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    a \eqn{npar \times npar} Hessian matrix of the log-likelihood function of the (E)CCC-GARCH model
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


