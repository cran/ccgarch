\name{loglik.eccc}
\alias{loglik.eccc}
\title{ The log-likelihood function of the (E)CCC-GARCH model }

\description{
    This function computes a log-likelihood of the (E)CCC-GARCH(1,1) model.
}

\usage{
   loglik.eccc(param, dvar, model)
}

\arguments{
    \item{param}{a vector of all the parameters in the (E)CCC-GARCH model}
    \item{dvar}{a matrix of the data used for estimating the (E)DCC-GARCH model \eqn{(T \times N)}}
    \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
                 and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    the negative of the (E)CCC-GARCH log-likelihood
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

