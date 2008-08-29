\name{analytical.grad}
\alias{analytical.grad}
\title{Analytical gradient of the log-likelihood function of the (E)CCC-GARCH(1,1) model}

\description{
This function returns the analytical gradient of the log-likelihood function of the (E)CCC-GARCH(1,1) model.
}

\usage{
   analytical.grad(a, A, B, R, u, model)
}

\arguments{
  \item{a}{a vector of constants \eqn{(N \times 1)}}
  \item{A}{an ARCH parameter matrix \eqn{(N \times N)}}
  \item{B}{a GARCH parameter matrix \eqn{(N \times N)}}
  \item{R}{a constant conditional correlation \eqn{(N \times N)}}
  \item{u}{a matrix of the observed residuals \eqn{(T \times N)}}
  \item{model}{a character string describing the model. "diagonal" for the diagonal model 
               and "extended" for the extended (full ARCH and GARCH parameter matrices) model}
}
\value{
    a \eqn{(npar \times T)} matrix of gradients
}

\note{
    In the output, each column (not row) corresponds to the gradient at observation \eqn{t}. 
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

