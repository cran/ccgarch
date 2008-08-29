\name{p.mat}
\alias{p.mat}
\title{ Re-arranging a vector into parameter matrices }
\description{
A utility function that re-arranges a vector of parameters into parameter matrices in the CC-GARCH(1,1) model. 
}

\usage{
   p.mat(para, model, ndim)
}
\arguments{
  \item{para}{a vector of parameters to be re-arranged into matrices of parameters }
  \item{model}{a character string describing the model. "diagonal" for the diagonal model 
               and "extended" for the extended (full ARCH and GARCH parameter matrices) model}
  \item{ndim}{the number of dimension of the model}
}
\value{
    A list with components:
    \item{a}{a vector of constants in the vector GARCH equation}
    \item{A}{an ARCH parameter matrix}
    \item{B}{a GARCH parameter matrix}
    \item{R}{a constant conditional correlation matrix}
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

\keyword{array},
\keyword{utilities}


