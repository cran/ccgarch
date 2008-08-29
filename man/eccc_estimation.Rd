\name{eccc.estimation}
\alias{eccc.estimation}
\title{ Estimating an (E)CCC-GARCH model }

\description{
This function estimates an (E)CCC-GARCH(1,1) model and returns estimates, estimated volatility and 
various diagnostic statistics.
}

\usage{
      eccc.estimation(a, A, B, R, dvar, model)
}

\arguments{
    \item{a}{a vector of constants \eqn{(N \times 1)}}
    \item{A}{an ARCH parameter matrix \eqn{(N \times N)}}
    \item{B}{a GARCH parameter matrix \eqn{(N \times N)}}
    \item{R}{a CCC matrix \eqn{(N \times N)}}
    \item{dvar}{a matrix of data \eqn{(T \times N)}}
    \item{model}{a character string describing the model. "diagonal" for the diagonal model 
               and "extended" for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    A list with components:
    \item{out}{a \eqn{(4 \times npar)} matrix. The estimates are contained in the first row. 
               The remaining rows report standard errors based on three different methods of 
               estimating the asymptotic covariance matrix}
    \item{h}{the estimated volatility \eqn{(T \times N)}}
    \item{opt}{the detailed results of the optimisation}
    \item{para.mat}{vectorised parameter estimates}
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
\keyword{models},
\keyword{multivariate}
