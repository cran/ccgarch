\name{eccc.estimation}
\alias{eccc.estimation}
\title{ Estimating an (E)CCC-GARCH model }

\description{
    This function estimates an (E)CCC-GARCH(1,1) model and returns estimates, estimated volatility and 
    various diagnostic statistics.
}

\usage{
      eccc.estimation(a, A, B, R, dvar, model, method="BFGS")
}

\arguments{
    \item{a}{initial values for constants \eqn{(N \times 1)}}
    \item{A}{initial values for an ARCH parameter matrix \eqn{(N \times N)}}
    \item{B}{initial values for a GARCH parameter matrix \eqn{(N \times N)}}
    \item{R}{initial values a constant conditional correlation matrix \eqn{(N \times N)}}
    \item{dvar}{a matrix of data used for (E)CCC-GARCH estimation \eqn{(T \times N)}}
    \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
               and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
    \item{method}{a character string specifying the optimisation method in \code{optim}. 
               There are three choices, namely, \code{Nelder-Mead}, \code{BFGS} (default) and \code{CG}.}
}

\value{
    A list with components:
    \item{out}{a \eqn{(4 \times npar)} matrix. The estimates are contained in the first row. 
               The remaining rows report standard errors based on three different methods of 
               estimating the asymptotic covariance matrix}
    \item{h}{the estimated conditional variances \eqn{(T \times N)}}
    \item{std.resid}{a matrix of the standardised residuals \eqn{(T \times N}). See \emph{Note}.}
    \item{opt}{the detailed results of the optimisation}
    \item{para.mat}{vectorised parameter estimates}
}

\note{
    The standardised residuals are calculated through dividing the
    original series by the estimated conditional standard deviations.
    See, for instance, p.303 of Bollerslev (1990) for details.
}


\references{
  Bollerslev, T. (1990), \dQuote{Modelling the Coherence in Short-run Nominal Exchange Rates: A Multivariate Generalized ARCH Model}, 
  \emph{Review of Economics and Statistics}, \bold{20}, 498--505.
  
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
\keyword{optimize}
\keyword{multivariate}
