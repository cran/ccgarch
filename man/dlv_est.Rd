\name{dlv.est}
\alias{dlv.est}
\title{Gradient of the GARCH part of the log-likelihood function of an (E)DCC GARCH model }
\description{
    This function returns the gradient of the volatility part of the log-likelihood function of the DCC.
}

\usage{
    dlv.est(par, dvar, model)
}

\arguments{
  \item{par}{a vector of the parameters in the vector GARCH equation}
  \item{dvar}{a matrix of the data used for estimating an (E)DCC-GARCH model \eqn{(T \times N)}}
  \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
               and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}
\value{
    A vector of the gradient. \eqn{(3N \times 1)} for \code{"diagonal"} 
    and \eqn{(2N^{2}+N \times 1)} for \code{"extended"}.
}

\note{
    The function can be called from \code{optim} in \code{\link{dcc.estimation1}}. For obtaining 
    the gradient for all \eqn{t}, use \code{\link{dlv}} instead.
}

\references{
    Engle, R.F. and K. Sheppard (2001),
    \dQuote{Theoretical and Empirical Properties of Dynamic
    Conditional Correlation Multivariate GARCH.}
    \emph{Stern Finance Working Paper Series}
    FIN-01-027 (Revised in Dec. 2001),
    New York University Stern School of Business.

    Engle, R.F. (2002),
    \dQuote{Dynamic Conditional Correlation: A Simple Class of
    Multivariate Generalized Autoregressive Conditional
    Heteroskedasticity Models.}
    \emph{Journal of Business and Economic Statistics}
    \bold{20}, 339--350.

    Hafner, C.M. and H. Herwartz (2008),
    \dQuote{Analytical Quasi Maximum Likelihood Inference in Multivariate Volatility Models.}
    \emph{Metrika}
    \bold{67}, 219--239.
}

\seealso{
    \code{\link{dcc.estimation1}},
    \code{\link{dlv}}
}


\keyword{ts}
\keyword{models}
\keyword{multivariate}

