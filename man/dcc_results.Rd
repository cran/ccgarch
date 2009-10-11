\name{dcc.results}
\alias{dcc.results}
\title{Computing robust standard errors of the estimates in the (E)DCC-GARCH model }
\description{
    This function computes the robust standard errors of the estimates of a DCC-GARCH model.
}
\usage{

    dcc.results(u, garch.para, dcc.para, h, model)

}
\arguments{
    \item{u}{a matrix of the data used for estimating the (E)DCC-GARCH model \eqn{(T \times N)}}
    \item{garch.para}{a vector of the estimates of the volatility parameters}
    \item{dcc.para}{a vector of the estimates of the DCC parameters \eqn{(2 \times 1)}}
    \item{h}{a matrix of the estimated conditional variances \eqn{(T \times N)}}
    \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
               and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}
\value{
    A matrix with the estimates in the first row, and the standard errors in the second row.
}

\note{
   \code{dcc.results} is called from \code{\link{dcc.estimation}}.    
    When \code{model="diagonal"}, only the diagonal entries in A and B are used.
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
}

\seealso{
    \code{\link{dcc.estimation}}
}

\keyword{ts}
\keyword{models}
\keyword{multivariate}

