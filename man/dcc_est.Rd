\name{dcc.est}
\alias{dcc.est}
\title{Dynamic conditional correlations }
\description{
    This function returns dynamic conditional correlations based on the parameters specified.
}
\usage{

    dcc.est(dvar, param)

}
\arguments{
  \item{dvar}{a matrix of the standardised residuals \eqn{(T \times N)}}
  \item{param}{a vector of the DCC parameters \eqn{(2 \times 1)}}
}
\value{
    a list with components:
    \item{DCC}{a matrix of the dynamic conditional correlations \eqn{(T \times N^{2})}}
    \item{Q}{a matrix of the \eqn{\mathbf{Q}_{t}} \eqn{(T \times N^{2})}}
}

\note{
    a constant matrix \eqn{\mathbf{Q}} in the DCC equation is computed by \eqn{\mathbf{Q}=cov(dvar)}.
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

\keyword{ts}
\keyword{models}
\keyword{multivariate}

