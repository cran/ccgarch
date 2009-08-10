\name{dcc.estimation2}
\alias{dcc.estimation2}
\title{ Maximising the second stage log-likelihood function of the (E)DCC-GARCH model }
\description{
    This function carries out the second stage (DCC part) estimation of the (E)DCC-GARCH model.
}
\usage{
    dcc.estimation2(dvar, para, gradient=0)
}
\arguments{
    \item{dvar}{a matrix of the standardised residuals \eqn{(T \times N)}}
    \item{para}{a vector of the DCC parameters \eqn{(2 \times 1)}}
    \item{gradient}{a switch variable whether to use the gradient in the constraint optimisation. passed to 
                    \code{constrOptim}}
}
\value{
    a list of the estimation results. See the explanations for \code{constrOptim}.
}

\note{
    \code{dcc.estimation2} is a wrapper to \code{constrOptim}. The restrictions are 
    \eqn{\alpha + \beta \leq 1} and \eqn{\alpha, \beta \geq 0} in the
    DCC equation. 
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
    \code{constrOptim}, 
    \code{\link{dcc.estimation1}}, 
    \code{\link{dcc.estimation}}
}

\keyword{ts}
\keyword{multivariate}
\keyword{optimize}


