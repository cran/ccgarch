\name{dcc.estimation}
\alias{dcc.estimation}
\title{ Estimating (E)DCC-GARCH model }
\description{

This function carries out the two step estimation of the (E)DCC-GARCH model and returns 
estimates, standardised residuals, the estimated volatility, and the dynamic conditional 
correlation. The details of the first and second stage estimation are also reported.

}
\usage{
    dcc.estimation(inia, iniA, iniB, ini.dcc, dvar, model)
}
\arguments{
    \item{inia}{a vector of initial values for the constants in the GARCH equation \eqn{(N \times 1)}}
    \item{iniA}{a matrix of initial values for the ARCH parameter \eqn{(N \times N)}}
    \item{iniB}{a matrix of initial values for the GARCH parameter \eqn{(N \times N)}}
    \item{ini.dcc}{a vector of initial values for the DCC parameters \eqn{(2 \times 1)}}
    \item{dvar}{a matrix of the observed residuals \eqn{(T \times N)}}
    \item{model}{a character string describing the model. "diagonal" for the diagonal model 
               and "extended" for the extended (full ARCH and GARCH parameter matrices) model}
}
\value{
    a list with components:
    \item{out}{the estimates and their standard errors}
    \item{h}{a matrix of the estimated volatilities \eqn{(T \times N)}}
    \item{DCC}{a matrix of DCC estimates \eqn{(T \times N^{2})}}
    \item{first}{the results of the first stage estimation}
    \item{second}{the results of the second stage estimation}
}

\note{
    \code{dcc.estimation} calls \code{\link{dcc.estimation1}} and \code{\link{dcc.estimation2}} for the 
    first and second stage estimation respectively. 
}

\references{
  Engle, R.F. and K. Sheppard (2001),
  \dQuote{Theoretical and Empirical Properties of Dynamic
  Conditional Correlation Multivariate GARCH.}
  \emph{Stern Finance Working Paper Series}
  {FIN}-01-027 (Revised in Dec. 2001),
  New York University Stern School of Business.
  
  Engle, R.F. (2002),
  \dQuote{Dynamic Conditional Correlation: A Simple Class of
  Multivariate Generalized Autoregressive Conditional
  Heteroskedasticity Models.}
  \emph{Journal of Business and Economic Statistics}
  \bold{20}, 339-350.
}

\seealso{
    \code{\link{dcc.estimation1}}, 
    \code{\link{dcc.estimation2}}, 
    \code{\link{loglik.dcc1}},
    \code{\link{loglik.dcc2}},
    \code{\link{vector.garch}},
    \code{\link{dcc.est}}
}

\keyword{ts},
\keyword{multivariate},
\keyword{optimize}


