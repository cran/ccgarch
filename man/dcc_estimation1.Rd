\name{dcc.estimation1}
\alias{dcc.estimation1}
\title{ Estimating the first stage (E)DCC-GARCH }
\description{
This function carries out the first stage (volatility part) estimation of the (E)DCC-GARCH model.
}

\usage{
    dcc.estimation1(dvar, a, A, B, model)
}

\arguments{
    \item{dvar}{a matrix of the observed residuals \eqn{(T \times N)}}
    \item{a}{a vector of constants \eqn{(N \times 1)}}
    \item{A}{an ARCH parameter matrix \eqn{(N \times N)}}
    \item{B}{a GARCH parameter matrix \eqn{(N \times N)}}
    \item{model}{a character string describing the model. "diagonal" for the diagonal model 
                 and "extended" for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    a list of the estimation results. See the explanations in \code{optim}.
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
    \code{optim}, 
    \code{\link{dcc.estimation2}}, 
    \code{\link{dcc.estimation}}
}


\keyword{ts},
\keyword{multivariate},
\keyword{optimize}


