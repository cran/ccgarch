\name{dlv}
\alias{dlv}
\title{Gradient of the GARCH part of the log-likelihood function of the DCC GARCH model }
\description{
This function returns the analytical partial derivatives of the volatility part of 
the log-likelihood function of the DCC-GARCH model.
The function is called from "dcc.results".
}
\usage{
    dlv(u, a, A, B, model)
}
\arguments{
  \item{u}{a matrix of the observed residuals \eqn{(T \times N)} }
  \item{a}{a vector of the constants in the volatility part \eqn{(N \times 1)}}
  \item{A}{an ARCH parameter matrix \eqn{(N \times N)}}
  \item{B}{a GARCH parameter matrix \eqn{(N \times N)}}
  \item{model}{a character string describing the model. "diagonal" for the diagonal model 
               and "extended" for the extended (full ARCH and GARCH parameter matrices) model}
}
\value{
    A matrix of partial derivatives. \eqn{(T \times npar.h)} where \eqn{npar.h} stand for 
    the number of parameters in the GARCH part, \eqn{npar.h = 3N} for "diagonal" and 
    \eqn{npar.h = 2N^{2}+N} for "extended".
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

    Hafner, C.M. and H. Herwartz (2008),
    \dQuote{Analytical Quasi Maximum Likelihood Inference in Multivariate Volatility Models}
    \emph{Metrika}
    \bold{67}, 219--239.
}

\seealso{
    \code{\link{dcc.estimation}}
}


\keyword{ts},
\keyword{models},
\keyword{multivariate}

