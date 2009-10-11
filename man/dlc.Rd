\name{dlc}
\alias{dlc}
\title{Various partial derivatives of the DCC part of the log-likelihood function }
\description{
    This function computes various analytical derivatives of the second stage log-likelihood 
    function (the DCC part) of the (E)DCC-GARCH model. 
}
\usage{

    dlc(dcc.para, B, u, h, model)

}
\arguments{
  \item{dcc.para}{the estimates of the (E)DCC parameters  \eqn{(2 \times 1)}}
  \item{B}{the estimated GARCH parameter matrix  \eqn{(N \times N)}}
  \item{u}{a matrix of the used for estimating the (E)DCC-GARCH model  \eqn{(T \times N)}}
  \item{h}{a matrix of the estimated conditional variances  \eqn{(T \times N)}}
  \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
               and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}
\value{
    a list with components:
    \item{dlc}{the gradient of the DCC log-likelihood function w.r.t. the DCC parameters \eqn{(T \times 2)}}
    \item{dvecP}{the partial derivatives of the DCC matrix, \eqn{P_t} w.r.t. the DCC parameters \eqn{(T \times N^{2})}}
    \item{dvecQ}{the partial derivatives of the \eqn{Q_t} matrices w.r.t. the DCC parameters \eqn{(T \times N^{2})}}
    \item{d2lc}{the Hessian of the DCC log-likelihood function w.r.t. the DCC parameters \eqn{(T \times 4)}}
    \item{dfdwd2lc}{the cross derivatives of the DCC log-likelihood function \eqn{(T \times npar.h+2)} 
                    \eqn{npar.h} stand for the number of parameters in the GARCH part, \eqn{npar.h = 3N} 
                    for \code{"diagonal"} and \eqn{npar.h = 2N^{2}+N} for \code{"extended"}.}
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

\keyword{ts}
\keyword{models}
\keyword{multivariate}

