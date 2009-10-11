\name{d2lv}
\alias{d2lv}
\title{Hessian of the DCC log-likelihood function}
\description{
    This function returns the analytical Hessian of the volatility part of the DCC log-likelihood function.
}
\usage{

    d2lv(u, B, h, model)

}
\arguments{
  \item{u}{a matrix of the data data used for estimating the (E)DCC-GARCH(1,1) model \eqn{(T \times N)}}
  \item{B}{a GARCH parameter matrix \eqn{(N \times N)}}
  \item{h}{a matrix of the conditional variances \eqn{(T \times N)}}
  \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
               and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}
\value{
    the Hessian of the volatility part of the DCC log-likelihood function \eqn{(T \times N^{2})}
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

