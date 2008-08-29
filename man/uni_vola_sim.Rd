\name{uni.vola.sim}
\alias{uni.vola.sim}
\title{ Simulating a series with univariate GARCH(1,1) volatility}
\description{

This function simulates an univariate time series with a GARCH(1,1) volatility.

}

\usage{
    uni.vola.sim(a, nobs, d.f=Inf, cut=1000)
}

\arguments{
    \item{a}{a vector of parameters \eqn{(3 \times 1)} }
    \item{nobs}{a number of observations simulated \eqn{(T)}}
    \item{d.f}{degrees of freedom parameter for \eqn{t}-distribution}
    \item{cut}{a number of observations to be removed to minimise the initial effects}

}
\value{
    A list with components:
    \item{h}{a GARCH(1,1) volatility \eqn{(T \times 1)}}
    \item{eps}{a series of error term with the conditional variances "h" \eqn{(T \times 1)}}
}

\note{
    When "d.f=Inf", the innovations (the standardised residuals) follow the standard 
    normal distribution. Otherwise, they follow a student's \eqn{t}-distribution with 
    the degrees of freedom equal to "d.f".
}

\references{
  Bollerslev, T. (1986),
  \dQuote{Generalized Autoregressive Conditional Heteroskedasticity},
  \emph{Journal of Econometrics},
  \bold{31},
  307--327.
  
  Fiorentini, G., G. Calzolari and L. Panattoni (1996),
  \dQuote{Analytic Derivatives and the Computation of GARCH Estimates},
  \emph{Journal of Applied Econometrics},
  \bold{11},
  399--417.

}


\seealso{ \code{\link{uni.vola}} }

\examples{

nobs <- 1000
nu = 8
a = c(0.1,0.2,0.7)           #  a <- c(a constant, ARCH parameter, GARCH parameter) 
# with normal innovations
eps = uni.vola.sim(a, nobs)
# with t innovations
eps.t = uni.vola.sim(a, nobs, d.f = df)

}

\keyword{ts},
\keyword{datagen},
\keyword{models}
