\name{uni.vola}
\alias{uni.vola}
\title{ Computing univariate GARCH(1,1) conditional variances }

\description{
    This function returns an univariate GARCH(1,1) conditional variances.
}
\usage{
   uni.vola(a,u)
}

\arguments{
  \item{a}{a vector of parameters in the GARCH(1,1) equation \eqn{(3 \times 1)}}
  \item{u}{a vector of the data \eqn{(T \times 1)}}
}

\value{
    a vector of GARCH(1,1) conditional variances \eqn{(T \times 1)}
}

\references{
  Bollerslev, T. (1986): 
  \dQuote{Generalized Autoregressive Conditional Heteroskedasticity},
  \emph{Journal of Econometrics},
  \bold{31}, 307--327.
}

\seealso{
  \code{\link{uni.vola.sim}}
}

\examples{
a <- c(0.01, 0.04, 0.95)    #  a <- c(a constant, ARCH parameter, GARCH parameter) 
u <- rnorm(1000)
h <- uni.vola(a, u)
}

\keyword{ts}
\keyword{univar}
