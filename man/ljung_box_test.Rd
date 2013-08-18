\name{ljung.box.test}
\alias{ljung.box.test}
\title{ The Ljung-Box Test statistic }

\description{
    This function performs the Ljung-Box Test for a univariate time series.
}

\usage{
    ljung.box.test(x)
}

\arguments{
    \item{x}{a vector of variables to be tested}
}

\value{
    LB test statistics and associated p-values for lags 5, 10,..., 50.
}

\note{
    Argument \code{x} must be a vector. When \code{x} is squared
    residuals, the test is equivalent to the McLeord and Li (1983) test.  
}

\references{
    Ljung, G.M. and G.E.P. Box (1978):
    \dQuote{On a Measure of Lack of Fit in Time-Series Models},
    \emph{Biometrika}, \bold{65}, 297--303.
    
    McLeod, A.I., and W.K. Li (1983):
    \dQuote{Diagnostic checking ARMA time series models using
      squared-residual autocorrelations}, 
    \emph{Journal of Time Series Analysis}, \bold{4}, 269--273.
}

\seealso{
    \code{\link{rob.sk}}, 
    \code{\link{rob.kr}}, 
    \code{\link{jb.test}}
}

\examples{
x <- rnorm(1000)
ljung.box.test(x)           # returns the LB Test statistic
ljung.box.test(x^2)         # returns the McLeord-Li Test for no-ARCH effect
}

\keyword{htest}
\keyword{ts}
