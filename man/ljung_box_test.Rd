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
    \item{x}{a vector of variables }
}

\value{
    LB test statistics and associated p-values for lags 5, 10,..., 50.
}

\note{
    This function is not valid for multiple series.
}

\references{
    Ljung, G.M. and G. E. P. Box (1978):
    \dQuote{On a Measure of Lack of Fit in Time-Series Models},
    \emph{Biometrika}, \bold{65}, 297-303.
}

\seealso{
    \code{\link{rob.sk}}, 
    \code{\link{rob.kr}}, 
    \code{\link{jb.test}}
}

\examples{
x <- rnorm(1000)
ljung.box.test(x)           # returns the LB Test statistic
ljung.box.test(x^2)         # returns the McLord-Li Test for no-ARCH effect
}

\keyword{htest},
\keyword{ts}
