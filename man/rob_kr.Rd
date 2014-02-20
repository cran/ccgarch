\name{rob.kr}
\alias{rob.kr}
\title{ Computing standard and robustified excess kurtosis }

\description{
    This function computes standard and robustified excess kurtosis of a vector or matrix of variables.
}

\usage{
    rob.kr(x)
}

\arguments{
    \item{x}{ vector or matrix of variables}
}

\value{
    Vector of excess kurtosis and robustified excess kurtosis
}

\references{
    Kim, T-H. and H. White (2004),
    \dQuote{On More Robust Estimation of Skewness and Kurtosis}, 
    \emph{Finance Research Letters}, \bold{1}, 56--73.
}

\seealso{
    \code{\link{rob.sk}}, 
    \code{\link{ljung.box.test}}, 
    \code{\link{jb.test}} 
}

\examples{
x <- matrix(rnorm(1000), 100, 10)
rob.kr(x)
}

\keyword{htest}
\keyword{ts}
