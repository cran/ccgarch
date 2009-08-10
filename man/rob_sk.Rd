\name{rob.sk}
\alias{rob.sk}
\title{ Computing standard and robustified skewness }

\description{
    This function computes standard and robustified skewness measures of a vector or matrix of variables.
}

\usage{
    rob.sk(x)
}

\arguments{
    \item{x}{ a vector or matrix of variables }
}

\value{
Vector of skewness and robustified skewness
}

\references{
   Kim, T-H. and H. White (2004),
   \dQuote{On More Robust Estimation of Skewness and Kurtosis}, 
   \emph{Finance Research Letters}, \bold{1}, 56--73.
}

\seealso{ 
    \code{\link{rob.kr}}, 
    \code{\link{ljung.box.test}}, 
    \code{\link{jb.test}} 
}
\examples{
x <- matrix(rnorm(1000), 100, 10)
rob.sk(x)
}

\keyword{htest}
\keyword{ts}
