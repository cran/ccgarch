\name{jb.test}
\alias{jb.test}
\title{ The Lomnicki-Jarque-Bera Test of normality (JB test) }

\description{
    This function performs the Lomnicki-Jarque-Bera Test of normality and returns test statistics and 
    associated p-values.
}

\usage{
    jb.test(x)
}

\arguments{
    \item{x}{ a vector or matrix of variables to be tested}
}

\value{
    Vector of test statistics and p-value
}
\references{
    Jarque, C.M. and A.K. Bera (1987),
    \dQuote{A Test for Normality of Observations and Regression Residuals},
    \emph{International Statistical Review}, \bold{55}, 163--172.

    Lomnicki, Z.A. (1961),
    \dQuote{Tests for Departure from Normality in the Case of Linear Stochastic Processes},
    \emph{Metrika}, \bold{4}, 37--62.
}

\seealso{
    \code{\link{rob.sk}}, 
    \code{\link{rob.kr}}, 
    \code{\link{ljung.box.test}}
}

\examples{
    # for a vector
    x <- rnorm(1000)
    jb.test(x)

    # for a matrix
    X <- matrix(rnorm(10000), 5000,2)
    jb.test(X)
}

\keyword{htest}
\keyword{ts}
