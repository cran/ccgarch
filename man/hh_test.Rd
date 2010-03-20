\name{hh.test}
\alias{hh.test}
\title{Carrying out the test of Hafner and Herwartz }
\description{
    This function computes the test statistic and the associated p-value of the 
    test for causality in conditiona variance in the CC-GARCH models. 
}
\usage{

    hh.test(dvar)

}
\arguments{
  \item{dvar}{ \eqn{(T \times N)}}
}
\value{
    A vector containing the test statistic and the associated p-value
}

\references{
    Hafner, C.M. and H. Herwartz (2006),
    \dQuote{A Lagrange Multiplier Test for Causality in Variance.}
    \emph{Economics Letters}
    \bold{93}, 137--141.
}

\seealso{
    \code{\link{nt.test}}
}

\keyword{ts}
\keyword{models}
\keyword{multivariate}

