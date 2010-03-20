\name{nt.test}
\alias{nt.test}
\title{Carrying out the test of Nakatani and Ter\"asvirta }
\description{
    This function computes the test statistic and the associated p-value of the 
    test for causality in conditiona variance in the CC-GARCH models. 
}
\usage{
    nt.test(dvar)
}
\arguments{
  \item{dvar}{ \eqn{(T \times N)}}
}
\value{
    A matrix containing the test statistics of the standard (non-robust) test and the robust version, and the associated p-values
}

\references{
    Nakatani, T and T. Ter\"asvirta (2010),
    \dQuote{An Alternative Test for Causality in Variance in the Conditional Correlation GARCH models.}
    \emph{mimeo}, 
    Stockholm School of Economics.
}

\seealso{
    \code{\link{hh.test}}
}


\keyword{ts}
\keyword{models}
\keyword{multivariate}

