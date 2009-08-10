\name{tr.func}
\alias{tr.func}
\title{ Logistic transition function }

\description{
    This function computes values from a Logistic transition function.
}

\usage{
    tr.func(tr.par, tr.var)
}

\arguments{
    \item{tr.par}{a vector of parameters \eqn{(2 \times 1)}}
    \item{tr.var}{a vector of transition variable \eqn{(T \times 1)}}
}

\value{
    a vector of transition function \eqn{(T \times 1)}
}

\note{
    this function is used in \code{\link{stcc.sim}}
}


\references{
  Ter\"asvirta, T. (1994):
  \dQuote{Specification, Estimation, and Evaluation of Smooth Transition Autoregressive Models},
  \emph{Journal of the American Statistical Association},
  \bold{89}, 208--218.
}

\seealso{
    \code{\link{stcc.sim}}
}

\keyword{ts}


