\name{grad.dcc.full}
\alias{grad.dcc.full}
\title{ Numerical gradient of the full log-likelihood function of the (E)DCC-GARCH model }

\description{
    This function computes numerical gradient of the full log-likelihood 
    function of the (E)DCC-GARCH(1,1) model with respect to its parameters.
}

\usage{
    grad.dcc.full(a, A, B, dcc.para, dvar, d=1e-5, model)
}

\arguments{
    \item{a}{a constant vector in the vector GARCH equation \eqn{(N \times 1)}}
    \item{A}{an ARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
    \item{B}{a GARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
    \item{dcc.para}{a vector of the DCC parameters \eqn{(2 \times 1)}}
    \item{dvar}{a matrix of the data used for estimating the (E)DCC-GARCH model \eqn{(T \times N)}}
    \item{d}{a step size for computing numerical gradient}
    \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
                 and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    a matrix of partial derivatives \eqn{(T \times npar)}
}

\note{
    this function is currently not in use.
}

\keyword{ts}
\keyword{models}
\keyword{multivariate}


