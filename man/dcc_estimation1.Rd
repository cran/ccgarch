\name{dcc.estimation1}
\alias{dcc.estimation1}
\title{ Maximising the first stage log-likelihood function of the (E)DCC-GARCH model}
\description{
    This function carries out the first stage (volatility part) estimation of the (E)DCC-GARCH model.
}
\usage{
    dcc.estimation1(dvar, a, A, B, model, method="BFGS")
}

\arguments{
    \item{dvar}{a matrix of the data used for estimating the (E)DCC-GARCH(1,1) model \eqn{(T \times N)}}
    \item{a}{a vector of constants in the vector GARCH equation \eqn{(N \times 1)}}
    \item{A}{an ARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
    \item{B}{a GARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
    \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
                 and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
    \item{method}{a character string specifying the optimisation method in \code{optim}. 
               There are three choices, namely, \code{"Nelder-Mead"}, \code{"BFGS"} (default) and \code{"CG"}.}
}
\value{
    a list of the estimation results. See the explanations in \code{optim}.
}
\references{
  Engle, R.F. and K. Sheppard (2001),
  \dQuote{Theoretical and Empirical Properties of Dynamic
  Conditional Correlation Multivariate GARCH.}
  \emph{Stern Finance Working Paper Series}
  FIN-01-027 (Revised in Dec. 2001),
  New York University Stern School of Business.
  
  Engle, R.F. (2002),
  \dQuote{Dynamic Conditional Correlation: A Simple Class of
  Multivariate Generalized Autoregressive Conditional
  Heteroskedasticity Models.}
  \emph{Journal of Business and Economic Statistics}
  \bold{20}, 339--350.
}
\seealso{
    \code{optim}, 
    \code{\link{dcc.estimation2}}, 
    \code{\link{dcc.estimation}}
}
\keyword{ts}
\keyword{multivariate}
\keyword{optimize}
