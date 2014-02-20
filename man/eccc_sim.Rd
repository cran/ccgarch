\name{eccc.sim}
\alias{eccc.sim}
\title{Simulating an (E)CCC-GARCH(1,1) process}

\description{
    This function simulates data either from the original CCC-GARCH by Bollerslev (1990) or
    from the Extended CCC-GARCH that has non-zero off-diagonal entries in the parameter 
    matrices in the GARCH equation.  The innovations (the standardised residuals) can be 
    either a normal or student's $t$ distribution.

The dimension \eqn{(N)} is determined by the number of elements in the \eqn{\mathbf{a}} vector. 
}

\usage{
   eccc.sim(nobs, a, A, B, R, d.f=Inf, cut=1000, model)
}

\arguments{
    \item{nobs}{a number of observations to be simulated \eqn{(T)}}
    \item{a}{a vector of constants in the GARCH equation \eqn{(N \times 1)}}
    \item{A}{an ARCH parameter matrix in the GARCH equation. \eqn{\mathbf{A}} can be a diagonal matrix for the 
             original CCC-GARCH model or a full matrix for the extended model \eqn{(N \times N)}}
    \item{B}{a GARCH parameter matrix in the GARCH equation. \eqn{\mathbf{B}} can be a diagonal matrix for the 
             original CCC-GARCH model or a full matrix for the extended model \eqn{(N \times N)}}
    \item{R}{a constant conditional correlation matrix \eqn{(N \times N)}}
    \item{d.f}{the degrees of freedom parameter for the \eqn{t}-distribution}
    \item{cut}{the number of observations to be thrown away for removing initial effects of simulation}
    \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
               and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}
\value{
    A list with components:
    \item{h}{a matrix of the simulated conditional variances \eqn{(T \times N)}}
    \item{eps}{a matrix of the simulated time series with (E)CCC-GARCH process \eqn{(T \times N)}}
}

\note{
    When \code{d.f=Inf}, the innovations (the standardised residuals) follow the standard 
    normal distribution. Otherwise, they follow a student's \eqn{t}-distribution with 
    \code{d.f} degrees of freedom equal.
    
    When \code{model="diagonal"}, only the diagonal entries in \eqn{\mathbf{A}} and \eqn{\mathbf{B}} are used. If the 
    ARCH and GARCH matrices do not satisfy the stationarity condition, the simulation is 
    terminated.
}


\references{
  Bollerslev, T. (1990),
  \dQuote{Modeling the Coherence in Short-Run Nominal Exchange 
          Rates: A Multivariate Generalized ARCH Approach},
  \emph{Review of Economics and Statistics},
  \bold{72}, 498--505.
  
  Nakatani, T. and T. Ter\"asvirta (2009),
  \dQuote{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model},
  \emph{Econometrics Journal}, \bold{12}, 147--163.

  Nakatani, T. and T. Ter\"asvirta (2008),
  \dQuote{Appendix to \emph{Testing for Volatility Interactions in the Constant Conditional Correlation GARCH Model}}
  Department of Economic Statistics, Stockholm School of Economics, 
  available at \url{http://swopec.hhs.se/hastef/abs/hastef0649.htm}.
}

\seealso{
    \code{\link{dcc.sim}}, 
    \code{\link{stcc.sim}}
}

\examples{
# Simulating data from the original CCC-GARCH(1,1) process
nobs <- 1000; cut <- 1000; nu <- 10
a <- c(0.003, 0.005, 0.001)
A <- diag(c(0.2,0.3,0.15))
B <- diag(c(0.79, 0.6, 0.8))
R <- matrix(c(1.0, 0.4, 0.3, 0.4, 1.0, 0.12, 0.3, 0.12, 1.0),3,3) 

\dontrun{
ccc.data   <- eccc.sim(nobs,a, A, B, R, model="diagonal")
ccc.data.t <- eccc.sim(nobs,a, A, B, R, d.f=nu, model="diagonal")
}
}

\keyword{ts}
\keyword{datagen}
\keyword{models}
\keyword{multivariate}

