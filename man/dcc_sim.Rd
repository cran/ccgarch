\name{dcc.sim}
\alias{dcc.sim}
\title{Simulating an (E)DCC-GARCH(1,1) process}

\description{
    This function simulates data either from the original DCC-GARCH by Engle (2002) or 
    from the Extended DCC-GARCH that has non-zero off-diagonal entries in the parameter
    matrices in the GARCH equation, with multivariate normal or student's \eqn{t} distributions. 

The dimension (\eqn{N}) is determined by the number of elements in the \eqn{a} vector. 
}

\usage{
    dcc.sim(nobs, a, A, B, R, dcc.para, d.f=Inf, cut=1000, model)
}

\arguments{
    \item{nobs}{a number of observations to be simulated (\eqn{T})}
    \item{a}{a vector of constants in the vector GARCH equation \eqn{(N \times 1)}}
    \item{A}{an ARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
    \item{B}{a GARCH parameter matrix in the vector GARCH equation \eqn{(N \times N)}}
    \item{R}{an unconditional correlation matrix \eqn{(N \times N)}}
    \item{dcc.para}{a vector of the DCC parameters \eqn{(2 \times 1)}}
    \item{d.f}{the degrees of freedom parameter for the \eqn{t}-distribution}
    \item{cut}{the number of observations to be thrown away for removing initial effects of simulation}
    \item{model}{a character string describing the model. "\code{diagonal}" for the diagonal model 
               and "\code{extended}" for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    A list with components:
    \item{z}{a matrix of random draws from \eqn{N(\mathbf{0}, \mathbf{I})}.  \eqn{(T \times N)}}
    \item{std.z}{a matrix of the standardised residuals. \eqn{\mathnormal{std.z}_{t} \sim N(0, \mathbf{R}_{t})} 
                 where \eqn{\mathbf{R}_{t}} is the DCC matrix at \eqn{t}.
             If d.f is set to a finite positive real number, \eqn{\mathbf{z}_{t} \sim t_{d.f}(0, \mathbf{R}_{t})} \eqn{(T \times N)}}
    \item{dcc}{a matrix of the simulated dynamic conditional correlations \eqn{(T \times N^2)}}
    \item{h}{a matrix of the simulated conditional variances \eqn{(T \times N)}}
    \item{eps}{a matrix of the simulated time series with DCC-GARCH process \eqn{(T \times N)}}

}

\note{
    When \code{d.f=Inf}, the innovations (the standardised residuals) follow the standard 
    normal distribution. Otherwise, they follow a student's \eqn{t}-distribution with 
    \code{d.f} degrees of freedom.

    When \code{model="diagonal"}, only the diagonal entries in A and B are used. If the 
    ARCH and GARCH matrices do not satisfy the stationarity condition, the simulation is 
    terminated.
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
    \code{\link{stcc.sim}}, 
    \code{\link{eccc.sim}}
}

\examples{

# Simulating data from the original DCC-GARCH(1,1) process
nobs <- 1000; cut <- 1000; nu <- 8
a <- c(0.003, 0.005, 0.001)
A <- diag(c(0.2,0.3,0.15))
B <- diag(c(0.75, 0.6, 0.8))
uncR <- matrix(c(1.0, 0.4, 0.3, 0.4, 1.0, 0.12, 0.3, 0.12, 1.0),3,3)
dcc.para <- c(0.01,0.98)

\dontrun{

# for normally distributed innovations
dcc.data <- dcc.sim(nobs, a, A, B, uncR, dcc.para, model="diagonal")

# for t distributed innovations
dcc.data.t <- dcc.sim(nobs, a, A, B, uncR, dcc.para, d.f=nu,
model="diagonal")
}
}


\keyword{ts}
\keyword{datagen}
\keyword{models}
\keyword{multivariate}

