\name{stcc.sim}
\alias{stcc.sim}
\title{Simulating Data from an STCC-GARCH$(1,1)$ process}

\description{
    This function simulates data either from the original STCC-GARCH by Silvennoinen and 
    Ter\"asvirta (2005) or from the Extended STCC-GARCH that has non-zero off-diagonal 
    entries in the parameter matrices in the GARCH equation, with multivariate normal or 
    student's \eqn{t} distribution. 

    The dimension \eqn{(N)} is determined by the number of elements in the \eqn{\mathbf{a}} vector. 
}

\usage{
    stcc.sim(nobs, a, A, B, R1, R2, tr.par, st.par, d.f=Inf, 
             cut=1000, model)
}

\arguments{
    \item{nobs}{a number of observations to be simulated \eqn{(T)}}
    \item{a}{a vector of constants in the vector GARCH equation \eqn{(N \times 1)}}
    \item{A}{an ARCH parameter matrix in the vector GARCH equation. \eqn{(N \times N)}}
    \item{B}{a GARCH parameter matrix in the vector GARCH equation. \eqn{(N \times N)}}
    \item{R1}{a conditional correlation matrix in regime 1 \eqn{(N \times N)}}
    \item{R2}{a conditional correlation matrix in regime 2 \eqn{(N \times N)}}
    \item{tr.par}{a vector of scale and location parameters in the transition function \eqn{(2 \times 1)}}
    \item{st.par}{a vector of parameters for the GARCH(1,1) transition variable \eqn{(3 \times 1)}}
    \item{d.f}{the degrees of freedom parameter for the \eqn{t}-distribution}
    \item{cut}{the number of observations to be thrown away for removing initial effects of simulation}
    \item{model}{a character string describing the model. \code{"diagonal"} for the diagonal model 
               and \code{"extended"} for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    A list with components:
    \item{h}{a matrix of conditional variances \eqn{(T \times N)}}
    \item{eps}{a matrix of time series with DCC-GARCH process \eqn{(T \times N)}}
    \item{tr.var}{a vector of the transition variable}
    \item{st}{a vector of time series of the transition function}
    \item{vecR}{a \eqn{(T \times N^{2})} matrix of Smooth Transition Conditional Correlations}
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
  Silvennoinen, A. and T. Ter\"asvirta (2005),
  \dQuote{Multivariate Autoregressive Conditional Heteroskedasticity with 
          Smooth Transitions in Conditional Correlations.}
  \emph{SSE/EFI Working Paper Series in Economics and Finance} No. 577, 
  Stockholm School of Economics, 
  available at \url{http://swopec.hhs.se/hastef/abs/hastef0577.htm}.
}

\seealso{
    \code{\link{dcc.sim}}, \code{\link{eccc.sim}}
}

\examples{
# Simulating data from the original STCC-GARCH(1,1) process
nobs <- 1000; cut <- 1000
a <- c(0.003, 0.005, 0.001)
A <- diag(c(0.2,0.3,0.15))
B <- diag(c(0.79, 0.6, 0.8))
# Conditional Correlation Matrix for regime 1
R1 <- matrix(c(1.0, 0.4, 0.3, 0.4, 1.0, 0.12, 0.3, 0.12, 1.0),3,3) 
# Conditional Correlation Matrix for regime 2
R2 <- matrix(c(1.0, 0.01, -0.3, 0.01, 1.0, 0.8, -0.3, 0.8, 1.0),3,3)
# a parameter vector for the scale and location parameters 
# in the logistic function
tr.para <- c(5,0)               
# a parameter vector for a GARCH(1,1) transition variable
st.para <- c(0.02,0.04, 0.95)   
nu <- 15

\dontrun{
stcc.data <- stcc.sim(nobs, a, A, B, R1, R2, 
                      tr.par=tr.para, st.par=st.para, model="diagonal")
stcc.data.t. <- stcc.sim(nobs, a, A, B, R1, R2, 
                      tr.par=tr.para, st.par=st.para, d.f=nu, model="diagonal")
}
}

\keyword{ts}
\keyword{datagen}
\keyword{models}
\keyword{multivariate}
