\name{loglik.dcc}
\alias{loglik.dcc}
\title{ The log-likelihood function for the (E)DCC GARCH model  }

\description{
This function returns a log-likelihood of the (E)DCC-GARCH model.
}

\usage{
   loglik.dcc(param, dvar, model)
}

\arguments{
  \item{param}{a vector of all the parameters in the (E)DCC-GARCH model. }
  \item{dvar}{a matrix of the observed residuals \eqn{(T \times N)}}
  \item{model}{a character string describing the model. "diagonal" for the diagonal model 
               and "extended" for the extended (full ARCH and GARCH parameter matrices) model}
}

\value{
    the negative of the full log-likelihood of the (E)DCC-GARCH model
}

\note{
    this function is currently not in use.
}

\references{
  Engle, R.F. and K. Sheppard (2001),
  \dQuote{Theoretical and Empirical Properties of Dynamic
  Conditional Correlation Multivariate GARCH.}
  \emph{Stern Finance Working Paper Series}
  {FIN}-01-027 (Revised in Dec. 2001),
  New York University Stern School of Business.
  
  Engle, R.F. (2002),
  \dQuote{Dynamic Conditional Correlation: A Simple Class of
  Multivariate Generalised Autoregressive Conditional
  Heteroskedasticity Models.}
  \emph{Journal of Business and Economic Statistics}
  \bold{20}, 339--350.
}

\keyword{ts},
\keyword{models},
\keyword{multivariate}

