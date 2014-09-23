\name{summary.addt.fit}
\alias{summary.addt.fit}
\title{
Summary of an "addt.fit" object.
}
\description{
Provides a brief summary of an \emph{addt.fit} object. 
}
\usage{
\method{summary}{addt.fit}(object, ...)
}
\arguments{
  \item{object}{
An \code{addt.fit} object.
}
   \item{\dots}{
Optional arguments.
}
}
\value{
Returns a list whose items vary depending on which procedure was used in \link{addt.fit}. If \emph{proc="ML"}, then the output has the following items:
\item{coef.mle.mat}{A matrix giving the estimated coefficients, their standard errors, and the lower and upper bounds of a 95\% confidence interval.}
\item{ti.CI}{The estimated thermal index of the material, its standard error, and 95\% confidence interval.}
\item{logLik}{The maximum log-likelihood value achieved.}
If \emph{proc="LS"}, the output has the following items:
\item{coefs}{A vector giving the estimated coefficients for the time-temperature line.}
\item{TI}{The estimated thermal index of the material.}
\item{interp.time}{A matrix contains the interpolated time at the \emph{failure.threshold} specified for each temperature level.}
If \emph{proc="Both"}, the output will give the preceeding values for both the \emph{LS.obj} and \emph{ML.obj}.
}
\seealso{
\code{\link{addt.fit}}
}
