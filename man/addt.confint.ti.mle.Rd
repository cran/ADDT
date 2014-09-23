\name{addt.confint.ti.mle}
\alias{addt.confint.ti.mle}
\title{
Confidence Interval for Thermal Index (TI) Estimate.
}
\description{
Computes a confidence interval for the TI estimate.
}
\usage{
addt.confint.ti.mle(obj, conflevel)
}
\arguments{
  \item{obj}{
An \code{addt.fit} object.
}
  \item{conflevel}{
Confidence level in decimal form (i.e., 95\% is 0.95)
}
}
\note{
This function only works with an object resulted from the \emph{ML} approach.
}
\value{
Returns a vector containing the estimated TI, standard error, and lower and upper confidence limits.
}
\seealso{
\code{\link{addt.fit}}
}
