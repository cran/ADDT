\name{addt.predint.ybar.mle}
\alias{addt.predint.ybar.mle}
\title{
Prediction of the mean of future observations
}
\description{
Given observations for one temperature level up to some time point, computes a prediction interval for the mean degradation level at some future time point.
}
\usage{
addt.predint.ybar.mle(obj, conflevel, num.fut.obs = 5, temp, tt)
}
\arguments{
  \item{obj}{
An \code{addt.fit} object.
}
  \item{conflevel}{
The confidence level of the prediction interval. This argument is in decimal form (i.e., 95\% is 0.95).
}
  \item{num.fut.obs}{
The number of future observations within a batch at the future time point.
}
  \item{temp}{
The temperature level at which predictions are to be made.
}
  \item{tt}{
The future time point where prediction is desired.
}
}
\value{
Returns a vector containing the lower and upper bounds of the prediction interval.
}
\note{
This function only works with an object resulted from \emph{ML} approach.
}
\references{
Hong, Y., King, C. B., Xie, Y., Van Mullekom, J. H., Dehart, S. P. and DeFeo, P. A. (2014). A Comparison of Least Squares and Maximum Likelihood Approaches to Estimating Thermal Indices for Polymeric Materials. \emph{Technical Report}.
}
\seealso{
\code{\link{addt.fit}}
}
