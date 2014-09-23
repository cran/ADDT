\name{addt.mean.summary}
\alias{addt.mean.summary}
\title{
ADDT Batch Means Summary
}
\description{
  Function that returns the average for each time-temperature batch.
}
\usage{
addt.mean.summary(dat)
}
\arguments{
  \item{dat}{
     A dataframe containing the measurements from the ADDT. The dataframe contains temperature, time and response from left to right.
}
}
\value{
  Returns a dataframe giving the mean response for each time-temperature batch.
}
\references{
Hong, Y., King, C. B., Xie, Y., Van Mullekom, J. H., Dehart, S. P. and DeFeo, P. A. (2014). A Comparison of Least Squares and Maximum Likelihood Approaches to Estimating Thermal Indices for Polymeric Materials. \emph{Technical Report}.
}
\examples{
data(AdhesiveBondB)
addt.mean.summary(AdhesiveBondB)
}
