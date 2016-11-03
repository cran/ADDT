\name{addt.mean.summary}
\alias{addt.mean.summary}
\title{
ADDT Batch Means Summary
}
\description{
  Function that returns the averaged responses for each time-temperature batch.
}
\usage{
addt.mean.summary(dat)
}
\arguments{
  \item{dat}{
     A dataframe contains the measurements from the ADDT. The dataframe contains temperature, time, and response from left to right.
}
}
\value{
  Returns a dataframe giving the mean response for each time-temperature batch.
}
\references{
Y. Hong, C. B. King, Y. Xie, J. H. Van Mullekom, S. P. Dehart, and P. A. DeFeo (2014). ``A Comparison of Least Squares and Maximum Likelihood Approaches to Estimating Thermal Indices for Polymeric Materials,''  Journal of Quality Technology, in press, 2016.
}
\examples{
data(AdhesiveBondB)
addt.mean.summary(AdhesiveBondB)
}
