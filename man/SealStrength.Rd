\name{SealStrength}
\alias{SealStrength}
\docType{data}
\title{
Data for Seal Strength
}
\description{
Seal Strength dataset presented in Li and Doganaksoy (2014).

In the R package ADDT, Seal Strength dataset has minor modifications where the temperature at time point 0 is changed to 200 degree where in the original dataset, it is 100 degree. This is a computing trick that will not affect the model results.
}
\usage{data(SealStrength)}
\format{
  A data frame with 210 observations on the following 3 variables.
  \describe{
    \item{\code{TempC}}{Temperature.}
    \item{\code{TimeH}}{Time in hours.}
    \item{\code{Response}}{Strength (Newtons)}
  }
}
\details{
Seal strength dataset is collected under an experiment of testing the strength of a new seal.
}
\source{
M. Li, and N. Doganaksoy (2014). Batch variability in accelerated-degradation testing. \emph{Journal of Quality Technology 46, 171-180}.
}

