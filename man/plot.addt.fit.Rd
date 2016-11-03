\name{plot.addt.fit}
\alias{plot.addt.fit}
\title{
ADDT Plotting
}
\description{
Provides graphical tools for ADDT analysis.
}
\usage{
\method{plot}{addt.fit}(x, type, ...)
}
\arguments{
  \item{x}{
An \emph{addt.fit} object.
}
  \item{type}{
Type of plot ("data", "ML", "LS", "SEMI") to be generated. If \emph{type="data"}, a scatter plot of the dataset in \emph{addt.fit} object is generated. If \emph{type="ML"}, a plot of data along with the results of ML approach is produced. If \emph{type="LS"}, the polynominal fit to the data is plotted. If \emph{type="SEMI"}, the semi-parametric fit to the data is plotted.
}
   \item{\dots}{
Optional arguments.
}
}
\seealso{
\code{\link{addt.fit}}
}
