\name{addt.fit}
\alias{addt.fit}
\title{
ADDT Model Fitting
}
\description{
Performs the least squares and maximum likelihood procedures for calculating thermal indices (TI) for polymeric materials.
}
\usage{
addt.fit(formula, data, initial.val = 100, proc = "Both",     
    failure.threshold, time.rti = 1e+05, method = "Nelder-Mead",
    subset, na.action, starts = NULL, fail.thres.vec = c(70, 80),            
    ...)
}
\arguments{
  \item{formula}{
A formula of type \emph{Response~Time+Temperature}. The order of \emph{Time} and \emph{Temperature} cannot be changed.
}
  \item{data}{
A data frame containing the ADDT data to be analyzed.
}
  \item{initial.val}{
If the data don't contain measurements at time 0, the user needs to specify the initial value for the response. Otherwise, the measurements at time 0 will be used. Default is 100.
}
  \item{proc}{
The type of analysis to be performed. Can be "LS" (least squares), "ML" (maximum likelihood), or "Both" (least squares and maximum likelihood).
}
  \item{failure.threshold}{
The value below which a soft failure occurs. Currently must be in the form of a percent of the initial value. See the example for usage.
}
  \item{time.rti}{
The time associated with the thermal index (TI) or relative thermal index (RTI). Typically 100,000 hours.
}
  \item{method}{
An argument passed to \link{optim} specifying the optimization procedure. Default is \code{"Nelder-Mead"}.
}
  \item{subset}{
An optional statement that allows the use of only part of the dataset.
}
  \item{na.action}{
Indicates the action required when data contain missing values.
}
  \item{starts}{
A vector of starting values for the maximum likelihood procedure. 
See \code{fail.thres.vec} for alternative.
}
  \item{fail.thres.vec}{
If the user does not specify \emph{"starts"} argument, the user may instead specify a vector of two different \emph{failure.thresholds}. The least square procedure is then used for the two different \emph{failure.thresholds} to produce starting values for the maximum likelihood procedure.
}
  \item{\dots}{
Optional arguments.
}
}
\value{
An object of class "addt.fit", which is a list containing: 
\item{LS.obj}{If least square approach is used, a LS.obj will returned. It contains estimates of coefficients in the TI formula, estimated TI, a matrix contains the failure time by polynomial interpolation.}
 \item{ML.obj}{A ML.obj object is returned if maximum likelihood approach is specified.}
 \item{dat}{The data set used in least square/maximum likelihood approaches.}
 \item{time.rti}{An argument stored to be used for functions related to "addt.fit" object. }
 \item{initial.val}{An argument stored to be used for functions related to "addt.fit" object.}
 \item{failure.threshold}{An argument stored to be used for functions related to "addt.fit" object.}
}
\details{
A thermal index (TI) or relative thermal index (RTI) is often used to evaluate long-term performance of polymeric materials. Accelerated destructive degradation testing (ADDT) is widely used to calculate the TI of certain polymeric materials. The least squares and maximum likelihood procedures are the most common procedures used to estimate TI. 

The dataset considered in \link{addt.fit} function contain repeated measurements of a response, say tensile strength, at some combinations of time and temperature. 

The least squares procedure aggregates data into the average of measurements at each combination of time and temperature. Then, polynominal regression is used to interpolate the failure time for each combination. A least squares line is fitted to the failure time data and the TI is then obtained by 
\deqn{TI=\frac{beta1}{log10(time.rti)-beta0}-273.16.} It is important to note that observations are required after failure in order for this procedure to be successful. 

The maximum likelihood procedure assumes a degradation path dependent on time and temperature. An example of a parametric form for this path can be found in Vaca-Trigo and Meeker (2009) and is the form currently used here. The error term is assumed to follow a multivariate normal distribution. A TI can be directly estimated from the parameter estimates for the degradation path. The \code{addt.fit} function will be generalized to allow other parametric forms of the mean function, and/or other non-Guassian distribution in later versions.
}
\references{
Hong, Y., King, C. B., Xie, Y., Van Mullekom, J. H., Dehart, S. P. and DeFeo, P. A. (2014). A Comparison of Least Squares and Maximum Likelihood Approaches to Estimating Thermal Indices for Polymeric Materials. \emph{Technical Report}.

Vaca-Trigo, I. and W. Q. Meeker (2009). A statistical model for linking field and laboratory exposure results for a model coating. \emph{NY: New York: Springer.}
}
\seealso{
\link{plot.addt.fit}, \link{summary.addt.fit}
}
\examples{
data(AdhesiveBondB)

## Least Squares
addt.fit.lsa<-addt.fit(Response~TimeH+TempC,data=AdhesiveBondB,proc="LS",
failure.threshold=70)

## Maximum Likelihood
addt.fit.mla<-addt.fit(Response~TimeH+TempC,data=AdhesiveBondB,proc="ML",
failure.threshold=70)

## Both LS and ML procedures
addt.fit.both<-addt.fit(Response~TimeH+TempC,data=AdhesiveBondB,proc="Both",
failure.threshold=70)

summary(addt.fit.lsa)
summary(addt.fit.mla)
summary(addt.fit.both)

plot(addt.fit.both, type="data")
plot(addt.fit.both, type="LS")
plot(addt.fit.both, type="ML")

addt.confint.ti.mle(addt.fit.mla,conflevel=0.95)
}
