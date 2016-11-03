\name{addt.fit}
\alias{addt.fit}
\title{
ADDT Model Fitting
}
\description{
Fits degradation data using the least-squares, maximum likelihood, and semiparametric methods and estimates the thermal indices (TI).
}
\usage{
addt.fit(formula, data, initial.val = 100, proc = "All",
    failure.threshold, time.rti = 1e+05, method = "Nelder-Mead",
    subset, na.action, starts = NULL, fail.thres.vec = c(70, 80),
    semi.control = list(cor = F,...),...)
}
\arguments{
  \item{formula}{
A formula of type \emph{Response~Time+Temperature}. The order of \emph{Time} and \emph{Temperature} cannot be changed.
}
  \item{data}{
A data frame contains the ADDT data to be analyzed.
}
  \item{initial.val}{
We need response measurements at time point 0 to compute the initial degradation level in the model. If the data does not contain that information, user must supply the initial.value. Otherwise, the function will give an error message.
}
  \item{proc}{
The type of analysis to be performed which can be "LS" (least squares), "ML" (maximum likelihood), "SemiPara" (semiparametric), or "All" (the least-squares, maximum likelihood and semiparametric methods).
}
  \item{failure.threshold}{
The value below which a soft failure occurs. Must be in the form of a percent of the initial value.
}
  \item{time.rti}{
The time associated with the thermal index (TI) or relative thermal index (RTI). Typically 100,000 hours.
}
  \item{method}{
An argument passed to \link{optim} specifying the optimization procedure. Default is \code{Nelder-Mead}.
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
If the user does not specify \emph{``starts''} argument, the user may specify a vector of two different \emph{failure.thresholds}. The least-squares procedure is then used for the two different \emph{failure.thresholds} to produce starting values for the maximum likelihood procedure.
}
  \item{semi.control}{\emph{list=(cor=F,...)}, control parameters related to the semiparametric method. If not specified, default is no correlation. When \emph{cor = T} in the list is specified, model assumes a correlation term.

}

  \item{\dots}{Optional arguments.}
}

\value{
An object of class "addt.fit", which is a list containing:
\item{LS.obj}{If least-squares approach is used, a LS.obj will returned. It contains estimates of coefficients in the TI formula, estimated TI, a matrix contains the failure time by polynomial interpolation.}
 \item{ML.obj}{A ML.obj object is returned if maximum likelihood approach is specified.}
 \item{SemiPara.obj}{A SemiPara.obj is returned if SemiPara approach is specified.}
 \item{dat}{The data set used in least square/maximum likelihood approaches.}
 \item{time.rti}{An argument stored to be used for functions related to "addt.fit" object. }
 \item{initial.val}{An argument stored to be used for functions related to "addt.fit" object.}
 \item{failure.threshold}{An argument stored to be used for functions related to "addt.fit" object.}
}
\details{
A thermal index (TI) or relative thermal index (RTI) is often used to evaluate long-term performance of polymeric materials. Accelerated destructive degradation testing (ADDT) is widely used to calculate the TI of certain polymeric materials.

The dataset considered in \link{addt.fit} function contain repeated measurements of a response, say tensile strength, at some combinations of time and temperature.

The least squares procedure aggregates data into the average of measurements at each combination of time and temperature. Then, polynominal regression is used to interpolate the failure time for each combination. A least squares line is fitted to the failure time data and the TI is then obtained by
\deqn{TI=\frac{beta1}{log10(time.rti)-beta0}-273.16.} It is important to note that observations are required after failure in order for this procedure to be successful.

The maximum likelihood procedure assumes a degradation path dependent on time and temperature. An example of a parametric form for this path can be found in Vaca-Trigo and Meeker (2009) and is the form currently used here. The error term is assumed to follow a multivariate normal distribution. A TI can be directly estimated from the parameter estimates for the degradation path. The \code{addt.fit} function will be generalized to allow other parametric forms of the mean function, and/or other non-Guassian distribution in later versions.

The semiparametric model assembles both parametric model like Arrhenius model for degradation variable extrapolation as well as non-parametric model in order to be more compatible for various materials.
}
\references{
Y. Hong, C. B. King, Y. Xie, J. H. Van Mullekom, S. P. Dehart, and P. A. DeFeo (2014). ``A Comparison of Least Squares and Maximum Likelihood Approaches to Estimating Thermal Indices for Polymeric Materials,''  Journal of Quality Technology, in press, 2016.

I. Vaca-Trigoand, W. Q. Meeker, ``A statistical model for linking field and laboratory exposure results for a model coating,'' in Service Life Prediction of Polymeric Materials, J. Martin, R. A. Ryntz, J. Chin, and R. A. Dickie, Eds. NY: New York: Springer, 2009, ch. 2.

Y. Xie, C. B. King, Y. Hong, and Q. Yang, ``Semi-parametric models for accelerated destructive
degradation test data analysis,'' Preprint: arXiv:1512.03036, 2015.

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

## Semiparametric
\dontrun{addt.fit.semi<-addt.fit(Response~TimeH+TempC,data=AdhesiveBondB,proc="SemiPara",
failure.threshold=70)}

## All LS, ML and Semi-Parametric procedures
\dontrun{addt.fit.all<-addt.fit(Response~TimeH+TempC,data=AdhesiveBondB,proc="All",
failure.threshold=70)}

summary(addt.fit.lsa)
summary(addt.fit.mla)
\dontrun{summary(addt.fit.semi)}
\dontrun{summary(addt.fit.all)}

\dontrun{plot(addt.fit.all, type="data")}
\dontrun{plot(addt.fit.all, type="LS")}
\dontrun{plot(addt.fit.all, type="ML")}
\dontrun{plot(addt.fit.semi, type="SEMI")}

\dontrun{addt.confint.ti.mle(addt.fit.mla,conflevel=0.95)}
}
