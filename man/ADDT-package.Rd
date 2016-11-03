\name{ADDT-package}
\alias{ADDT-package}
\alias{ADDT}
\docType{package}
\title{
Accelerated Destructive Degradation Testing
}
\description{
Accelerated destructive degradation tests (ADDT) are often used to collect necessary data for assessing the long-term properties of polymeric materials. Based on the collected data, a thermal index (TI) is estimated. The TI can be useful for material rating and comparisons. This package implements the traditional method based on the least-squares method, the parametric method based on maximum likelihood estimation, and the semiparametric method based on spline methods, and their corresponding methods for estimating TI for polymeric materials. The traditional approach is a two-step approach that is currently used in industrial standards, while the parametric method is widely used in the statistical literature. The semiparametric method is newly developed. The parametric and semiparametric approaches allow one to do statistical inference such as quantifying uncertainties in estimation, hypothesis testing, and predictions. Publicly available datasets are provided for illustrations. More details can be found in Jin et al. (2017).
}
\details{
\tabular{ll}{
Package: \tab ADDT\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2016-10-08\cr
License: \tab GPL-2\cr
}
}
\author{
Yili Hong, Yimeng Xie, Zhongnan Jin, and Caleb King

Maintainer: Yili Hong <yilihong@vt.edu>
}
\references{
C. B. King, Y. Xie, Y. Hong, J. H. Van Mullekom, S. P. DeHart, and P. A. DeFeo, ``A comparison of
traditional and maximum likelihood approaches to estimating thermal indices for polymeric
materials,'' Journal of Quality Technology, in press, 2016.

L. A. Escobar, W. Q. Meeker, D. L. Kugler, and L. L. Kramer, ``Accelerated destructive degradation
tests: Data, models, and analysis,'' in Mathematical and Statistical Methods in Reliability,
B. H. Lindqvist and K. A. Doksum, Eds. River Edge, NJ: World Scientific Publishing
Company, 2003, ch. 21.

M. Li and N. Doganaksoy, ``Batch variability in accelerated-degradation testing,'' Journal of
Quality Technology, vol. 46, pp. 171-180, 2014.

Y. Xie, C. B. King, Y. Hong, and Q. Yang, ``Semi-parametric models for accelerated destructive
degradation test data analysis,'' Preprint: arXiv:1512.03036, 2015.

Y. Xie, Z. Jin, Y. Hong, and J. H. Van Mullekom, ``Statistical methods for thermal index estimation based on accelerated destructive degradation test data,'' in Statistical Modeling for Degradation Data, D. G. Chen, Y. L. Lio, H. K. T. Ng, and T. R. Tsai, Eds. NY: New York: Springer, 2017, ch. 12.

Z. Jin, Y. Xie, Y. Hong, and J. H. Van Mullekom, ``ADDT: An R package for analysis of accelerated destructive degradation test data,'' in Statistical Modeling for Degradation Data, D. G. Chen, Y. L. Lio, H. K. T. Ng, and T. R. Tsai, Eds. NY: New York: Springer, 2017, ch. 14.
}
\keyword{package}
