\name{ECFsup-package}
\alias{ECFsup-package}
\alias{ECFsup}
\docType{package}
\title{
\packageTitle{ECFsup}
}
\description{
\packageDescription{ECFsup}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{ECFsup}
\packageIndices{ECFsup}
Testing equality of several covariance functions of functional data. Four different methods are implemented: L2-norm using W-S naive method, L2-norm using W-S bias-reduced method, L2-norm (ZHANG 2013, GUO et al. 2016), and sup-norm using resampling method (GUO et al. 2017).  The input functional data should have been registered and presmoothed. See Ramsay and Silverman (2005) Ch.7 for registration, Zhang (2013) Ch.3 for presmoothing. Tools for preprocessing raw functional data are available
in R package \pkg{fda}, see also Ramsay et al. (2009).
}
\author{
\packageAuthor{ECFsup}

Maintainer: \packageMaintainer{ECFsup}
}
\references{
[GUO et al. (2016)], Guo, J., Zhou, B., Zhang, J.-T. (2016). A further study of an L2-norm based test for the equality of several covariance
functions. \url{https://arxiv.org/abs/1609.04231}

[GUO et al. (2017)], Guo, J., Zhou, B., Zhang, J.-T. (2017). Testing the equality of several covariance functions for functional data:
a supremum-norm based test. \url{https://arxiv.org/abs/1609.04232}

[PAPARODITIS and SAPATINAS (2016)], Paparoditis, E., Sapatinas, T. (2016). Bootstrap-based testing of equality of mean functions or equality of covariance operators for functional data. \emph{Biometrika}, \bold{103}, 727--733. \doi{10.1093/biomet/asw033}

[RAMSAY and SILVERMAN (2005)], Ramsay, J. O., Silverman, B. W. (2005). \emph{Functional Data Analysis}. Springer. \doi{10.1007/b98888}

[RAMSAY et al. (2009)], Ramsay, J. O., Hooker, G., Graves, S. (2009). \emph{Functional data analysis with R and MATLAB}. Springer. \doi{10.1007/978-0-387-98185-7}

[ZHANG (2013)], Zhang, J.-T. (2013). \emph{Analysis of variance for functional data}. CRC Press.
}
\keyword{ package }
\keyword{ htest }
\examples{
fdata<-list();
fdata[[1]]<-matrix(rnorm(200),20,10);
fdata[[1]]<-matrix(rnorm(300),20,15);
KSCovL2(fdata);KSCovsup(fdata);
}
