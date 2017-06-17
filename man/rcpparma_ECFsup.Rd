\name{RcppArmadillo-Functions}
\alias{rcpparma_L2stat}
\alias{rcpparma_fKSCovL2}
\alias{rcpparma_maxstat}
\alias{rcpparma_fKSCovsup}
\title{Set of functions in ECPsup package}
\description{
  Tests for equal covariance functions problem, implemented in C++.
}
\usage{
rcpparma_L2stat(V, sample_num, k, p, n);
rcpparma_fKSCovL2(data, sample_num, k, p, n, Nsim);
rcpparma_maxstat(V, sample_num, k, p, n);
rcpparma_fKSCovsup(data, sample_num, k, p, n, Nsim);
}
\arguments{
  \item{V}{centered data matrix}
  \item{data}{data matrix}
  \item{sample_num}{sample sizes vector}
  \item{k}{number of groups}
  \item{p}{number of time points}
  \item{n}{total number of samples}
  \item{Nsim}{number of pseudo samples by resampling}
}
\value{

  \code{rcpparma_L2stat} returns a numeric value computed as the
  test statistic of L2-norm test.

  \code{rcpparma_fKSCovL2} returns a double computed as the p-value of the L2-norm based test.

  \code{rcpparma_maxstat} returns a numeric value computed as the
  test statistic of sup-norm test.
  
  \code{rcpparma_fKSCovsup} returns a double computed as the p-value of the sup-norm based test.
}
\details{
  These are cpp versions of the tests for the ECF problem.
}
\references{
  ZHANG (2013), GUO et al. (2016), PAPARODITIS and SAPATINAS (2016), GUO et al. (2017).
}
\examples{
  vn <- c(20,30,30); k <- length(vn); n <- sum(vn);
  p <- 100; Nsim <- 500;
  datamx <- matrix(rnorm(p*n),p,n,Nsim);
  rcpparma_fKSCovL2(datamx,vn,k,p,n,Nsim);
  rcpparma_fKSCovsup(datamx,vn,k,p,n,Nsim);
}
\author{Bu Zhou}
