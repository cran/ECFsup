#' @title L2-norm test using W-S approximation for equality of several covariance functions
#' @export
#' 
#' @description
#' L2-norm test of equality of several covariance functions, using the 
#' naive or bias-reduced method (Welch--Satterthwaite approximation)  
#' to approximate the null distribution.
#' 
#' @details
#' L2-norm test of equality of several covariance functions. The null distribution will be
#' approximated by a scaled chi-squared random variable. Two approximation methods are implemented:
#' naive method and bias-reduced method, which work for Gaussian data only.
#' The bias-reduced method is more accurate than the naive method for Gaussian data. 
#' The input functional data should have been registered and presmoothed.
#' See Ramsay and Silverman (2005) Ch.7 for registration, and Zhang (2013) Ch.3 for presmoothing.
#' Tools for preprocessing raw functional data are available in R package \pkg{fda}, 
#' see also Ramsay et al. (2009).
#' 
#' @param data The list variable containing k groups of presmoothed functional observations. 
#' Each element of the list is a p (number sampling points) by n (sample size) matrix.
#' @param apprflag Approximation method, 0: naive method, 1: bias-reduced method.
#' @param method placeholder for L2 resampling method, for testing purpose, should not use currently.
#' @param Nsim placeholder for L2 resampling method, for testing purpose, should not use currently.
#' @return The p-value of the test.
#' @examples
#' fdata<-list();
#' fdata[[1]]<-matrix(rnorm(200),20,10);
#' fdata[[1]]<-matrix(rnorm(300),20,15);
#' KSCovL2WS(fdata, 0)
#' KSCovL2WS(fdata, 1)
#' @references ZHANG (2013), GUO et al. (2016), RAMSAY and SILVERMAN (2005), RAMSAY et al. (2009).
#' @seealso \code{\link{KSCovL2}}, \code{\link{KSCovsup}}.
KSCovL2WS <- function(data,apprflag = 0,method = 1,Nsim = 1000)
{
  #require("foreach");
  dim <- dim(data[[1]]);
  p <- dim[1];
  k <- length(data);
  sample_num <- sapply(data, ncol);
  n <- sum(sample_num);
  V <- matrix(,p,n);
  ind <- 0;
  for (i in 1:k)
  {
    mui <- rowMeans(data[[i]]);
    V[,(ind + 1):(ind + sample_num[i])] <-
      data[[i]] - mui %*% t(rep(1,sample_num[i]));
    ind <- ind + sample_num[i];
  }
  if (n > p)
  {
    Sigma <- (V %*% t(V)) / (n - k);
  }else
  {
    Sigma <- (t(V) %*% V) / (n - k);
  }
  
  stat <- 0;
  ind <- 0;
  for (i in 1:k)
  {
    ni <- sample_num[i];
    Vi <- V[,(ind + 1):(ind + ni)];
    if (n > p)
    {
      Si <- Vi %*% t(Vi) / (ni - 1);
      temp <- sum(diag((Si - Sigma) %*% (Si - Sigma)));
    }else{
      Si <- t(Vi) %*% Vi / (ni - 1);
      temp <-
        sum(diag(Si %*% Si)) - 2 * sum(diag(t(Vi) %*% V %*% t(V) %*% Vi)) / (n - k) / (ni - 1) + sum(diag(Sigma %*% Sigma));
    }
    stat <- stat + (ni - 1) * temp;
    ind <- ind + ni;
  }
  
  if (method == 1)
    # L2-norm based test
  {
    A <- sum(diag(Sigma));B <- sum(diag(Sigma %*% Sigma));
    if (apprflag == 0)
      # the naive method
    {
      A2 <- A %*% A;B2 <- B;
    }
    else if (apprflag == 1)
      #the bias-reduced method
    {
      A2 <- (n - k) * (n - k + 1) / (n - k - 1) / (n - k + 2) * (A %*% A - 2 * B / (n - k + 1));
      B2 <- (n - k) ^ 2 / (n - k - 1) / (n - k + 2) * (B - A %*% A / (n - k));
    }
    A3 <- B2 + A2;
    B3 <- 2 * sum(diag(Sigma %*% Sigma %*% Sigma %*% Sigma)) + 2 * B2 %*% B2;
    
    hk <- (k - 1) * A3 %*% A3 / B3;hb <- B3 / A3;
    pvalue <- pchisq(
      q = stat / hb, df = hk, ncp = 0, lower.tail = FALSE, log.p = FALSE
    );
    
    pstat <- c(stat,pvalue); #print(stat);
  }
  else
  {
    stats <- foreach(i=1:Nsim, .combine = rbind) %dopar% oneSPL(V, sample_num, k, p, n);
    pvalue <- mean(stats>stat);
  }
  return(pvalue);
}
#' @title oneSPL of L2-norm test
#' 
#' @description
#' Generate one pseudo sample by resampling and compute the value of the L2-norm test statistic.
#' 
#' @details
#' The input data should have been centered. This function is obsolete and implemented in R, 
#' for testing purpose only. Should use \code{\link{oneSPLL2}} instead.
#' 
#' @param odata The list variable containing k groups of presmoothed functional observations. 
#' Each element of the list is a p (number sampling points) by n (sample size) matrix.
#' @param sample_num vector of group sizes.
#' @param k number of groups.
#' @param p number of time points.
#' @param n total number of samples.
#' @return The value of test statistic.
#' @examples
#' p <- 100; sample_num <- c(40,60); k <- length(sample_num); n <- sum(sample_num);
#' odata <- matrix(rnorm(p*n),p,n);
#' oneSPL(odata, sample_num, k, p, n);
#' @references ZHANG (2013), GUO et al. (2016).
#' @seealso \code{\link{oneSPLL2}}, \code{\link{oneSPLmax}}.
oneSPL <- function(odata, sample_num, k, p, n)
{
  flag <- sample.int(n, size = n);
  V <- odata[,flag];

  if (n > p)
  {
    Sigma <- (V %*% t(V)) / (n - k);
  }else
  {
    Sigma <- (t(V) %*% V) / (n - k);
  }
  
  stat <- 0;
  ind <- 0;
  for (i in 1:k)
  {
    ni <- sample_num[i];
    Vi <- V[,(ind + 1):(ind + ni)];
    if (n > p)
    {
      Si <- Vi %*% t(Vi) / (ni - 1);
      temp <- sum(diag((Si - Sigma) %*% (Si - Sigma)));
    }else{
      Si <- t(Vi) %*% Vi / (ni - 1);
      temp <-
        sum(diag(Si %*% Si)) - 2 * sum(diag(t(Vi) %*% V %*% t(V) %*% Vi)) / (n - k) / (ni - 1) + sum(diag(Sigma %*% Sigma));
    }
    stat <- stat + (ni - 1) * temp;
    ind <- ind + ni;
  }
  return(stat);
}

#' @title L2-norm based test of equality of several covariance functions
#' @export
#' 
#' @description
#' L2-norm test of equality of several covariance functions, using resampling  
#' to approximate the null distribution.
#' 
#' @details
#' L2-norm test of equality of several covariance functions, see Zhang (2013), Guo et al. (2016).
#' 
#' @param data The list variable containing k groups of presmoothed functional observations. 
#' Each element of the list is a p (number sampling points) by n (sample size) matrix.
#' @param Nsim Number of pseudo samples by resampling, default=1000.
#' @return The p-value of the test.
#' @examples
#' fdata<-list();
#' fdata[[1]]<-matrix(rnorm(200),20,10);
#' fdata[[1]]<-matrix(rnorm(300),20,15);
#' KSCovL2(fdata)
#' KSCovL2(fdata, 500)
#' @references ZHANG (2013), GUO et al. (2016), PAPARODITIS and SAPATINAS (2016).
#' @seealso \code{\link{KSCovL2WS}}, \code{\link{KSCovsup}}.
KSCovL2 <- function(data,Nsim = 1000)
{
  #require("foreach");
  dim <- dim(data[[1]]);
  p <- dim[1];
  k <- length(data);
  sample_num <- sapply(data, ncol);
  n <- sum(sample_num);
  V <- matrix(,p,n);
  ind <- 0;
  for (i in 1:k)
  {
    mui <- rowMeans(data[[i]]);
    V[,(ind + 1):(ind + sample_num[i])] <-
      data[[i]] - mui %*% t(rep(1,sample_num[i]));
    ind <- ind + sample_num[i];
  }
  stat <- rcpparma_L2stat(V, sample_num, k, p, n);
  stats <- foreach(i=1:Nsim, .combine = rbind) %dopar% oneSPLL2(V, sample_num, k, p, n);
  pvalue <- mean(stats>stat);
  return(pvalue);
}
#' @title oneSPL of L2-norm test
#' 
#' @description
#' Generate one pseudo sample by resampling and compute the value of the L2-norm test statistic.
#' 
#' @details
#' The input data should have been centered.
#' 
#' @param odata The list variable containing k groups of presmoothed functional observations. 
#' Each element of the list is a p (number sampling points) by n (sample size) matrix.
#' @param sample_num vector of group sizes.
#' @param k number of groups.
#' @param p number of time points.
#' @param n total number of samples.
#' @return The value of test statistic.
#' @examples
#' p <- 100; sample_num <- c(40,60); k <- length(sample_num); n <- sum(sample_num);
#' odata <- matrix(rnorm(p*n),p,n);
#' oneSPL(odata, sample_num, k, p, n);
#' @references ZHANG (2013), GUO et al. (2016).
#' @seealso \code{\link{oneSPL}}, \code{\link{oneSPLmax}}.
oneSPLL2 <- function(odata, sample_num, k, p, n)
{
  flag <- sample.int(n, size = n);
  V <- odata[,flag];
  return(rcpparma_L2stat(V, sample_num, k, p, n));
}

#' @title Sup-norm based test of equality of several covariance functions
#' @export
#' 
#' @description
#' Sup-norm test of equality of several covariance functions, using resampling 
#' to approximate the null distribution.
#' 
#' @details
#' Sup-norm test of equality of several covariance functions, see GUO et al. (2017).
#' 
#' @param data The list variable containing k groups of presmoothed functional observations. 
#' Each element of the list is a p (number sampling points) by n (sample size) matrix.
#' @param Nsim Number of pseudo samples by resampling, default=1000.
#' @return The p-value of the test.
#' @examples
#' fdata<-list();
#' fdata[[1]]<-matrix(rnorm(200),20,10);
#' fdata[[1]]<-matrix(rnorm(300),20,15);
#' KSCovsup(fdata)
#' KSCovsup(fdata, 500)
#' @references GUO et al. (2017).
#' @seealso \code{\link{KSCovL2WS}}, \code{\link{KSCovL2}}.
KSCovsup <- function(data,Nsim = 1000)
{
  #require("foreach");
  dim <- dim(data[[1]]);
  p <- dim[1];
  k <- length(data);
  sample_num <- sapply(data, ncol);
  n <- sum(sample_num);
  V <- matrix(,p,n);
  ind <- 0;
  for (i in 1:k)
  {
    mui <- rowMeans(data[[i]]);
    V[,(ind + 1):(ind + sample_num[i])] <-
      data[[i]] - mui %*% t(rep(1,sample_num[i]));
    ind <- ind + sample_num[i];
  }
  stat <- rcpparma_maxstat(V, sample_num, k, p, n);
  stats <- foreach(i=1:Nsim, .combine = rbind) %dopar% oneSPLmax(V, sample_num, k, p, n);
  pvalue <- mean(stats>stat);
  return(pvalue);
}
#' @title oneSPL of sup-norm test
#' 
#' @description
#' Generate one pseudo sample by resampling and compute the value of the sup-norm test statistic.
#' 
#' @details
#' The input data should have been centered.
#' 
#' @param odata The list variable containing k groups of presmoothed functional observations. 
#' Each element of the list is a p (number sampling points) by n (sample size) matrix.
#' @param sample_num vector of group sizes.
#' @param k number of groups.
#' @param p number of time points.
#' @param n total number of samples.
#' @return The value of test statistic.
#' @examples
#' p <- 100; sample_num <- c(40,60); k <- length(sample_num); n <- sum(sample_num);
#' odata <- matrix(rnorm(p*n),p,n);
#' oneSPL(odata, sample_num, k, p, n);
#' @references ZHANG (2013), GUO et al. (2016).
#' @seealso \code{\link{oneSPL}}, \code{\link{oneSPLL2}}.
oneSPLmax <- function(odata, sample_num, k, p, n)
{
  flag <- sample.int(n, size = n);
  V <- odata[,flag];
  return(rcpparma_maxstat(V, sample_num, k, p, n));
}