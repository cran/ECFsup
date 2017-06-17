// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double rcpparma_L2stat(const arma::mat & V, const arma::uvec & sample_num, int k, int p, int n)
{
  arma::mat Sigma;
  if (n > p)
  {
    Sigma = (V * V.t()) / double(n - k);
  }else
  {
    Sigma = (V.t() * V) / double(n - k);
  }
  
  double stat = 0.0;
  int ind = 0;
  for (int i=0; i<k; i++)
  {
    int ni = sample_num(i);
    arma::mat Vi = V.cols(ind, ind+ni-1);
    arma::mat Si;
    double temp;
    if (n > p)
    {
      Si = Vi * Vi.t() / double(ni - 1);
      temp = arma::trace((Si - Sigma) * (Si - Sigma));
    }else{
      Si = Vi.t() * Vi / double(ni - 1);
      temp = arma::trace(Si * Si) - 2 * arma::trace(Vi.t() * V * V.t() * Vi) / double(n - k) / double(ni - 1) 
        + arma::trace(Sigma * Sigma);
    }
    stat = stat + (ni - 1) * temp;
    ind = ind + ni;
  }
  return(stat);
}

// [[Rcpp::export]]
double rcpparma_fKSCovL2(const arma::mat & data, const arma::uvec & sample_num, int k, int p, int n, int Nsim = 1000) {
  
  arma::mat V(p,n);//centered data
  int ind=0;
  for (int i=0; i<k; i++)
  {
    arma::mat datai = data.cols(ind,ind + sample_num(i)-1);
    arma::colvec mui = arma::mean(datai, 1);
    arma::mat one1ni(1,sample_num(i)); one1ni.ones();
    V.cols(ind,ind + sample_num(i)-1) = datai - mui * one1ni;
    ind = ind + sample_num(i);
  }
  
  double stat = rcpparma_L2stat(V, sample_num, k, p, n);
  
  arma::colvec stats(Nsim);//stats.fill(0.0);
  
  //int cnt=0;
  for(int i=0; i<Nsim; i++)
  {
    arma::mat pdata = arma::shuffle(V, 1);
    stats(i) = rcpparma_L2stat(pdata, sample_num, k, p, n);
    //if(stats(i)>stat)cnt++;
  }
  arma::uvec result=(stats>stat);
  double pvalue = arma::sum(result)/double(Nsim);//mean gives 0!
  return pvalue;
}

// [[Rcpp::export]]
double rcpparma_maxstat(const arma::mat & V, const arma::uvec & sample_num, int k, int p, int n)
{
  arma::mat Sigma = (V * V.t()) / double(n - k);
  
  arma::mat stat(p,p);stat.zeros();
  int ind = 0;
  for (int i=0; i<k; i++)
  {
    int ni = sample_num(i);
    arma::mat Vi = V.cols(ind, ind+ni-1);
    arma::mat Si = Vi * Vi.t() / double(ni - 1);
    arma::mat temp = (Si - Sigma) % (Si - Sigma);
    stat = stat + (ni - 1) * temp;
    ind = ind + ni;
  }
  return(stat.max());
}
// [[Rcpp::export]]
double rcpparma_fKSCovsup(const arma::mat & data, const arma::uvec & sample_num, int k, int p, int n, int Nsim = 1000) {
  
  arma::mat V(p,n);//centered data
  int ind=0;
  for (int i=0; i<k; i++)
  {
    arma::mat datai = data.cols(ind,ind + sample_num(i)-1);
    arma::colvec mui = arma::mean(datai, 1);
    arma::mat one1ni(1,sample_num(i)); one1ni.ones();
    V.cols(ind,ind + sample_num(i)-1) = datai - mui * one1ni;
    ind = ind + sample_num(i);
  }
  
  double stat = rcpparma_maxstat(V, sample_num, k, p, n);
  
  arma::colvec stats(Nsim);//stats.fill(0.0);
  
  //int cnt=0;
  for(int i=0; i<Nsim; i++)
  {
    arma::mat pdata = arma::shuffle(V, 1);
    stats(i) = rcpparma_maxstat(pdata, sample_num, k, p, n);
    //if(stats(i)>stat)cnt++;
  }
  arma::uvec result=(stats>stat);
  double pvalue = arma::sum(result)/double(Nsim);//mean gives 0!
  return pvalue;
}
