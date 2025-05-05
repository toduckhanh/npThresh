#define ARMA_NO_DEBUG
#define ARMA_USE_BLAS

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;
using namespace arma;


// define the k(u) kernel function
static inline double kernel(int Ktype, double x){
  if(Ktype == 0) return exp(-x*x/2)/sqrt(2*M_PI); // Gauss kernel
  else{ // if(Ktype == 1)
    if(fabs(x) <= 1) return 3*(1 - x*x)/4; // Epanechnikov kernel
    else return 0;
  }
  // else return 0;
}

// define the K(u) cumulative kernel function
static inline double Kernel(int Ktype, double x){
  if(Ktype == 0) return R::pnorm(x, 0.0, 1.0, 1, 0); // Gauss kernel
  else{ // if(Ktype == 1)
    if(fabs(x) <= 1) return 0.5 + 3*(x - x*x*x/3)/4; // Epanechnikov kernel
    else if(x > 1) return 1;
    else return 0;
  }
  // else return 0;
}

// define the F_kernel
// [[Rcpp::export]]
NumericVector cdf_kernel_C(NumericVector x, NumericVector X, int Ktype,
                         double bwd){
  // x: points to compute CDF
  // X: data
  // Ktype: type of kernel, 0 - Gaussian, 1 - Epanechnikov
  // bwd: bandwidth
  int n = X.size(), m = x.size();
  double uij = 0;
  NumericVector num(m);
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      uij = (x[i] - X[j])/bwd;
      num[i] += Kernel(Ktype, uij);
    }
  }
  return num/n;
}

