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
