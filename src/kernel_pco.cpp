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
  if(Ktype == 0) return exp(-x*x/2.0)/sqrt(2.0*M_PI); // Gauss kernel
  else{ // if(Ktype == 1)
    if(fabs(x) <= 1) return 3.0*(1.0 - x*x)/4.0; // Epanechnikov kernel
    else return 0;
  }
  // else return 0;
}

// PCO
// [[Rcpp::export]]
List PCO_bwd_C(NumericVector X, int Ktype, NumericVector bwd_seq){
  arma::vec bwds = as<vec>(bwd_seq);
  int n = X.size(), bwd_pts = bwds.size();
  double h_min = bwds.min();
  double h_min2 = pow(h_min, 2);
  arma::vec hw(bwd_pts, fill::zeros);
  arma::vec hw2(bwd_pts, fill::zeros);
  arma::vec pen(bwd_pts, fill::zeros);
  hw = sqrt(pow(bwds, 2) + h_min2);
  hw2 = sqrt(2.0)*bwds;
  pen = sqrt(2.0)/(n*sqrt(M_PI)*hw);
  arma::vec f_target(bwd_pts, fill::zeros);
  double ujl = 0.0, ujl_1 = 0.0, ujl_2 = 0.0;
  for(int i = 0; i < bwd_pts; i++){
    for(int j = 0; j < n; j++){
      for(int l = 0; l < n; l++){
        ujl = X[j] - X[l];
        ujl_1 = ujl/hw2[i];
        ujl_2 = ujl/hw[i];
        f_target[i] += kernel(Ktype, ujl_1)/hw2[i] -
          2.0*kernel(Ktype, ujl_2)/hw[i];
      }
    }
  }
  f_target = f_target/pow(n, 2) + pen;
  double id_h_pco = f_target.index_min();
  List out;
  out["bwd_seq"] = wrap(bwd_seq);
  out["f_target"] = wrap(f_target);
  out["h_opt"] = bwd_seq[id_h_pco];
  return out;
}
