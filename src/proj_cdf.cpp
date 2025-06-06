#define ARMA_NO_DEBUG
#define ARMA_USE_BLAS

//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>


using namespace Rcpp;
using namespace arma;

static inline double my_mean(NumericVector x, int n, int type){
  // type: 0 - cos, 1 - sin
  double sum = 0.0;
  if(type == 0){
    for(int i = 0; i < n; i++){
      sum += cos(x[i]);
    }
  } else{
    for(int i = 0; i < n; i++){
      sum += sin(x[i]);
    }
  }
  return sum/n;
}

// define the projection estimator F
// [[Rcpp::export]]
NumericVector cdf_proj_C(NumericVector x, NumericVector X, double a, double b,
                         int M){
  int n = X.size(), m = x.size();
  double rt = b - a;
  NumericVector X_tt = 2.0*M_PI*(X - a)/rt;
  // obtain coefficient B
  NumericVector Bj(2*M + 1);
  Bj[0] = 1/sqrt(b - a);
  NumericVector ui(n);
  int id_even, id_odd;
  for (int i = 1; i < M + 1; i++){
    ui = i*X_tt;
    id_even = 2*i;
    id_odd = id_even - 1;
    Bj[id_even] = my_mean(ui, n, 1);
    Bj[id_odd] = my_mean(ui, n, 0);
  }
  NumericVector out_cdf(m);
  NumericVector cst(m);
  cst = (x - a)/rt;
  double temp1 = 0.0, temp2 = 0.0;
  for(int i = 0; i < m; i++){
    if(x[i] < a) out_cdf[i] = 0.0;
    else if(x[i] > b) out_cdf[i] = 1.0;
    else{
      temp1 = 2*M_PI*(x[i] - a)/rt;
      for(int j = 1; j < M + 1; j++){
        temp2 = temp1*j;
        id_even = 2*j;
        id_odd = id_even - 1;
        out_cdf[i] += (Bj[id_even]*(1 - cos(temp2)) + Bj[id_odd]*sin(temp2))/j;
      }
      out_cdf[i] = out_cdf[i]/M_PI + cst[i];
    }
  }
  return out_cdf;
}

static inline double B2_coef(NumericVector u, double rt, int M, int n){
  int k = 2*M + 1;
  NumericVector Bj(k);
  Bj[0] = 1/sqrt(2);
  NumericVector ui(n);
  int id_even, id_odd;
  for (int i = 1; i < M + 1; i++){
    ui = i*u;
    id_even = 2*i;
    id_odd = id_even - 1;
    Bj[id_even] = my_mean(ui, n, 1);
    Bj[id_odd] = my_mean(ui, n, 0);
  }
  Bj = Bj*sqrt(2)/sqrt(rt);
  double Bj2 = 0.0;
  for(int i = 0; i < k; i++){
    Bj2 += pow(Bj[i], 2);
  }
  return Bj2;
}

// define the projection penalty
// [[Rcpp::export]]
NumericVector proj_pen_C(NumericVector X, double a, double b, NumericVector Mj,
                         double kappa){
  int n = X.size(), n_m = Mj.size();
  double rt = b - a;
  double kappa1 = kappa/(n*rt);
  NumericVector X_tt = 2.0*M_PI*(X - a)/rt;
  // obtain the penalty
  NumericVector loss(n_m);
  for(int i = 0; i < n_m; i++){
    loss[i] = (-1.0)*B2_coef(X_tt, rt, Mj[i], n) + kappa1*(2*Mj[i] + 1);
  }
  return loss;
}



