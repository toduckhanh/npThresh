#' @import stats
#' @import utils

#' @export
YI_kern <- function(X1, X2, h1, h2, type_ker){
  ff <- function(tau, X1, X2, h1, h2, type_ker){
    return(Sen_kern(X2, tau, h2, type_ker) + Spe_kern(X1, tau, h1, type_ker))
  }
  a1 <- min(c(X1, X2))
  a2 <- max(c(X1, X2))
  a0 <- mean(mean(X1), mean(X2))
  out <- optim(a0, ff, method = "L-BFGS-B", lower = a1, upper = a2,
               control = list(fnscale = -1), X1 = X1, X2 = X2, h1 = h1,
               h2 = h2, type_ker = type_ker)$par
  return(out)
}

## Youden index approach for Normal distribution
#' @export
YI_normal <- function(mu, sigma, crit.var = 0.01){
  sigma2 <- sigma^2
  sigma_rt <- sigma2[1]/sigma2[2]
  var.check <- abs(sigma_rt - 1) < crit.var
  if(var.check){
    cpt <- (mu[1] + mu[2])/2
  } else{
    dif_var <- sigma2[1] - sigma2[2]
    cpt <- ((mu[2]*sigma2[1] - mu[1]*sigma2[2]) -
              sigma[1]*sigma[2]*sqrt((mu[1] - mu[2])^2 +
                                       dif_var*log(sigma_rt)))/dif_var
  }
  return(cpt)
}

