#' @export
empi_llike <- function(se, sp, tau, X1, X2){
  ll <- Inf
  a1 <- max(min(X1), min(X2))
  a2 <- min(max(X1), max(X2))
  if(tau >= a1 & tau < a2){
    Fx_tau <- mean(X1 <= tau)
    Fy_tau <- mean(X2 <= tau)
    ll1 <- 2*length(X1)*(Fx_tau*log(Fx_tau/sp) +
                           (1 - Fx_tau)*log((1 - Fx_tau)/(1 - sp)))
    ll2 <- 2*length(X2)*(Fy_tau*log(Fy_tau/(1 - se)) +
                           (1 - Fy_tau)*log((1 - Fy_tau)/se))
    ll <- ll1 + ll2
  }
  return(ll)
}
