YI_logGamma <- function(mu, sigma, start = NULL, method.optim = "L-BFGS-B",
                        maxit = 200, lower = 0, upper = Inf){
  if(is.null(start)) start <- mean(mu)
  YI_logGamma_ff <- function(par, mu, sigma){
    sp <- plnorm(par, meanlog = mu[1], sdlog = sigma[1])
    se <- 1 - pgamma(par, shape = mu[2], rate = sigma[2])
    return(se + sp - 1)
  }
  YI_optim <- optim(par = start, fn = YI_logGamma_ff, mu = mu, sigma = sigma,
                    method = method.optim,
                    control = list(maxit = maxit, fnscale = -1),
                    lower = lower, upper = upper)
  return(YI_optim)
}

CtP_logGamma <- function(mu, sigma, start = NULL, method.optim = "L-BFGS-B",
                         maxit = 200, lower = 0, upper = Inf){
  if(is.null(start)) start <- mean(mu)
  CtP_logGamma_ff <- function(par, mu, sigma){
    sp <- plnorm(par, meanlog = mu[1], sdlog = sigma[1])
    se <- 1 - pgamma(par, shape = mu[2], rate = sigma[2])
    return((1 - sp)^2 + (1 - se)^2)
  }
  CtP_optim <- optim(par = start, fn = CtP_logGamma_ff, mu = mu, sigma = sigma,
                    method = method.optim,
                    control = list(maxit = maxit),
                    lower = lower, upper = upper)
  return(CtP_optim)
}

MA_logGamma <- function(mu, sigma, start = NULL, method.optim = "L-BFGS-B",
                         maxit = 200, lower = 0, upper = Inf){
  if(is.null(start)) start <- mean(mu)
  MA_logGamma_ff <- function(par, mu, sigma){
    sp <- plnorm(par, meanlog = mu[1], sdlog = sigma[1])
    se <- 1 - pgamma(par, shape = mu[2], rate = sigma[2])
    return(log(sp) + log(se))
  }
  MA_optim <- optim(par = start, fn = MA_logGamma_ff, mu = mu, sigma = sigma,
                     method = method.optim,
                     control = list(maxit = maxit, fnscale = -1),
                     lower = lower, upper = upper)
  return(MA_optim)
}


mu_logG_true_2 <- c(0, 3.6)
sigma_logG_true_2 <- c(0.5, 1)

tau_logG_true_2 <- YI_logGamma(mu = mu_logG_true_2,
                               sigma = sigma_logG_true_2)$par
sp_logG_true_2 <- plnorm(tau_logG_true_2, meanlog = mu_logG_true_2[1],
                         sdlog = sigma_logG_true_2[1])
se_logG_true_2 <- 1 - pgamma(tau_logG_true_2, shape = mu_logG_true_2[2],
                             rate = sigma_logG_true_2[2])

tau_CtP_logG_true_2 <- CtP_logGamma(mu = mu_logG_true_2,
                                    sigma = sigma_logG_true_2)$par
sp_CtP_logG_true_2 <- plnorm(tau_CtP_logG_true_2, meanlog = mu_logG_true_2[1],
                             sdlog = sigma_logG_true_2[1])
se_CtP_logG_true_2 <- 1 - pgamma(tau_CtP_logG_true_2, shape = mu_logG_true_2[2],
                                rate = sigma_logG_true_2[2])

tau_MA_logG_true_2 <- MA_logGamma(mu = mu_logG_true_2,
                                  sigma = sigma_logG_true_2)$par
sp_MA_logG_true_2 <- plnorm(tau_MA_logG_true_2, meanlog = mu_logG_true_2[1],
                             sdlog = sigma_logG_true_2[1])
se_MA_logG_true_2 <- 1 - pgamma(tau_MA_logG_true_2, shape = mu_logG_true_2[2],
                                rate = sigma_logG_true_2[2])


n0 <- n1 <- 200
X1 <- rlnorm(n0, meanlog = mu_logG_true_2[1], sdlog = sigma_logG_true_2[1])
X2 <- rgamma(n1, shape = mu_logG_true_2[2], rate = sigma_logG_true_2[2])

h1_RT <- bwd_plugin(x = X1, bwd_method = "RT")
h1_AZZ <- bwd_plugin(x = X1, bwd_method = "AZZ")

h2_RT <- bwd_plugin(x = X2, bwd_method = "RT")
h2_AZZ <- bwd_plugin(x = X2, bwd_method = "AZZ")

system.time({
  t_YI_RT <- YI_kern(X1 = X1, X2 = X2, h1 = h1_RT, h2 = h2_RT, type_ker = 0)
})

cdf_kernel(x = t_YI_RT, X = X1, kernel_type = 0, bwd = h1_RT)
1 - cdf_kernel(x = t_YI_RT, X = X2, kernel_type = 0, bwd = h2_RT)

system.time({
  t_YI_AZZ <- YI_kern(X1 = X1, X2 = X2, h1 = h1_AZZ, h2 = h2_AZZ, type_ker = 0)
})

cdf_kernel(x = t_YI_AZZ, X = X1, kernel_type = 0, bwd = h1_AZZ)
1 - cdf_kernel(x = t_YI_AZZ, X = X2, kernel_type = 0, bwd = h2_AZZ)


t_CtP_RT <- CtP_kern(X1 = X1, X2 = X2, h1 = h1_RT, h2 = h2_RT, type_ker = 0)
cdf_kernel(x = t_CtP_RT, X = X1, kernel_type = 0, bwd = h1_RT)
1 - cdf_kernel(x = t_CtP_RT, X = X2, kernel_type = 0, bwd = h2_RT)

t_CtP_AZZ <- CtP_kern(X1 = X1, X2 = X2, h1 = h1_AZZ, h2 = h2_AZZ, type_ker = 0)
cdf_kernel(x = t_CtP_AZZ, X = X1, kernel_type = 0, bwd = h1_AZZ)
1 - cdf_kernel(x = t_CtP_AZZ, X = X2, kernel_type = 0, bwd = h2_AZZ)

t_MA_RT <- MA_kern(X1 = X1, X2 = X2, h1 = h1_RT, h2 = h2_RT, type_ker = 0)
t_MA_AZZ <- MA_kern(X1 = X1, X2 = X2, h1 = h1_AZZ, h2 = h2_AZZ, type_ker = 0)

cdf_kernel(x = t_MA_RT, X = X1, kernel_type = 0, bwd = h1_RT)
1 - cdf_kernel(x = t_MA_RT, X = X2, kernel_type = 0, bwd = h2_RT)

