# test 1 ----
x <- rnorm(100, mean = 0, sd = 0.25)
x_grid <- seq(-1, 1, length.out = 151)

system.time({
  out1 <- npThresh:::cdf_proj_C(x = x_grid, X = x, a = min(x), b = max(x), M = m)
})

system.time({
  out2 <- cdf_proj(x = x_grid, X = x, M = m, ab = range(x))
})

x_F <- ecdf(x)
plot(x_grid, x_F(x_grid), type = "l", col = "gray60")
lines(x_grid, pnorm(x_grid, mean = 0, sd = 1), col = "green")
lines(x_grid, out2, col = "blue")

system.time({
  out3 <- npThresh:::proj_pen_C(X = x, a = min(x), b = max(x), Mj = c(1:100),
                                kappa = 2)
})
out3

plot(c(1:100), out3, type = "l")
kappa_grid <- seq(0.1, 2, by = 0.1)

system.time({
  out3 <- sapply(kappa_grid, function(t){
    npThresh:::proj_pen_C(X = x, a = min(x), b = max(x), Mj = c(1:50),
                          kappa = t)
    })
})
any(out3 > 0)


id <- apply(out3, 2, which.min)
id
out3[id[9], 9]
plot(c(1:50), out3[, 9], type = "l")

kappa_grid[9]

system.time({
  out2_1 <- cdf_proj(x = x_grid, X = x, M = 2, ab = range(x))
})

x_F <- ecdf(x)
plot(x_grid, x_F(x_grid), type = "l", col = "gray60")
lines(x_grid, pnorm(x_grid, mean = 0, sd = 0.25), col = "green")
lines(x_grid, out2_1, col = "blue")


a <- min(x)
b <- max(x)

m <- 40
bj <- numeric(2*m + 1)
bj[1] <- 1/sqrt(2)

for(i in 1:m){
  uu <- 2*pi*i*(x - a)/(b - a)
  bj[2*i] <- mean(cos(uu))
  bj[2*i + 1] <- mean(sin(uu))
}
bj <- bj*sqrt(2)/sqrt(b - a)

-sum(bj^2) + 1.5*(2*m + 1)/(100*(b - a))


# test CV ----
x <- rnorm(60, mean = 0, sd = 1)
x_grid <- seq(-4, 4, length.out = 151)
m <- c(1:50)

system.time({
  out_cv_proj <- sapply(m, function(y){
    out2 <- sapply(c(1:100), function(i){
      ((x_grid - x[i] >= 0) - cdf_proj(x = x_grid, X = x[-i], M = y,
                                       ab = range(x)))^2
    })
    return(mean(sapply(c(1:100), function(i){
      caTools::trapz(x_grid, out2[,i])
    })))
  })
})
plot(m, out_cv_proj, type = "l")
m[which.min(out_cv_proj)]

system.time({
  out2_1 <- cdf_proj(x = x_grid, X = x, M = m[which.min(out_cv_proj)], ab = range(x))
})

x_F <- ecdf(x)
plot(x_grid, x_F(x_grid), type = "l", col = "gray60")
lines(x_grid, pnorm(x_grid, mean = 0, sd = 1), col = "green")
lines(x_grid, out2_1, col = "blue")

id_5 <- matrix(sample(x = 1:length(x), size = length(x), replace = FALSE),
               nrow = 20, ncol = 5)
n <- length(x)
k_fold <- 10
nk_fold <- ceiling(n/k_fold)
id_cv <- lapply(1:(k_fold - 1), function(i){
  (1 + (i - 1)*nk_fold) : (i*nk_fold)
})
id_cv[[k_fold]] <- (1 + (k_fold - 1)*nk_fold) : n

system.time({
  out_cv_5_proj <- sapply(m, function(y){
    out2 <- sapply(1:k_fold, function(i){
      temp1 <- sapply(x_grid, function(u) mean(u - x[id_cv[[i]]] >= 0))
      temp2 <- cdf_proj(x = x_grid, X = x[-id_cv[[i]]], M = y, ab = range(x))
      temp3 <- (temp1 - temp2)^2
      return(caTools::trapz(x_grid, temp3))
    })
    return(mean(out2))
  })
})
plot(m, out_cv_5_proj, type = "l")
m[which.min(out_cv_5_proj)]


system.time({
  out2_2 <- cdf_proj(x = x_grid, X = x, M = m[which.min(out_cv_5_proj)],
                     ab = range(x))
})

x_F <- ecdf(x)
plot(x_grid, x_F(x_grid), type = "l", col = "gray60")
lines(x_grid, pnorm(x_grid, mean = 0, sd = 1), col = "green")
lines(x_grid, out2_1, col = "blue")
lines(x_grid, out2_2, col = "red")

n <- 68
k_fold <- 10
nk_fold <- ceiling(n/k_fold)
id_cv <- lapply(1:(k_fold - 1), function(i){
  (1 + (i - 1)*nk_fold) : (i*nk_fold)
})
id_cv[[k_fold]] <- (1 + (k_fold - 1)*nk_fold) : n

sapply(id_cv, length)

system.time({
  out_kcv_proj <- kcv_proj(X = x, k_fold = 5)
})
which.min(out_kcv_proj$cv_out)
plot(out_kcv_proj$M, out_kcv_proj$cv_out, type = "l")

system.time({
  out2_2 <- cdf_proj(x = x_grid, X = x, M = out_kcv_proj$m_opt,
                     ab = range(x))
})
