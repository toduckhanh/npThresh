# test 1 ----
x <- rnorm(100)
x_F <- ecdf(x)
x_grid <- seq(-1, 2, by = 0.01)

plot(x_grid, x_F(x_grid), type = "b")

system.time({
  u1 <- cdf_kernel(x = x_grid, X = x, kernel_type = 0, bwd = 0.1)
})

system.time({
  u2 <- sapply(x_grid, function(y) mean(pnorm((y - x)/0.1)))
})

all.equal(u1, u2)

system.time({
  u3 <- cdf_kernel(x = x_grid, X = x, kernel_type = 0, bwd = 0.1589813)
})
u3

# test leave-one-out cdf ----
require(caTools)
x <- rnorm(100)
x_grid <- seq(min(x), max(x), length.out = 151)

# system.time({
#   w_K <- npThresh:::cdf_loocv_kernel_C(x = x_grid, X = x, Ktype = 0, bwd = 0.1)
#   nF_kernel <- rowSums(w_K)
#   cv_est <- sapply(1:length(x), function(i){
#     nF_kernel - w_K[,i]
#   })/(length(x) - 1)
# })
# dim(cv_est)

system.time({
  cv_out <- npThresh:::cdf_loocv_kernel_C(x = x_grid, X = x, Ktype = 0,
                                          bwd = 0.1, hx = diff(x_grid)[1])
})

cv_out

out_trap <- sapply(1:length(x), function(i){
  cdf_locv <- cdf_kernel(x = x_grid, X = x[-i], kernel_type = 0, bwd = 0.1)
  res <- trapz(x = x_grid, y = ((x_grid - x[i] >= 0) - cdf_locv)^2)
})
mean(out_trap)


simp_int <- function (x, fx, n.pts = 256, ret = FALSE){
  if (class(fx) == "function") fx = fx(x)
  n.x = length(x)
  if (n.x != length(fx)) stop("Unequal input vector lengths")
  if (n.pts < 64) n.pts = 64
  ap = approx(x, fx, n = 2 * n.pts + 1)
  h = diff(ap$x)[1]
  integral = h * (ap$y[2 * (1:n.pts) - 1] + 4 * ap$y[2 * (1:n.pts)] +
                    ap$y[2 * (1:n.pts) + 1])/3
  return(list(value = sum(integral),
              cdf = list(x = ap$x[2 * (1:n.pts)], y = cumsum(integral))))
}

fy_test <- ((x_grid - x[5] >= 0) - cv_out[,5])^2
npThresh:::simpson(fx = fy_test, nx = 100, hx = diff(x_grid)[1])
trapz(x = x_grid, fy_test)
simp_int(x = x_grid, fx = fy_test, n.pts = 101)$value

# cdf_kernel(x = x_grid, X = x, kernel_type = 0, bwd = 0.1)

tt <- (1:1001)*pi/1000
f_tt <- sin(tt)

npThresh:::simpson(fx = f_tt, nx = 1000, hx = diff(tt)[1])
trapz(tt, f_tt)
simp_int(x = tt, fx = f_tt, n.pts = 1001)$value

system.time({
  cdf_locv <- sapply(1:length(x), function(y){
    cdf_kernel(x = x_grid, X = x[-y], kernel_type = 0, bwd = 0.1)
  })
})

cdf_locv

all.equal(cv_est, cdf_locv)
all.equal(cv_out, cdf_locv)

out_trap <- sapply(1:length(x), function(i){
  trapz(x = x_grid, y = ((x_grid - x[i] >= 0) - cdf_locv[, i])^2)
})
mean(out_trap)

h_grid <- seq((max(x) - min(x))/200, (max(x) - min(x))/2, length.out = 101)
system.time({
  cv_out_bwd <- sapply(h_grid, function(h){
    cv_out <- npThresh:::cdf_loocv_kernel_C(x = x_grid, X = x, Ktype = 0,
                                            bwd = h, hx = diff(x_grid)[1])
    return(cv_out)
  })
})
plot(h_grid, cv_out_bwd)
h_grid[which.min(cv_out_bwd)]

system.time({
  out_cv_bwd <- npThresh:::cv_bwd_C(x = x_grid, X = x, Ktype = 0,
                                    hx = diff(x_grid)[1], bwd_pts = 101)
})

out_cv_bwd$bwd_seq
out_cv_bwd$cv_out
out_cv_bwd$h_opt

system.time({
  out_cv_bwd_f <- bwd_cv(X = x, kernel_type = "gauss")
})
out_cv_bwd_f$h_opt

system.time({
  out_cv <- sapply(h_grid, function(h){
    out_trap <- sapply(1:length(x), function(i){
       cdf_locv <- cdf_kernel(x = x_grid, X = x[-i], kernel_type = 0, bwd = h)
       res <- trapz(x = x_grid, y = ((x_grid - x[i] >= 0) - cdf_locv)^2)
    })
    return(mean(out_trap))
  })
})

plot(h_grid, out_cv)
h_grid[which.min(out_cv)]

system.time({
  out_cv <- sapply(h_grid, function(h){
    cdf_locv <- sapply(1:length(x), function(y){
      cdf_kernel(x = x_grid, X = x[-y], kernel_type = 0, bwd = h)
    })
    out_trap <- sapply(1:length(x), function(i){
      trapz(x = x_grid, y = ((x_grid - x[i] >= 0) - cdf_locv[, i])^2)
    })
    return(mean(out_trap))
  })
})

plot(h_grid, out_cv)
h_grid[which.min(out_cv)]

require(kerdiest)
system.time({
  out_kerd <- CVbw(type_kernel = "n", vec_data = x)
})

# LOOCV bandwidth ----
x <- rnorm(100)
system.time({
  out_cv_bwd_f <- bwd_cv(X = x, kernel_type = "gauss")
})
out_cv_bwd_f$h_opt
bwd_plugin(X = x, bwd_method = "RT")
bwd_plugin(X = x, bwd_method = "AZZ")

x_F <- ecdf(x)
x_grid <- seq(min(x), max(x), by = 0.01)

u1 <- cdf_kernel(x = x_grid, X = x, kernel_type = "gauss",
                 bwd = out_cv_bwd_f$h_opt)
u2 <- cdf_kernel(x = x_grid, X = x, kernel_type = "gauss",
                 bwd = bwd_plugin(X = x, bwd_method = "RT"))
u3 <- cdf_kernel(x = x_grid, X = x, kernel_type = "gauss",
                 bwd = bwd_plugin(X = x, bwd_method = "AZZ"))

plot(x_grid, x_F(x_grid), type = "l")
lines(x_grid, u1, col = "blue")
lines(x_grid, pnorm(x_grid), col = "green")
lines(x_grid, u2, col = "red", lty = 2)
lines(x_grid, u3, col = "pink", lty = 2)

