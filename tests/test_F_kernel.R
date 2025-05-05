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
x_grid <- seq(min(x), max(x), length.out = 101)

system.time({
  cdf_locv <- sapply(1:length(x), function(y){
    cdf_kernel(x = x_grid, X = x[-y], kernel_type = 0, bwd = 0.32)
  })
})

out_trap <- sapply(1:length(x), function(i){
  trapz(x = x_grid, y = ((x_grid - x[i] >= 0) - cdf_locv[, i])^2)
})
mean(out_trap)

h_grid <- seq((max(x) - min(x))/200, (max(x) - min(x))/2, length.out = 50)
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

