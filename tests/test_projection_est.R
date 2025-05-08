# test 1 ----
x <- rnorm(100)
x_grid <- seq(-4, 4, length.out = 151)
m <- 10

a <- min(x)
b <- max(x)

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

bj <- numeric(2*m + 1)
bj[1] <- 1/sqrt(b - a)

for(i in 1:m){
  uu <- 2*pi*i*(x - a)/(b - a)
  bj[2*i] <- mean(cos(uu))
  bj[2*i + 1] <- mean(sin(uu))
}
all.equal(out1, bj)

bj[id_even] <- colMeans(matrix(cos(2*pi*(1:length(id_even))*(x - a)/(b - a)),
                               ncol = 10, byrow = TRUE))

mean(sin(2*pi*(x - a)/(b - a)))

mean(cos(2*pi*2*(x - a)/(b - a)))
mean(sin(2*pi*2*(x - a)/(b - a)))
