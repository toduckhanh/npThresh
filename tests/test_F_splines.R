library(tidyverse)
library(splines)

x <- rnorm(30)
ecdf_x <- ecdf(x)
y <- ecdf_x(x)

qlogis(y)

data_test <- tibble(x = x, y = y)
ggplot(data_test, mapping = aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_bw()

knots_x <- quantile(data_test$x, probs = c(0.25, 0.5, 0.75))
cdf_x_spline <- lm(y ~ bs(x, knots = knots_x, degree = 3), data = data_test)
summary(cdf_x_spline)


predict(cdf_x_spline)

ggplot(data_test, mapping = aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm", 
              formula = y ~ bs(x, knots = knots_x, degree = 3)) +
  theme_bw()

