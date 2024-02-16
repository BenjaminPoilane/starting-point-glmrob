rm(list = ls())
library(robustbase)
library(ggplot2)

source("poissonglm_function.R")
set.seed(16022024) # for reproducibility

# sample size
n <- 20
# generative parameter
beta_0 <- 1

# generate data
x <- round(rnorm(n), digit = 2)
lambdas <- exp(beta_0 * x) # no intercept
y <- rpois(n, lambdas)

# contaminated data
# xcont <- c(x, -4, 4)
# ycont <- c(y, 30, 0)
xcont <- c(x, -3.5, 3.5)
ycont <- c(y, 25, 0)
# plot data
plot(xcont, ycont, 
     xlab = "x", ylab = "y", 
     pch = 16, col = c(rep(1, n), rep(2, 2)), 
     ylim = c(0, 30))



# add the y = exp(beta_0 * x) line
tplot <- seq(-5, 5, length.out = 200)
points(tplot, exp(tplot), type = "l", col = "grey50", lty = 2)




bplot <- seq(-1.1, 0.8, by = 0.005)
plot(bplot, 
     score_poisson(t(bplot), ycont, xcont, c_rob = 1.345), 
     type = "l", 
     xlab = "beta", 
     ylab = "score function", 
     col = "blue", 
     main = "robust poisson glm (c_rob = 1.345)")
abline(h = 0)


#
fit0 <- glmrob(ycont ~ xcont - 1, family = poisson())
fit1 <- glmrob(ycont ~ xcont - 1, family = poisson(), start = 1)

abline(v = fit0$coefficients, lty = 2)
abline(v = fit1$coefficients, lty = 2)
text(x = fit0$coefficients, y = 40, paste("sum w =", signif(sum(fit0$w.r), 2)), pos = 4)
text(x = fit1$coefficients, y = 40, paste("sum w =", signif(sum(fit1$w.r), 2)), pos = 4)
sum(fit1$w.r)
