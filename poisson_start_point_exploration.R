library(robustbase)

n = 100 # sample size

# generate a dataset with p = 1 covariate + intercept

x = cbind(rep(1, n),
          rnorm(n, 0, 1))
beta = c(1, 0.5)
lambdas = exp(x %*% beta)
y = rpois(n, lambdas)

# plot data
plot(x[, 2], y)

n_start = 200 # number of starting points to try

# list of starting points to try
start_list = cbind(rep(0, n_start), # Keep the intercept at 0 
                   seq(from = -10, to = 10, length.out = n_start))

# initialize empty vectors
est = 0 * start_list
cv = vector('logical', length = n_start)


for (i in 1:n_start) {
  start = t(start_list[i, ])
  
  fit = glmrob(y ~ 1 + I(x[, 2]), 
                       family = poisson(), 
                       method = "Mqle", 
                       start = start, 
                       control = glmrobMqle.control(maxit = 1000, 
                                                    tcc = 0.5))
  cv[i] = fit$converged # check for convergence
  est[i, ] = fit$coefficients # save coefficient
}


# plot
plot(start_list[, 2], est[, 2], col = 2 - cv, 
     pch = 19,
     xlab = "starting point", ylab = "estimate")
legend("topleft", 
       legend = c("converded", "didn't converge (after 1000 iter)"),
       col = c("black", "red"), pch = 19, bty = "n")