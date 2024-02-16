psi <- function(z, c_rob = 1.345) {
  out_up <- z > c_rob
  out_low <- z < -c_rob
  
  pz <- z
  pz[out_up]  <- c_rob
  pz[out_low] <- -c_rob
  
  return(pz)
}

# Fisher consistency constant
a_cst <- function(lambda, c_rob = 1.345) {
  
  l_up <- ceiling(lambda + c_rob * sqrt(lambda)) - 1
  l_lo <- floor(lambda - c_rob * sqrt(lambda))
  
  a_1 <- (1 - ppois(l_lo, lambda) - ppois(l_up, lambda)) * c_rob * sqrt(lambda)
  a_2 <- (dpois(l_lo, lambda) - dpois(l_up, lambda)) * lambda
  
  return(a_1 + a_2)
  #return(0)
}

# lambda : vector of size m 
# y      : vector of size n

# score_poisson <- function(lambda, y, c_rob = 1.345) {
#   
#   y <- as.vector(y)
#   lambda <- as.vector(lambda)
#   n <- length(y)
#   m <- length(lambda)
#   sc <- numeric(m)
#   for (i in 1:m) {
#     res <- (y - lambda[i]) / sqrt(lambda[i])
#     sc[i] <- sum(sqrt(lambda[i]) * psi(res, c_rob)) - a_cst(lambda[i])
#   }
#   return(sc - n * a_cst(lambda, c_rob = c_rob))
# }

score_poisson <- function(beta, y, x = NULL, c_rob = 1.345) {
  
  y <- as.vector(y)
  beta <- as.matrix(beta)
  n <- length(y)
  
  if (is.null(x)) {
    x <- rep(1, n)
    print("no x")
  }
  
  x <- as.matrix(x)
  p <- ncol(x)
  m <- ncol(beta)
  
  sc <- matrix(ncol = m, nrow = p)
  
  for (i in 1:m) {
    beta_tmp = c(beta[, i])
    loglambda <- x %*% beta_tmp
    lambda <- c(exp(loglambda ))
    res <- (y - lambda) / sqrt(lambda)
    s <- (sqrt(lambda) * psi(res, c_rob = c_rob) - a_cst(lambda, c_rob)) %x% t(rep(1, p))
    sc[, i] <- colSums(s * x)
  }
  return(sc)
}