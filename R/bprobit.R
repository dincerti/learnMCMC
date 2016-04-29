bprobitR <- function(sims = 1000, x, y, b0, B0inv, beta.start){
  # initialize
  n <- length(y)
  n0 <- sum(y == 0)
  n1 <- sum(y == 1)
  ystar <- rep(0, n)
  beta.store <- matrix(NA, nrow = sims, ncol = ncol(x))
  colnames(beta.store) <- colnames(x)

  # starting values
  beta <- beta.start

  # crossproducts
  B <- solve(B0inv + crossprod(x, x))

  # gibbs sampler
  for (i in 1:sims){
    # draw ystar
    ystar[y == 0] <- rtnormR(n0, upper = 0, mean = x[y == 0, ] %*% beta)
    ystar[y==1] <- rtnormR(n1, lower = 0, mean = x[y == 1, ] %*% beta)

    # draw beta
    b <- B %*% (B0inv %*% b0  + crossprod(x, ystar))
    beta <- t(rmvnorm(1, b, B))

    # store
    beta.store[i, ] <- beta
  }
  return(list(beta = beta.store))
}
