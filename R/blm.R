blmR <- function(sims = 1000, x, y, b0, B0inv, c0, d0, beta.start, sigma2.start){
  # initialize 
  n <- length(y)
  sigma2.store <- rep(NA, sims)
  beta.store <- matrix(NA, nrow = sims, ncol = ncol(x))
  colnames(beta.store) <- colnames(x)
  
  # starting values
  beta <- beta.start
  sigma2 <- sigma2.start
  
  # crossproducts
  xtx <- crossprod(x, x)
  xty <- crossprod(x, y)
  
  # gibbs sampler
  for (i in 1:sims){
    
    # draw beta 
    B <- solve(B0inv + xtx/sigma2)
    b <- B %*% (B0inv %*% b0  + xty/sigma2)
    beta <- rmvnorm(1, b, B)
    
    # draw sigma2
    epsilon <- y - tcrossprod(x, beta)
    sigma2 <- 1/rgamma(1, shape = (n + c0)/2, 
                       rate = (d0 + crossprod(epsilon, epsilon))/2)
    
    # store
    beta.store[i, ] <- beta
    sigma2.store[i] <- sigma2
  }
  return(list(beta = beta.store, sigma2 = sigma2.store))
}  