rtnormInvR <- function(n, lower = -Inf, upper = Inf, mean = 0, sd = 1){
  if (lower == -Inf){
      pl <- 0
  } else {
      pl <- pnorm(lower, mean, sd)
  }
  if (upper == Inf){
      pu <- 1
  } else {
      pu <- pnorm(upper, mean, sd)
  }
  u <- runif(n, pl, pu)
  return(qnorm(u, mean, sd))
}

rtnormR <- function(n, lower = -Inf, upper = Inf, mean = 0, sd = 1){
  rand <- rtnormInvR(n, lower, upper, mean, sd)
  return(rand)
}
