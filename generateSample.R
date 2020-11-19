
library(mvtnorm)
library(abind)  # combine matrices or list to a array

X.phi <- function(NumberOfphi, xgrid) {
  phi <- matrix(1, NumberOfphi, length(xgrid))
  for (j in 2 : NumberOfphi)
    phi[j, ] <- 4 * (2 * j - 1)^(-1) * sqrt(2) * cos((j-1) * pi * xgrid)
  return(phi)
}

Wfun <- function(NumberOfphi, xgrid){
  xsi <- rnorm(NumberOfphi, 0, 1)
  W <- colSums(xsi * X.phi(NumberOfphi, xgrid)) 
  return(W)
}

#################################################################
############# generate data for p = 4 ###########
geSample <- function(n, ygrid, xgrid, sigma, rho, tau){
  NumberOfphi = 50
  data <- list()
  w1 <- t(replicate(n, Wfun(NumberOfphi, xgrid))) # n * length(xgrid)
  w2 <- t(replicate(n, Wfun(NumberOfphi, xgrid)))
  w3 <- t(replicate(n, Wfun(NumberOfphi, xgrid)))
  w4 <- t(replicate(n, Wfun(NumberOfphi, xgrid)))
  x1 <- w1 + tau * (w2 + w3) # n * length(xgrid)
  x2 <- w2 + tau * (w1 + w3)
  x3 <- w3 + tau * (w1 + w2)
  x4 <- w4
  X <- abind(x1, x2, x3, x4, along = 3)
  data$X <- X
  
  beta1Par <- matrix(0, 4, 4)
  for(k in 1 : 4){
    for(l in 1 : 4){
      beta1Par[k, l] = 0.1 * (k + l)
    }
  }
  betaMa1 <- t(X.phi(4, ygrid)) %*% beta1Par %*% X.phi(4, xgrid)  # length(ygrid) * length(xgrid)
  
  beta2Par <- matrix(0, 50, 50)
  for(k in 1 : 50){
    for(l in 1 : 50){
      beta2Par[k, l] = 2 * (-1)^(k + l) * k^(-1) * l^(-2)
    }
  }
  betaMa2 <- t(X.phi(50, ygrid)) %*% beta2Par %*% X.phi(50, xgrid)   # length(ygrid) * length(xgrid)
  betaMa3 <- matrix(0, length(ygrid), length(xgrid))
  betaMa4 <- matrix(0, length(ygrid), length(xgrid))
  betaMa <- abind(betaMa1, betaMa2, betaMa3, betaMa4, along = 3)
  data$beta <- betaMa
  
  fx1 <- x1 %*% t(betaMa1) * (xgrid[2] - xgrid[1])
  fx2 <- x2 %*% t(betaMa2) * (xgrid[2] - xgrid[1])
  fx3 <- x3 %*% t(betaMa3) * (xgrid[2] - xgrid[1])
  fx4 <- x4 %*% t(betaMa4) * (xgrid[2] - xgrid[1])
  
  SigMa <- matrix(0, length(ygrid), length(ygrid))
  for (i in 1 : length(ygrid)){
    for (j in 1 : length(ygrid)){
      SigMa[i, j] <- sigma^2 * rho^(10 * abs(ygrid[i] - ygrid[j]))
    }
  }
  eps <- rmvnorm(n, mean = rep(0, length(ygrid)), sigma = SigMa)
  y <- fx1 + fx2 + fx3 + fx4 + eps
  data$Y <- y
  data$ygrid <- ygrid
  data$xgrid <- xgrid
  data$sigma2 <- sigma^2
  data$tau <- tau
  return(data)
}


#################################################################
############# generate data for relatively large p ###########

geSample2 <- function(n, ygrid, xgrid, p, sigma, rho, tau){
  NumberOfphi = 50
  data <- list()
  w.array <- array(NA, c(n, length(xgrid), p + 2))
  for (j in 1 : (p +2)){
    w.array[ , , j] <- t(replicate(n, Wfun(NumberOfphi, xgrid)))
  }
  x.array <- array(NA, c(n, length(xgrid), p))
  for (j in 2 : 11){
    x.array[ , , j-1] <- w.array[ , , j] + tau * (w.array[ , , j-1] + w.array[ , , j+1])
  }
  for (j in 11 : p) {
    x.array[ , , j] <- w.array[ , , j + 2]
  }
  data$X <- x.array
  
  betaPar <- array(0, c(50, 50, p))
  for (j in 1 : 6){
    for (k in 1 : 50) {
      for (l in 1 : 50){
        betaPar[k, l, j] <- (k + j)^(-1) * l^(-1) * (-1) ^(k + l)
      }
    }
  }
  beta.array <- array(NA, c(length(ygrid), length(xgrid), p))
  fx.array <- array(NA, c(n, length(ygrid), p))
  for (j in 1 : p){
    beta.array[ , , j] <- t(X.phi(50, ygrid)) %*% betaPar[ , , j] %*% X.phi(50, xgrid)  
    fx.array[ , , j] <- x.array[ , , j] %*% t(beta.array[ , , j]) * (xgrid[2] - xgrid[1])   
  }
  data$beta <- beta.array
  
  SigMa <- matrix(0, length(ygrid), length(ygrid))
  for (i in 1 : length(ygrid)){
    for (j in 1 : length(ygrid)){
      SigMa[i, j] <- sigma^2 * rho^(10 * abs(ygrid[i] - ygrid[j]))
    }
  }
  eps <- rmvnorm(n, mean = rep(0, length(ygrid)), sigma = SigMa)
  y <- apply(fx.array, c(1, 2), sum) + eps
  data$Y <- y
  data$ygrid <- ygrid
  data$xgrid <- xgrid
  data$sigma2 <- sigma^2
  data$tau <- tau
  return(data)
}
