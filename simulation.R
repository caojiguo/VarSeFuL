
library(fda.usc)
library(grpreg)

source("CVAIC.R")
source("estimation.R")

##############################################################
###### proposed method #######
##### common and fixed truncation parameters for predictors  #####
simulation1 <- function(data, testdat, Xtrun){
  ygrid <- data[, 1]$ygrid
  xgrid <- data[, 1]$xgrid
  betaTrue <- data[, 1]$beta
  R <- dim(data)[2]
  ISE <- rep(NA, R)
  PSR <- rep(NA, R)  ## positive selection rate of nonzeoro coefficients
  NSR <- rep(NA, R)   ## noncausal selection rate of zeoro coefficients
  RPE <- rep(NA, R)
  for (r in 1 : R){
    Y <- data[, r]$Y
    X <- data[, r]$X
    k <- iniTrun(Y, X, ygrid, xgrid)[1]
    trun <- c(2, Xtrun)  ## Xtrun is the given parameters of the predcitors
    result <- GCDfun(Y, X, ygrid, xgrid, trun)
    estbeta <- result$estbeta
    ISE[r] <- sum((estbeta - betaTrue)^2) * (xgrid[2] - xgrid[1]) * (ygrid[2] - ygrid[1])
    est_ind <- apply(abs(estbeta), 3, sum)
    tru_ind <- apply(abs(betaTrue), 3, sum) # true index of zero and nonzero beta
    nzero_loca <- which(tru_ind != 0)   # true locations of nonzero beta
    zero_loca <- which(tru_ind == 0)  # true locations of zero beta
    PSR[r] <- sum(est_ind[nzero_loca] != 0)/length(nzero_loca)
    NSR[r] <- sum(est_ind[zero_loca] == 0)/length(zero_loca)
    
    Ytest <- testdat[, r]$Y
    Xtest <- testdat[, r]$X
    p <- dim(Xtest)[3]
    fxtest <- array(NA, c(dim(Ytest), p))
    for (j in 1 : p){
      fxtest[ , , j] <- Xtest[ , , j] %*% t(estbeta[ , , j]) * (xgrid[2] - xgrid[1])
    }
    Ypred <- apply(fxtest, c(1, 2), sum)
    RPE[r] <- mean(apply((Ytest - Ypred)^2, 1, sum) * (ygrid[2] - ygrid[1])/
                     (apply(Ytest^2, 1, sum) * (ygrid[2] - ygrid[1])))
  }
  return(list(PSR = PSR, NSR = NSR, ISE = ISE, RPE = RPE))
}

##### data-driven truncation parameters for predictors  #####
simulation <- function(data, testdat){
  ygrid <- data[, 1]$ygrid
  xgrid <- data[, 1]$xgrid
  betaTrue <- data[, 1]$beta
  p <- dim(betaTrue)[3]
  R <- dim(data)[2]
  ISE <- rep(NA, R)
  PSR <- rep(NA, R)  ## positive selection rate of nonzeoro coefficients
  NSR <- rep(NA, R)   ## noncausal selection rate of zeoro coefficients
  RPE <- rep(NA, R)
  for (r in 1 : R){
    Y <- data[, r]$Y
    X <- data[, r]$X
    tun_result <- CVAIC(Y, X, ygrid, xgrid)  ## select truncation parameters by the proposed method
    trun <- tun_result$optrun
    indent_ind <- tun_result$x_select  # 1 represents nonzero, 0 represents zero 
    indent_nzloca <- tun_result$nzloca    # identified locations for nonzero predictors
    X_nonzero <- X[ , , indent_nzloca]
    result <- GCDfun(Y, X_nonzero, ygrid, xgrid, trun)
    estbeta <- array(0, dim(betaTrue))
    estbeta[ , , indent_nzloca] <- result$estbeta
    ISE[r] <- sum((estbeta - betaTrue)^2) * (xgrid[2] - xgrid[1]) * (ygrid[2] - ygrid[1])
  
    tru_ind <- apply(abs(betaTrue), 3, sum)
    nzero_loca <- which(tru_ind != 0)   # true locations of nonzero predictors
    zero_loca <- which(tru_ind == 0)    # true locations of zero predictors
    est_ind <- apply(abs(estbeta), 3, sum) # finial select results, 0 represents zero coefficient, otherwise nonzero
    PSR[r] <- sum(est_ind[nzero_loca] != 0)/length(nzero_loca)
    NSR[r] <- sum(est_ind[zero_loca] == 0)/length(zero_loca)
    
    Ytest <- testdat[, r]$Y
    Xtest <- testdat[, r]$X
    fxtest <- array(NA, c(dim(Ytest), p))
    for (j in 1 : p){
      fxtest[ , , j] <- Xtest[ , , j] %*% t(estbeta[ , , j]) * (xgrid[2] - xgrid[1])
    }
    Ypred <- apply(fxtest, c(1, 2), sum)
    RPE[r] <- mean(apply((Ytest - Ypred)^2, 1, sum) * (ygrid[2] - ygrid[1])/
                     (apply(Ytest^2, 1, sum) * (ygrid[2] - ygrid[1])))
  }
  return(list(PSR = PSR, NSR = NSR, ISE = ISE, RPE = RPE))
}

######################################################################
######### ordinary least method for full model  ############
#####################################################################
#### AIC criterion for refitting and selecting truncation parameters of X
AIC_ls <- function(Y, X, ygrid, xgrid, klpar){
  k <- klpar[1] # klpar is the initial parameters of Y and X
  if(is.matrix(X)) X = array(X, c(dim(X), 1)) else X = X 
  d <- dim(X)[3]
  
  trun_list <- list()
  trun_list[[1]] <- k
  for (j in 1 : d){
    trun_list[[j + 1]] <- 1 : klpar[j + 1]
  }
  trunMa <- as.matrix(expand.grid(trun_list))   # all permutations and combinations for the truncation parameters
  num <- nrow(trunMa)
  aic_value <- rep(NA, num)
  for (l in 1 : num){
    result <- olsfun(Y, X, ygrid, xgrid, trunMa[l, ])  # ordinary least squares estimates 
    res <- result$Intres2
    aic_value[l] <- log(sum(res)) + 2 * (length(res)^(-1)) * sum(trunMa[l, -1])
  }
  loca_opt <- which.min(aic_value)
  aic_min <- aic_value[loca_opt]
  trun_opt <- trunMa[loca_opt, ]
  return(list(aic_min = aic_min, trun = trun_opt))
}

##### data-driven truncation parameters for predictors  #####
olsimuful <- function(data, testdat){
  ygrid <- data[, 1]$ygrid
  xgrid <- data[, 1]$xgrid
  betaTrue <- data[, 1]$beta
  p <- dim(betaTrue)[3]
  R <- dim(data)[2]
  ISE <- rep(NA, R)
  RPE <- rep(NA, R)
  for (r in 1 : R){
    Y <- data[, r]$Y
    X <- data[, r]$X
    klpar <- iniTrun(Y, X, ygrid, xgrid)
    aic_result <- AIC_ls(Y, X, ygrid, xgrid, klpar)
    trun <- aic_result$trun
    
    result <- olsfun(Y, X, ygrid, xgrid, trun)
    estbeta <- result$estbeta
    ISE[r] <- sum((estbeta - betaTrue)^2) * (xgrid[2] - xgrid[1]) * (ygrid[2] - ygrid[1])
    
    Ytest <- testdat[, r]$Y
    Xtest <- testdat[, r]$X
    fxtest <- array(NA, c(dim(Ytest), p))
    for (j in 1 : p){
      fxtest[ , , j] <- Xtest[ , , j] %*% t(estbeta[ , , j]) * (xgrid[2] - xgrid[1])
    }
    Ypred <- apply(fxtest, c(1, 2), sum)
    RPE[r] <- mean(apply((Ytest - Ypred)^2, 1, sum) * (ygrid[2] - ygrid[1])/
                     (apply(Ytest^2, 1, sum) * (ygrid[2] - ygrid[1])))
  }
  return(list(ISE = ISE, RPE = RPE))
}

