
##### initial values of the truncation parameters for Y and X  #####
iniTrun <- function(Y, X, ygrid, xgrid){
  p <- dim(X)[3]
  yfdat <- fdata(Y, ygrid)
  ypc <- fdata2pc(yfdat)
  y.sd <- ypc$d
  y.var <- y.sd^2/sum(y.sd^2)  # explained variance for the all the yfpc scores
  y.cumvar <- cumsum(y.var)
  k <- which(y.cumvar >= 0.99)[1]
  
  trunPar <- rep(1, p)
  for (j in 1 : p){
    xfdat <- fdata(X[, , j], xgrid)
    xpc <- fdata2pc(xfdat)
    x.sd <- xpc$d
    x.var <- x.sd^2/sum(x.sd^2)
    x.cumvar <- cumsum(x.var)
    trunPar[j] <- which(x.cumvar >= 0.95)[1]
  }
  trun <- c(k, trunPar)
  return(trun)
}

#### AIC criterion for refitting and selecting truncation parameters of X

AIC <- function(Y, X, ygrid, xgrid, klpar){
  k <- klpar[1]
  if(is.matrix(X)) X = array(X, c(dim(X), 1)) else X = X 
  d <- dim(X)[3]
  
  lcomset <- matrix(rep(2 : max(klpar[-1]), d), ncol = d)  #common parameter sets for each X
  num0 <- nrow(lcomset)
  klcomset <- cbind(rep(k, num0), lcomset)
  aic_value0 <- rep(NA, num0)
  for (l in 1 : num0){
    result0 <- olsfun(Y, X, ygrid, xgrid, klcomset[l, ])  # ordinary least squares estimates
    res0 <- result0$Intres2
    aic_value0[l] <- log(sum(res0)) + 2 * (length(res0)^(-1)) * sum(klcomset[l, -1])
  }
  loca_opt0 <- which.min(aic_value0)
  opl_com <- klcomset[loca_opt0, ][2]
  
  trun_list <- list()
  trun_list[[1]] <- k
  lcom_min <- ifelse((opl_com - 2) > 1, opl_com - 2, 1)
  for (j in 1 : d){
    trun_list[[j + 1]] <- c(lcom_min : (opl_com + 2))
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

###### CVAIC method for choosing the finial parameters ######
CVAIC <- function(Y, X, ygrid, xgrid){
  p <- dim(X)[3]
  klPar <- iniTrun(Y, X, ygrid, xgrid)
  iniEst <- GCDfun(Y, X, ygrid, xgrid, klPar)
  
  estbeta <- iniEst$estbeta
  x_select <- rep(NA, p)          # 1 represents nonzro predictor, 0 corresponds to zero predictor
  for (j in 1 : p){
    x_select[j] <- 1 - (sum(abs(estbeta[, , j])) == 0)
  }
  nonzero_loca <- which(x_select == 1)    # locations of identified nonzero predictors
  X_nonzero <- X[ , , nonzero_loca]
  klpar_nonzero <- c(klPar[1], klPar[nonzero_loca + 1])
  
  aic_result <- AIC(Y, X_nonzero, ygrid, xgrid, klpar_nonzero)
  trun_opt <- aic_result$trun
  return(list(optrun = trun_opt, initun = klPar, x_select = x_select, nzloca = nonzero_loca))
}
