
##############################################################
### GCD algorithm for the proposed method ###
GCDfun <- function(Y, X, ygrid, xgrid, trunPar){
  n <- dim(Y)[1]
  if(is.matrix(X)) X = array(X, c(dim(X), 1)) else X = X 
  p <- dim(X)[3]
  k <- trunPar[1]      # the truncation parameter for y
  lpar <- trunPar[-1]  # the truncation parameter vector for p functional predictors
  suml <- sum(lpar)
  yfdat <- fdata(Y, ygrid)
  ypc <- fdata2pc(yfdat, ncomp = k)
  ybasis <- ypc$rotation$data   # the estimated fpc basis
  esteta <- ypc$x[, 1 : k]  # the firsd k estimated fpc scores for y
  
  xbasis <- list()
  estxsi_list <- list()
  group_list <- list()
  for (j in 1 : p){
    xfdat <- fdata(X[, , j], xgrid)
    xpc <- fdata2pc(xfdat, ncomp = lpar[j])
    xbasis[[j]] <- xpc$rotation$data
    estxsi_list[[j]] <- xpc$x[, 1 : lpar[j]]
    group_list[[j]] <- rep(j, lpar[j])
  }
  estxsi <- matrix(unlist(estxsi_list), ncol = suml) # desigen matrix of fpc scores for predictors
  group <- unlist(group_list)
  cvfit <- cv.grpreg(estxsi, esteta, penalty = "grSCAD", nfolds = 5, group = group, 
                     lambda.min = 0.01)
  
  if (is.matrix(coef(cvfit))) betacoef <- coef(cvfit)[, -1] else
      betacoef <- t(as.matrix(coef(cvfit)[-1]))
  betacoef_list <- lapply(split(t(betacoef), group), matrix, ncol = k) # split a matrix to some submatrices by row
  est_beta <- array(NA, c(length(ygrid), length(xgrid), p)) # array for the estimated beta
  for (j in  1 : p){
    est_beta[, , j] <- t(ybasis) %*% t(betacoef_list[[j]]) %*% xbasis[[j]]
  }
  yfit = estxsi %*% t(betacoef) %*% ybasis
  Integral.res <- (norm.fdata(yfdat - fdata(yfit, ygrid)))^2
  return(list(estbeta = est_beta, yfit = yfit, Intres2 = Integral.res))
}

########################################################################
####### ols method for the selected model  #######
olsfun <- function(Y, X, ygrid, xgrid, trunPar){
  n <- dim(Y)[1]
  if(is.matrix(X)) X = array(X, c(dim(X), 1)) else X = X 
  p <- dim(X)[3]
  k <- trunPar[1]
  lpar <- trunPar[-1]
  suml <- sum(lpar)
  yfdat <- fdata(Y, ygrid)
  ypc <- fdata2pc(yfdat, ncomp = k)
  ybasis <- ypc$rotation$data
  
  esteta <- ypc$x[, 1 : k]
  xbasis <- list()
  estxsi_list <- list()
  group_list <- list()
  for (j in 1 : p){
    xfdat <- fdata(X[, , j], xgrid)
    xpc <- fdata2pc(xfdat, ncomp = lpar[j])
    xbasis[[j]] <- xpc$rotation$data
    estxsi_list[[j]] <- xpc$x[, 1 : lpar[j]]
    group_list[[j]] <- rep(j, lpar[j])
  }
  estxsi <- matrix(unlist(estxsi_list), ncol = suml)
  group <- unlist(group_list)
  yy <- as.vector(t(esteta))
  xx <- kronecker(estxsi, diag(rep(1, k)))
  ols <- lm(yy ~ xx - 1)
  coefMa <- matrix(ols$coefficients, nrow = k)
  
  betacoef_list <- lapply(split(t(coefMa), group), matrix, ncol = k) # split a matrix to some submatrices by row
  est_beta <- array(NA, c(length(ygrid), length(xgrid), p))  # array for the estimated beta
  for (j in  1 : p){
    est_beta[ , , j] <- t(ybasis) %*% t(betacoef_list[[j]]) %*% xbasis[[j]]
  }
  yfit = estxsi %*% t(coefMa) %*% ybasis
  Integral.res <- (norm.fdata(yfdat - fdata(yfit, ygrid)))^2
  return(list(estbeta = est_beta, yfit = yfit, Intres2 = Integral.res))
}
