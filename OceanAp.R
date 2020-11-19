
library(fda.usc)
library(grpreg)
library(fields) # heatmaps  image.plot
library(abind)  # combine matrices or list to a array
source("estimation.R")

####################################################################
########### tuning parameters selection  #############
## initial values of the truncation parameters for Y and X##
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
    trunPar[j] <- which(x.cumvar >= 0.99)[1]
  }
  trun <- c(k, trunPar)
  return(trun)
}

#### AIC criterion for refitting and selecting truncation parameters of X
AIC <- function(Y, X, ygrid, xgrid, klpar){
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

##################################################################################


########### ocean dataset ################
Ocean <- read.table("Ocean1999-2018.txt", header = T, sep = ",", stringsAsFactors = FALSE)
id <- unique(Ocean$crn)
temper <- tapply(Ocean$temp, Ocean$crn, function(x)x)
Potemper <- tapply(Ocean$theta, Ocean$crn, function(x)x)
salinity <- tapply(Ocean$sal, Ocean$crn, function(x)x)
oxygen <- tapply(Ocean$oxy, Ocean$crn, function(x)x)
density <- tapply(Ocean$sigma, Ocean$crn, function(x)x) 
chlo <- tapply(Ocean$fluor, Ocean$crn, function(x)x) 
Nitrate <- tapply(Ocean$nit, Ocean$crn, function(x)x) 

salinityMa <- matrix(unlist(salinity), byrow = T, nrow = length(id))
chlomMa <- matrix(unlist(chlo), byrow = T, nrow = length(id))
id.na <- unique(c(which(salinityMa == -9, arr.ind = T)[, 1],
                  which(chlomMa == -9, arr.ind = T)[, 1]))

temperMa <- matrix(unlist(temper[-id.na]), byrow = T, nrow = length(id) - length(id.na))
salMa <- matrix(unlist(salinity[-id.na]), byrow = T, nrow = length(id) - length(id.na))
oxyMa <- matrix(unlist(oxygen[-id.na]), byrow = T, nrow = length(id) - length(id.na))
densityMa <- matrix(unlist(density[-id.na]), byrow = T, nrow = length(id) - length(id.na))
PotemMa <- matrix(unlist(Potemper[-id.na]), byrow = T, nrow = length(id) - length(id.na))
NitMa <- matrix(unlist(Nitrate[-id.na]), byrow = T, nrow = length(id) - length(id.na))
chloMa <- matrix(unlist(chlo[-id.na]), byrow = T, nrow = length(id) - length(id.na))

temfdat <- fdata(temperMa, argvals = seq(2, 200, by = 2))
salfdat <- fdata(salMa, argvals = seq(2, 200, by = 2))
oxyfdat <- fdata(oxyMa, argvals = seq(2, 200, by = 2))
densityfdat <- fdata(densityMa, argvals = seq(2, 200, by = 2))
Potemfdat <- fdata(PotemMa, argvals = seq(2, 200, by = 2))
Nitfdat <- fdata(NitMa, argvals = seq(2, 200, by = 2))
chlofdat <- fdata(chloMa, argvals = seq(2, 200, by = 2))


########### plots of the five variables  ############

pdf(file = "Temperature.pdf")
plot(temfdat, lty = 1, ylab = "Temperature", xlab = "Depth below the sea surface", 
     main = "Temperature", mgp = c(2, 0.5, 0))
dev.off()

pdf(file = "Salinity.pdf")
plot(salfdat, lty = 1, ylab = "Salinity", xlab = "Depth below the sea surface", 
     main = "Salinity", mgp = c(2, 0.5, 0))
dev.off()

pdf(file = "Oxygen.pdf")
plot(oxyfdat, lty = 1, ylab = "Oxygen", xlab = "Depth below the sea surface", 
     main = "Oxygen", mgp = c(2, 0.5, 0))
dev.off()

pdf(file = "Density.pdf")
plot(densityfdat, lty = 1, ylab = "Potential Density", xlab = "Depth below the sea surface", 
     main = "Potential Density", mgp = c(2, 0.5, 0))
dev.off()

pdf(file = "Chloropigment.pdf")
plot(chlofdat, lty = 1, ylab = "Chloropigment", xlab = "Depth below the sea surface", 
     main = "Chloropigment", mgp = c(2, 0.5, 0))
dev.off()


############### centralize the data ##################

y <- scale(temperMa, scale = F)
x1 <- scale(oxyMa, scale = F)
x2 <- scale(densityMa, scale = F)
x3 <- scale(salMa, scale = F)
x4 <- scale(chloMa, scale = F)
x <- abind(x1, x2, x3, x4, along = 3)
t <- seq(2, 200, by = 2)
s <- seq(2, 200, by = 2)
ygrid <- t / 200; xgrid <- s / 200


############## estimating coefficients by using gSCAD method #############

trun <- CVAIC(y, x, ygrid, xgrid)$optrun  # selected parameter values 5 1 5 5 1
est <- GCDfun(y, x, ygrid, xgrid, trun)
estbeta <- est$estbeta

pdf(file = "beta1.pdf")
image.plot(t, s, estbeta[, , 1], main = expression(hat(beta)[1](t, s)), 
           xlab = expression(t), ylab = expression(s), mgp = c(2, 0.5, 0))
dev.off()

pdf(file = "beta2.pdf")
image.plot(t, s, estbeta[, , 2], main = expression(hat(beta)[2](t, s)), 
           xlab = expression(t), ylab = expression(s), mgp = c(2, 0.5, 0))
dev.off()

pdf(file = "beta3.pdf")
image.plot(t, s, estbeta[, , 3], main = expression(hat(beta)[3](t, s)), 
           xlab = expression(t), ylab = expression(s), mgp = c(2, 0.5, 0))
dev.off()

pdf(file = "beta4.pdf")
image.plot(t, s, estbeta[, , 4], main = expression(hat(beta)[4](t, s)), 
           xlab = expression(t), ylab = expression(s), mgp = c(2, 0.5, 0))
dev.off()

# the proposed method select x2 and x3 as two relevant preditors

yfit <- est$yfit
res.sum <- apply((y - yfit)^2, 2, sum)
var.sum <- apply(y^2, 2, sum)
fR2 <- 1 - sum(res.sum/var.sum) * (ygrid[2] - ygrid[1]) # functional R square for the selected model 
#fR2 = 0.989682

### including x1 ### 
x.add1 <- x[, , c(1, 2, 3)]
trun.ini <- iniTrun(y, x.add1, ygrid, xgrid)
klpar <- AIC(y, x.add1, ygrid, xgrid, trun.ini)$trun # selected paramter values 5 1 5 5
result.add1 <- olsfun(y, x.add1, ygrid, xgrid, klpar)

yfit.add1 <- result.add1$yfit
res.sum <- apply((y - yfit.add1)^2, 2, sum)
fR2.add1 <- 1 - sum(res.sum/var.sum) * (ygrid[2] - ygrid[1]) 
# fR2 = 0.9896897

### including x1 and x4 ### 
trun.ini <- iniTrun(y, x, ygrid, xgrid)
klpar <- AIC(y, x, ygrid, xgrid, trun.ini)$trun # selected paramter values 5 1 5 5 1
result.add14 <- olsfun(y, x, ygrid, xgrid, klpar)

yfit.add14 <- result.add14$yfit
res.sum <- apply((y - yfit.add14)^2, 2, sum)
fR2.add14 <- 1 - sum(res.sum/var.sum) * (ygrid[2] - ygrid[1]) 
# fR2 = 0.989716

#calculate functional R square for the model with only non-selected predictors#
x.nonselect <- x[, , c(1, 4)]
trun.ini <- iniTrun(y, x.nonselect, ygrid, xgrid)
klpar <- AIC(y, x.nonselect, ygrid, xgrid, trun.ini)$trun # selected paramter values 5 6 8
result.nonse <- olsfun(y, x.nonselect, ygrid, xgrid, klpar)

yfit.nons <- result.nonse$yfit
res.sum <- apply((y - yfit.nons)^2, 2, sum)
fR2.nons <- 1 - sum(res.sum/var.sum) * (ygrid[2] - ygrid[1]) 
# fR2.nons = 0.5990233


############### prediction over 200 repeats #####################
predfun <- function(y, x, ygrid, xgrid){
  n = nrow(y)
  ntrain = floor(0.7 * n)
  ntest = n - ntrain
  ind.n = sample(c(1 : n), ntrain, replace = FALSE)
  ytrain <- y[ind.n, ]
  ytest <- y[-ind.n, ]
  xtrain <- x[ind.n, , ] 
  xtest <- x[-ind.n, , ]
  
  trun <- CVAIC(y, x1, x2, x3, x4, ygrid, xgrid)$optrun
  est <- GCDfun(ytrain, xtrain, ygrid, xgrid, trun)
  estbeta <- est$estbeta
  p <- dim(xtest)[3]
  fxtest <- array(NA, c(dim(ytest), p))
  for (j in 1 : p){
    fxtest[ , , j] <- xtest[ , , j] %*% t(estbeta[ , , j]) * (xgrid[2] - xgrid[1])
  }
  ypred <- apply(fxtest, c(1, 2), sum)
  RPE <- mean(apply((ytest - ypred)^2, 1, sum) * (ygrid[2] - ygrid[1])/
                   (apply(ytest^2, 1, sum) * (ygrid[2] - ygrid[1])))
  
  ###### marginal model containing only x1 and x4
  x.nonselect <- xtrain[, , c(1, 4)]
  trun.ini <- iniTrun(ytrain, x.nonselect, ygrid, xgrid)
  klpar <- AIC(ytrain, x.nonselect, ygrid, xgrid, trun.ini)$trun
  result.nonse <- olsfun(ytrain, x.nonselect, ygrid, xgrid, klpar)
  estbeta.nons <- result.nonse$estbeta
  
  ypred.nons <- (xtest[, , 1] %*% t(estbeta.nons[, , 1]) + 
                 xtest[, , 4] %*% t(estbeta.nons[, , 2])) * (xgrid[2] - xgrid[1])
  RPE.nons <- mean(apply((ytest - ypred.nons)^2, 1, sum) * (ygrid[2] - ygrid[1])/
                (apply(ytest^2, 1, sum) * (ygrid[2] - ygrid[1])))
  
  #### full model including all the four predictors
  trun.ini2 <- iniTrun(ytrain, xtrain, ygrid, xgrid)
  trun.full <- AIC(ytrain, xtrain, ygrid, xgrid, trun.ini2)$trun 
  olsest <- olsfun(ytrain, xtrain, ygrid, xgrid, trun.full)
  estbeta.full <- olsest$estbeta
  ypred.full <- (xtest[, , 1] %*% t(estbeta.full[, , 1]) + 
                 xtest[, , 2] %*% t(estbeta.full[, , 2]) + 
                 xtest[, , 3] %*% t(estbeta.full[, , 3]) + 
                 xtest[, , 4] %*% t(estbeta.full[, , 4])) * (xgrid[2] - xgrid[1])
  RPE.full <- mean(apply((ypred.full - ytest)^2, 1, sum) * (ygrid[2] - ygrid[1])/
                      (apply(ytest^2, 1, sum) * (ygrid[2] - ygrid[1])))
  
  PE <- c(RPE, RPE.nons, RPE.full)
  names(PE) <- c("RPE", "RPE.nons", "RPE.full")
  return(PE)
}

pred <- replicate(200, predfun(y, x, ygrid, xgrid))


apply(pred, 1, mean)
#RPE   RPE.nons   RPE.full 
#0.01599392 0.61703567 0.01606222
apply(pred, 1, sd)
#RPE    RPE.nons    RPE.full 
#0.002342174 0.124579814 0.002382175 
save.image("OceanPre.RData")



