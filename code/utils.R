### This file has all the required functions to generate datasets, fit the models, and to evaluate the models

###### Simulate datasets --------
## Note that we are using R-3.6.1. The results may not be reproducible if running pre-3.6.0 versions. 
## This is because R has changed the default random number generators. 
library(MASS)
# Orthogonal X
gen.data.orthx <- function(n, p, snr, type, nrep=1000, seed=66, print.r2=FALSE, center.y=TRUE){
  set.seed(seed)
  
  # Design matrix X, columns are orthnormal with zero mean
  x = matrix(seq(0,n-1),ncol=p,nrow=n,byrow=FALSE)
  x[,1:(p/2)] = sin((2*pi/n) * (x[,1:(p/2)] %*% diag(seq(1,p/2))))
  x[,(p/2+1):p] = cos((2*pi/n) * (x[,(p/2+1):p] %*% diag(seq(1,p/2))))
  colnorm = sqrt(colSums(x^2))
  x = scale(x, center=FALSE, scale=colnorm)
  # True beta
  if(type == 'Orth-Sparse-Ex1'){
    beta = c(rep(1,6), rep(0,p-6))
  }else if(type == 'Orth-Sparse-Ex2'){
    beta = c(c(1,-1,5,-5,10,-10), rep(0,p-6))
  }else if(type == 'Orth-Dense'){
    kappa = 10
    beta = (-1)^seq(1,p) * exp(-seq(1,p)/kappa)
  }else if(type == 'Null'){
    beta = rep(0, p)
  }
  # mu and sigma
  mu = x%*%beta
  if(type == 'Null'){
    sigma = 1
  }else{
    sigma = as.numeric(sqrt(var(mu) / snr))
  }
  # 1000 replications of response y
  epsilon = matrix(rnorm(n*nrep,mean=0,sd=sigma),nrow=n,ncol=nrep)
  if(center.y){
    epsilon = scale(epsilon, center=TRUE, scale=FALSE)
  }
  y = sweep(epsilon, 1, mu, '+')
  
  # Oracle R squares (regression y upon true predictors)
  if(print.r2){
    r.squares = 0
    for(rep in 1:nrep){
      r.squares = r.squares + summary(lm(y[,rep]~x[,beta!=0]-1))$r.squared
    }
    r.squares <- r.squares/nrep
    print(paste('Average R squares of regressing upon true predictors is ', as.character(round(r.squares,2)),sep=''))
  }
  
  return(list(x=x, y=y, beta=beta, sigma=sigma))
}
# test = gen.data.orthx(n=200, p=14, snr=1.5, type='Orth-Sparse-Ex1')

# General X
gen.data.generalx <- function(n, p, rho, snr, type, nrep=1000, seed=66, print.r2=FALSE, center.y=TRUE){
  set.seed(seed)
  
  # Covariance matrix and true beta
  if(type == 'Sparse-Ex1'){
    a = matrix(rep(1:p,each=p),ncol=p,byrow=TRUE)
    b = matrix(rep(1:p,each=p),nrow=p)
    covmatrix = rho^abs(a-b)
    beta = rep(0,p)
    beta[round(seq(1,p,length.out = 6))] = 1
  }else if(type == 'Sparse-Ex2'){
    covmatrix = diag(p)
    for(i in 1:3){
      covmatrix[2*i-1,2*i] = covmatrix[2*i,2*i-1] = rho
    }
    beta = c(rep(c(1,-1),3), rep(0,p-6))
  }else if(type == 'Sparse-Ex3'){
    covmatrix = diag(p)
    for(i in 1:6){
      covmatrix[i,i+6] = covmatrix[i+6,i] = rho
    }
    beta = c(rep(1,6), rep(0,p-6))
  }else if(type == 'Sparse-Ex4'){
    covmatrix = diag(p)
    for(i in 1:3){
      covmatrix[2*i-1,2*i] = covmatrix[2*i,2*i-1] = rho
    }
    beta = c(c(1,-1,5,-5,10,-10), rep(0,p-6))
  }else if(type == 'Dense'){
    a = matrix(rep(1:p,each=p),ncol=p,byrow=TRUE)
    b = matrix(rep(1:p,each=p),nrow=p)
    covmatrix = rho^abs(a-b)
    kappa = 10
    beta = (-1)^seq(1,p) * exp(-seq(1,p)/kappa)
  }else if(type == 'Null'){
    # same correlation structure as Sparse-Ex3
    covmatrix = diag(p)
    for(i in 1:6){
      covmatrix[i,i+6] = covmatrix[i+6,i] = rho
    }
    beta = rep(0, p)
  }
  
  # Design matrix X, columns are orthnormal with zero mean
  x = mvrnorm(n, mu=rep(0,p), Sigma=covmatrix)
  x = scale(x, center=TRUE, scale=FALSE)
  colnorm = sqrt(colSums(x^2))
  x = scale(x, center=FALSE, scale=colnorm)
  
  # mu and sigma
  mu = x%*%beta
  if(type == 'Null'){
    sigma = 1
  }else{
    sigma = as.numeric(sqrt(t(beta/colnorm)%*%covmatrix%*%(beta/colnorm) / snr))
  }
  # 1000 replications of response y
  epsilon = matrix(rnorm(n*nrep,mean=0,sd=sigma),nrow=n,ncol=nrep)
  if(center.y){
    epsilon = scale(epsilon, center=TRUE, scale=FALSE)
  }
  y = sweep(epsilon, 1, mu, '+')
  
  # Oracle R squares (regression y upon true predictors)
  if(print.r2){
    r.squares = 0
    for(rep in 1:nrep){
      r.squares = r.squares + summary(lm(y[,rep]~x[,beta!=0]-1))$r.squared
    }
    r.squares <- r.squares/nrep
    print(paste('Average R squares of regressing upon true predictors is ', as.character(round(r.squares,2)),sep=''))
  }
  
  return(list(x=x, y=y, beta=beta, sigma=sigma))
}
# test = gen.data.generalx(n=200, p=14, rho=0.5, snr=0.2, type='Sparse-Ex1')

###### Fit the models ---------
library(Matrix)
library(glmnet)
library(sparsenet)
library(gamlr)
library(relaxo)
library(leaps)
library(BOSSreg)

## BS on orthogonal X
bs.orthx <- function(x, y){
  p = dim(x)[2]
  z = t(x) %*% y
  tmp = order(-z^2)
  row_i = rep(tmp, times=seq(p,1))
  col_j = unlist(lapply(2:(p+1), function(xx){seq(xx,p+1)}))
  coef_bs_all = sparseMatrix(row_i, col_j, x=z[row_i], dims=c(p, p+1))
  return(coef_bs_all)
}

## BS on general X, 'leaps' algorithm
bs.generalx <- function(x, y, intercept=TRUE){
  p = dim(x)[2]
  colnames(x) = paste0(seq(1,p))
  # fit to the full data
  outs = regsubsets(x=x, y=y, nbest=1, nvmax=p, intercept=intercept, method='exhaustive')
  if(intercept){
    val = unlist(coef(outs, id=seq(1,p)))
    names(val)[names(val) == "(Intercept)"] = '0'
    row_i = as.numeric(names(val)) + 1
    col_j = rep(2:(p+1), times=2:(p+1))
    coef_bs_all = sparseMatrix(row_i, col_j, x=val, dims=c(p+1, p+1))
    coef_bs_all[1,1] = mean(y)
  }else{
    val = unlist(coef(outs, id=seq(1,p)))
    row_i = as.numeric(names(val))
    col_j = rep(2:(p+1), times=1:p)
    coef_bs_all = sparseMatrix(row_i, col_j, x=val, dims=c(p, p+1))
  }
  return(coef_bs_all)
}

## Cross-validation for BS
cv.bs <- function(x, y, n.folds=10, orthx, intercept.generalx=TRUE){
  n = dim(x)[1]
  p = dim(x)[2]

  fold.index = sample(rep(1:n.folds,length.out=n)) # randomly assign a fold to each observation
  cv_tmp = matrix(NA,nrow=n.folds,ncol=p+1)

  for(fold in 1:n.folds){
    # Split the training and testing sets
    test.index = which(fold.index==fold)
    x.test = x[test.index, , drop=FALSE]
    y.test = y[test.index]
    x.train = x[-test.index, , drop=FALSE]
    y.train = y[-test.index]
    # Fit the model on training set
    if(orthx){
      coef_bs = bs.orthx(x.train, y.train)
    }else{
      coef_bs = bs.generalx(x.train, y.train, intercept=intercept.generalx)
    }
    # Evaluate on the testing set
    if(orthx | !intercept.generalx){
      cv_tmp[fold,] = Matrix::colMeans( sweep(x.test%*%coef_bs, 1, y.test, '-')^2 )
    }else{
      cv_tmp[fold,] = Matrix::colMeans( sweep(cbind(rep(1,nrow(x.test)),x.test)%*%coef_bs, 1, y.test, '-')^2 )
    }
  }
  cv = Matrix::colMeans(cv_tmp)
  # Fit on the full sample
  if(orthx){
    coef_bs = bs.orthx(x, y)
  }else{
    coef_bs = bs.generalx(x, y, intercept=intercept.generalx)
  }
  return(list(betahat=coef_bs, i.min=which.min(cv)))
}

## Calculate edf based on simulations
## Input: muhat (nrep by n by p+1), y (n by nrep), sigma (a number)
calc.edf <- function(muhat, y, sigma){
  edf = c()
  for(k in 1:dim(muhat[[1]])[2]){
    tmp = do.call(rbind, lapply(muhat, function(xx){xx[,k]}))
    edf = c(edf,
            sum(unlist(Map(function(xx,yy){cov(xx,yy)}, split(tmp, rep(1:ncol(tmp), each=nrow(tmp))), split(y, rep(1:nrow(y), ncol(y)))))) / sigma^2)
  }
  edf
}

## Calculate bdf based on simulations
bdf.bs.orthx <- function(x, y, seed=66){
  set.seed(seed)
  nrep = dim(y)[2]
  n = dim(x)[1]
  p = dim(x)[2]
  # Bootstrap samples
  nrep_bootstrap = 100
  betahat_multireg = t(x) %*% y
  muhat = x %*% betahat_multireg
  resid = y - muhat
  hatmatrix_diag = diag(x %*% t(x))
  resid_modified = resid / sqrt(1-hatmatrix_diag)
  tmp_function <- function(xx){
    matrix(sample(x=xx, size=n*nrep_bootstrap, replace=TRUE), nrow=n, ncol=nrep_bootstrap)
  }
  resid_bootstrap = lapply(split(resid_modified, rep(1:nrep, each=n)), tmp_function)
  y_bootstrap = Map(function(xx,j){scale(sweep(xx, 1, muhat[,j], '+' ), center=TRUE, scale=FALSE)}, resid_bootstrap, 1:nrep)
  sigma_bootstrap = sqrt(colSums(resid^2) / (n-p))
  # calculate bdf
  bdf = list()
  for(rep in 1:nrep){
    betahat = lapply(split(y_bootstrap[[rep]], rep(1:nrep_bootstrap, each=n)), function(yy){bs.orthx(x, yy)})
    bdf[[rep]] = calc.edf(lapply(betahat, function(xx){x %*% xx}), y_bootstrap[[rep]], sigma_bootstrap[rep])
  }
  return(bdf)
}

## Index of the minimum value in a two-dimensional matrix
pick.best <- function(val){
  ind_minloss = which(val == min(val), arr.ind=TRUE)
  if(nrow(ind_minloss)>1){
    ind_minloss = ind_minloss[1,]
  }
  as.numeric(ind_minloss)
}

## Cross-validation for all methods, and all replications of y
run.cv <- function(x, y, seed=66, orthx){
  n = dim(x)[1]
  p = dim(x)[2]
  nrep = dim(y)[2]

  # Declare variables to save the results
  if(orthx){
    allmethods = c('bs', 'lasso', 'sparsenet', 'gamlr', 'relaxlasso')
  }else{
    allmethods = c('bs', 'boss', 'fs', 'lasso', 'sparsenet', 'gamlr', 'relaxlasso')
  }
  i.min.cv = betahat = replicate(length(allmethods), list(), simplify=FALSE)
  names(i.min.cv) = names(betahat) = allmethods
  i.min.gamlr.aicc = ic_hdf_boss = list()

  set.seed(seed)
  for(rep in 1:nrep){
    # if(rep %% 100 == 0){print(rep)}
    print(rep)
    # BS
    if(p <= 30 | orthx){
      bs_cv = cv.bs(x, y[,rep], orthx = orthx, intercept.generalx=FALSE)
      i.min.cv$bs[[rep]] = bs_cv$i.min
      betahat$bs[[rep]] = bs_cv$betahat
    }
    if(!orthx){
      # BOSS and FS
      boss_cv = cv.boss(x, y[,rep], intercept=FALSE)
      i.min.cv$boss[[rep]] = boss_cv$i.min.boss
      i.min.cv$fs[[rep]] = boss_cv$i.min.fs
      betahat$boss[[rep]] = boss_cv$boss$beta_boss
      betahat$fs[[rep]] = boss_cv$boss$beta_fs
      ic_hdf_boss[[rep]] = boss_cv$boss$IC_boss
    }
    # LASSO
    lasso_cv = cv.glmnet(x, y[,rep], intercept=FALSE)
    i.min.cv$lasso[[rep]] = which.min(lasso_cv$cvm)
    betahat$lasso[[rep]] = coef(lasso_cv$glmnet.fit)[-1,]
    # SparseNet
    sparsenet_cv = cv.sparsenet(x, y[,rep])
    i.min.cv$sparsenet[[rep]] = rev(sparsenet_cv$which.min)
    betahat$sparsenet[[rep]] = lapply(coef(sparsenet_cv$sparsenet.fit), function(xx){xx[-1,]})
    # Gamma LASSO
    gamma_seq = c(0,1,10)
    cvm = aicc_gamlr = matrix(NA, nrow=length(gamma_seq), ncol=100)
    coef_tmp = list()
    for(j in 1:length(gamma_seq)){
      gamlr_cv_gamma = cv.gamlr(x, y[,rep], gamma=gamma_seq[j], nlambda=100, nfold=10)
      cvm[j,] = gamlr_cv_gamma$cvm
      coef_tmp[[j]] = gamlr_cv_gamma$gamlr$beta
      aicc_gamlr[j,] = calc.ic(x%*%gamlr_cv_gamma$gamlr$beta, y[,rep], ic='aicc', df=gamlr_cv_gamma$gamlr$df)
    }
    i.min.cv$gamlr[[rep]] = pick.best(cvm)
    betahat$gamlr[[rep]] = coef_tmp
    i.min.gamlr.aicc[[rep]] = pick.best(aicc_gamlr)
    # Relaxed LASSO
    # Sometimes the function 'cvrelaxo' throws error messages, re-try for at most 3 times
    attempt = 1
    while( attempt <= 3 ){
      relaxlasso_cv = try(
        cvrelaxo(x, y[,rep], K=10, keep.data=FALSE), silent = TRUE
      )
      if(class(relaxlasso_cv) == 'try-error'){
        attempt = attempt + 1
      }else{
        attempt = 4
      }
    }
    if(class(relaxlasso_cv) == 'try-error'){
      stop('relaxed LASSO CV wrong')
    }
    relaxlasso_fullsample = relaxo(x, y[,rep], keep.data=FALSE, phi = seq(0, 1, length = 10))
    tmp = which(relaxlasso_fullsample$lambda == relaxlasso_cv$lambda & relaxlasso_fullsample$phi == relaxlasso_cv$phi)
    if(length(tmp) != 1){
      stop('relaxed LASSO CV wrong')
    }
    i.min.cv$relaxlasso[[rep]] = which(relaxlasso_fullsample$lambda == relaxlasso_cv$lambda & relaxlasso_fullsample$phi == relaxlasso_cv$phi)
    betahat$relaxlasso[[rep]] = Matrix(t(relaxlasso_fullsample$beta), sparse = TRUE)
  }

  return(list(i.min.cv = i.min.cv,
              i.min.gamlr.aicc = i.min.gamlr.aicc,
              ic_hdf_boss = ic_hdf_boss,
              betahat = betahat))
}

###### Evaluate the models --------

# betahat is p by nrep, selected coefficient vector
rmse.sparsistency.extravariable <- function(betahat, x, beta){
  if(is.null(betahat)){
    return(NULL)
  }
  n = dim(x)[1]
  rmse = sqrt(colSums(sweep(x%*%betahat, 1, x%*%beta,'-')^2)/n)
  sparsistency = colSums(betahat[beta!=0, , drop=FALSE] != 0)
  extravariable = colSums(betahat[beta==0, , drop=FALSE] != 0)
  return(list(rmse=rmse, sparsistency=sparsistency, extravariable=extravariable))
}

## Evaluate the selected subsets for all methods
eval.metrics <- function(x, y, beta, sigma, result.cv, bdf.bs, orthx){
  betahat = result.cv$betahat
  i.min.cv = result.cv$i.min.cv
  i.min.gamlr.aicc = result.cv$i.min.gamlr.aicc
  ic_hdf_boss = result.cv$ic_hdf_boss

  mu = x%*%beta
  n = dim(x)[1]
  p = dim(x)[2]
  nrep = dim(y)[2]

  allmethods = names(betahat)

  result = list()
  # Best possible BS or BOSS
  if(orthx){
    method = 'bs'
  }else{
    method = 'boss'
  }
  rmse_method_allrep = lapply(betahat[[method]], function(xx){sqrt(colSums(sweep(x%*%xx,1,mu,'-')^2)/n)})
  betahat_method_best = do.call(cbind, Map(function(xx, yy){yy[,which.min(xx)]}, rmse_method_allrep, betahat[[method]]))
  result[[method]][['best']] = rmse.sparsistency.extravariable(betahat_method_best, x, beta)

  if(orthx){
    # BS-IC-df
    df_bs = list()
    tmp = calc.edf(lapply(betahat$bs, function(xx){x %*% xx}), y, sigma)
    df_bs$edf = replicate(nrep, tmp, simplify = FALSE)
    tmp = lapply(split(y, rep(1:nrep, each=n)), function(yy){BOSSreg:::calc.hdf(x, yy)})
    df_bs$hdf = lapply(tmp, function(xx){xx$hdf})
    sigmahat = lapply(tmp, function(xx){xx$sigma})
    df_bs$ndf = replicate(nrep, 0:p, simplify = FALSE)
    df_bs$bdf = bdf.bs

    ic_bs = list()
    ic_bs = lapply(df_bs[c('hdf','ndf','bdf')], function(df){lapply(1:nrep, function(rep){BOSSreg:::calc.ic.all(betahat$bs[[rep]], x=x, y=y[,rep], df=df[[rep]], sigma=sigmahat[[rep]])})})
    ic_bs$edf = lapply(1:nrep, function(rep){BOSSreg:::calc.ic.all(betahat$bs[[rep]], x=x, y=y[,rep], df=df_bs$edf[[rep]], sigma=sigma)})

    ic_method = ic_bs
  }else{
    # BOSS-IC-hdf
    ic_boss = list()
    ic_boss$hdf = ic_hdf_boss
    ic_method = ic_boss
  }

  tmp = lapply(ic_method, function(xx){lapply(1:nrep, function(rep){lapply(xx[[rep]], function(yy){betahat[[method]][[rep]][,which.min(yy)]})})})
  df_names = names(tmp)
  ic_names = names(tmp$hdf[[1]])
  for(ic in ic_names){
    for(df in df_names){
      betahat_method_ic = do.call(cbind, lapply(tmp[[df]], function(xx){xx[[ic]]}) )
      result[[method]][['ic']][[ic]][[df]] = rmse.sparsistency.extravariable(betahat_method_ic, x, beta)
    }
  }

  # LASSO-AICc
  betahat_lasso_aicc = do.call(cbind, lapply(1:nrep, function(rep){betahat$lasso[[rep]][,which.min( calc.ic(x%*%betahat$lasso[[rep]], y[,rep], ic='aicc', df=colSums(betahat$lasso[[rep]]!=0)) )]} ))
  result$lasso$ic$aicc = rmse.sparsistency.extravariable(betahat_lasso_aicc, x, beta)

  # Gamma LASSO-AICc
  betahat_gamlr_aicc = do.call(cbind, lapply(1:nrep, function(rep){ betahat$gamlr[[rep]][[i.min.gamlr.aicc[[rep]][1]]][,i.min.gamlr.aicc[[rep]][2]] }))
  result$gamlr$ic$aicc = rmse.sparsistency.extravariable(betahat_gamlr_aicc, x, beta)

  # CV for all the methods
  for(method in allmethods){
    if(method %in% c('sparsenet', 'gamlr')){
      betahat_method_cv = do.call(cbind, lapply(1:nrep, function(rep){ betahat[[method]][[rep]][[i.min.cv[[method]][[rep]][1]]][,i.min.cv[[method]][[rep]][2]] }))
    }else if(method == 'bs' & p > 30 & !orthx){
      betahat_method_cv = NULL
    }else{
      betahat_method_cv = do.call(cbind, lapply(1:nrep, function(rep){ betahat[[method]][[rep]][,i.min.cv[[method]][[rep]]] }))
    }
    result[[method]][['cv']] = rmse.sparsistency.extravariable(betahat_method_cv, x, beta)
  }

  return(result)
}

###### Other functions to be used in reproducing plots and tables --------
## Degrees of freedom for LBS when X is orthogonal, based on definitions in Tibshirani(2015)
edf.lbs.orthx <- function(Q, sigma, mu, sqrt_2lambda){
  n = dim(Q)[1]
  p = dim(Q)[2]
  
  # if mu and sigma are not specified, use the full multiple regression
  xtmu = t(Q)%*%mu
  
  sqrt_2lambda_matrix = matrix(rep(sqrt_2lambda,each=p),nrow=p,byrow=F)
  xtmu_matrix = matrix(rep(xtmu,each=length(sqrt_2lambda)),ncol=length(sqrt_2lambda),byrow=T)
  
  a = stats::dnorm((sqrt_2lambda_matrix-xtmu_matrix) / sigma)
  b = stats::dnorm((-sqrt_2lambda_matrix-xtmu_matrix) / sigma)
  c = stats::pnorm((sqrt_2lambda_matrix-xtmu_matrix) / sigma)
  d = stats::pnorm((-sqrt_2lambda_matrix-xtmu_matrix) / sigma)
  
  size = colSums(1 - c + d)
  sdf = (sqrt_2lambda/sigma) * colSums(a + b)
  df = size + sdf
  return(df)
}

## Fit Lagrangian BS, X is orthogonal, here y is a matrix (n by nrep)
lbs.orthx <- function(x, y){
  p = dim(x)[2]
  nrep = dim(y)[2]
  # specify a fixed sequence of lambda
  z = t(x)%*%y
  nlambda = 200
  sqrt_2lambda = exp(seq(log(max(abs(z))+0.1), log(0.001*(max(abs(z))+0.1)), length.out=nlambda))
  sqrt_2lambda = matrix(rep(sqrt_2lambda, each=nrep), ncol=nrep, byrow=TRUE)
  
  # calculate the fit based on the analytical formula
  tmp_function <- function(ii){
    row_i = lapply(sqrt_2lambda[,ii], function(xx){which(abs(z[,ii])>=xx)})
    col_j = Map(function(xx, yy){rep(yy, length(xx))}, row_i, 1:length(row_i))
    val = lapply(row_i, function(xx){z[xx,ii]})
    sparseMatrix(unlist(row_i), unlist(col_j), x=unlist(val), dims=c(p, length(sqrt_2lambda[,ii])))
  }
  betahat = lapply(1:nrep, tmp_function)
  
  return(list(betahat=betahat, sqrt_2lambda=sqrt_2lambda))
}
## Evaluate the BS-Cp and LBS-Cp, 'data' is the object by running 'gen.data.orthx'
eval.metrics.bs.lbs.cp.orthx <- function(data){
  # Fit BS and LBS
  betahat = edf = betahat_cp = list()
  result_lbs = lbs.orthx(data$x, data$y)
  betahat$lbs = result_lbs$betahat
  betahat$bs = lapply(1:nrep, function(rep){bs.orthx(data$x, data$y[,rep])})
  muhat = lapply(betahat, function(xx){lapply(xx, function(yy){data$x %*% yy})})
  rss = lapply(muhat, function(xx){do.call(cbind, Map(function(yy, rep){colSums(sweep(yy, 1, data$y[,rep])^2)}, xx, 1:nrep))})
  # edf 
  edf$lbs = apply(result_lbs$sqrt_2lambda, 2, function(xx){edf.lbs.orthx(data$x, data$sigma, data$x %*% data$beta, xx)})
  edf$bs = calc.edf(muhat$bs, data$y, data$sigma)
  # Cp
  cp = Map(function(xx, yy){xx + 2*data$sigma^2*yy}, rss, edf)
  # Coefficients selected by Cp-edf
  betahat_cp = Map(function(aa, bb){do.call(cbind, Map(function(xx,yy){xx[,yy]}, aa, apply(bb, 2, which.min)))}, betahat, cp)
  # Evaluate the results
  result_cp = lapply(betahat_cp, function(xx){rmse.sparsistency.extravariable(xx, data$x, data$beta)})
  return(result_cp)
}
