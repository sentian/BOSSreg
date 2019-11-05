###### Simulate datasets --------
library(MASS)
# Orthogonal X
gen.data.orthx <- function(n, p, snr, type, nrep=1000, seed=66, print.r2=FALSE){
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
  }
  # mu and sigma
  mu = x%*%beta
  sigma = as.numeric(sqrt(var(mu) / snr))
  # 1000 replications of response y, each with zero mean
  epsilon = scale(matrix(rnorm(n*nrep,mean=0,sd=sigma),nrow=n,ncol=nrep), center=TRUE, scale=FALSE)
  y = sweep(epsilon, 1, mu, '+')

  # Oracle R squares (regression y upon true predictors)
  if(print.r2){
    r.squares = 0
    for(rep in 1:nrep){
      r.squares = r.squares + summary(lm(y[,rep]~x[,beta!=0]-1))$r.squared
    }
    r.squares <- r.squares/nrep
    print(paste('Average R squares of regressing upon true predictors is', as.character(round(r.squares,2)),sep=''))
  }

  return(list(x=x, y=y, beta=beta, sigma=sigma))
}

# n = c(200, 2000)
# p = c(14, 30, 60, 180)
# snr = c(0.2, 1.5, 7)
# names(snr) = c('lsnr', 'msnr', 'hsnr')
# type = c('Orth-Sparse-Ex1', 'Orth-Sparse-Ex2', 'Orth-Dense')
# test = gen.data.orthx(n=200, p=14, snr=1.5, type='Orth-Sparse-Ex1')

# General X
gen.data.generalx <- function(n, p, rho, snr, type, nrep=1000, seed=66, print.r2=FALSE){
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
  }

  # Design matrix X, columns are orthnormal with zero mean
  x = mvrnorm(n, mu=rep(0,p), Sigma=covmatrix)
  x = scale(x, center=TRUE, scale=FALSE)
  colnorm = sqrt(colSums(x^2))
  x = scale(x, center=FALSE, scale=colnorm)

  # mu and sigma
  mu = x%*%beta
  sigma = as.numeric(sqrt(t(beta/colnorm)%*%covmatrix%*%(beta/colnorm) / snr))
  # 1000 replications of response y, each with zero mean
  epsilon = scale(matrix(rnorm(n*nrep,mean=0,sd=sigma),nrow=n,ncol=nrep), center=TRUE, scale=FALSE)
  y = sweep(epsilon, 1, mu, '+')

  # Oracle R squares (regression y upon true predictors)
  if(print.r2){
    r.squares = 0
    for(rep in 1:nrep){
      r.squares = r.squares + summary(lm(y[,rep]~x[,beta!=0]-1))$r.squared
    }
    r.squares <- r.squares/nrep
    print(paste('Average R squares of regressing upon true predictors is', as.character(round(r.squares,2)),sep=''))
  }

  return(list(x=x, y=y, beta=beta, sigma=sigma))
}

# n = c(200, 2000)
# p = c(14, 30, 60, 180)
# rho = c(0, 0.5, 0.9)
# snr = c(0.2, 1.5, 7)
# names(snr) = c('lsnr', 'msnr', 'hsnr')
# type = c('Sparse-Ex1', 'Sparse-Ex2', 'Sparse-Ex4', 'Sparse-Ex5', 'Dense')
#
# test = gen.data.generalx(n=200, p=14, rho=0.5, snr=0.2, type='Sparse-Ex1')

###### Run the models ---------
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
## Input: betahat (nrep by n by p+1), y, sigma
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

## Cross-validation for all methods, and all replications of y
## Output: for each method, the index of subset that gives minimum CV error, and the entire coefficient matrix
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

  set.seed(seed)
  for(rep in 1:nrep){
    if(rep %% 100 == 0){print(rep)}
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
    }
    # LASSO
    lasso_cv = cv.glmnet(x, y[,rep], intercept=FALSE)
    i.min.cv$lasso[[rep]] = which.min(lasso_cv$cvm)
    betahat$lasso[[rep]] = coef(lasso_cv$glmnet.fit)[-1,]
    # SparseNet
    sparsenet_cv = cv.sparsenet(x, y[,rep])
    i.min.cv$sparsenet[[rep]] = sparsenet_cv$which.min
    betahat$sparsenet[[rep]] = lapply(coef(sparsenet_cv$sparsenet.fit), function(xx){xx[-1,]})
    # Gamma LASSO
    gamma_seq = c(0,1,10)
    cvm = matrix(NA, nrow=length(gamma_seq), ncol=100)
    coef_tmp = list()
    for(j in 1:length(gamma_seq)){
      gamlr_cv_gamma = cv.gamlr(x, y[,rep], gamma=gamma_seq[j], nlambda=100, nfold=10)
      cvm[j,] = gamlr_cv_gamma$cvm
      coef_tmp[[j]] = gamlr_cv_gamma$gamlr$beta
    }
    ind_minloss = which(cvm == min(cvm), arr.ind=TRUE)
    if(nrow(ind_minloss)>1){
      ind_minloss = ind_minloss[1,]
    }
    i.min.cv$gamlr[[rep]] = as.numeric(ind_minloss)
    betahat$gamlr[[rep]] = coef_tmp
    # Relaxed LASSO
    relaxlasso_cv = cvrelaxo(x, y[,rep], K=10, keep.data=FALSE)
    relaxlasso_fullsample = relaxo(x, y[,rep], keep.data=FALSE, phi = seq(0, 1, length = 10))
    tmp = which(relaxlasso_fullsample$lambda == relaxlasso_cv$lambda & relaxlasso_fullsample$phi == relaxlasso_cv$phi)
    if(length(tmp) != 1){
      stop('relaxed LASSO CV wrong')
    }
    i.min.cv$relaxlasso[[rep]] = which(relaxlasso_fullsample$lambda == relaxlasso_cv$lambda & relaxlasso_fullsample$phi == relaxlasso_cv$phi)
    betahat$relaxlasso[[rep]] = Matrix(t(relaxlasso_fullsample$beta), sparse = TRUE)
  }

  return(list(i.min.cv = i.min.cv,
              betahat = betahat))
}

# betahat is p by nrep, selected coefficient vector
rmse.sparsistency.extravariable <- function(betahat, x, beta){
  n = dim(x)[1]
  rmse = sqrt(colSums(sweep(x%*%betahat, 1, x%*%beta,'-')^2)/n)
  sparsistency = colSums(betahat[beta!=0, , drop=FALSE] != 0)
  extravariable = colSums(betahat[beta==0, , drop=FALSE] != 0)
  return(list(rmse=rmse, sparsistency=sparsistency, extravariable=extravariable))
}

eval.metrics <- function(x, y, beta, sigma, betahat, i.min.cv, orthx){
  allmethods = names(betahat)
  mu = x%*%beta
  n = dim(x)[1]
  p = dim(x)[2]
  nrep = dim(y)[2]

  # result$extravariable$bestsub$ic$cov$edf
  result = list()
  # BS: IC with various types of df, CV and oracle (best possible)
  rmse_bs_allrep = lapply(betahat$bs, function(xx){sqrt(colSums(sweep(x%*%xx,1,mu,'-')^2)/n)})
  betahat_bs_best = do.call(cbind, Map(function(xx, yy){yy[,which.min(xx)]}, rmse_bs_allrep, betahat$bs))
  result$bs$best = rmse.sparsistency.extravariable(betahat_bs_best, x, beta)

  df_bs = list()
  tmp = calc.edf(lapply(betahat$bs, function(xx){x %*% xx}), y, sigma)
  df_bs$edf = replicate(nrep, tmp, simplify = FALSE)
  tmp = lapply(split(y, rep(1:nrep, each=n)), function(yy){BOSSreg:::calc.hdf(x, yy)})
  df_bs$hdf = lapply(tmp, function(xx){xx$hdf})
  sigmahat = lapply(tmp, function(xx){xx$sigma})
  df_bs$ndf = replicate(nrep, 0:p, simplify = FALSE)
  df_bs$bdf = bdf.bs.orthx(x, y)

  ic_bs = list()
  ic_bs = lapply(df_bs[c('hdf','ndf','bdf')], function(df){lapply(1:nrep, function(rep){BOSSreg:::calc.ic.all(betahat$bs[[rep]], x=x, y=y[,rep], df=df[[rep]], sigma=sigmahat[[rep]])})})
  ic_bs$edf = lapply(1:nrep, function(rep){BOSSreg:::calc.ic.all(betahat$bs[[rep]], x=x, y=y[,rep], df=df_bs$edf[[rep]], sigma=sigma)})



}

# test run
x = unname(as.matrix(read.table('/home/st1864/model_selection/data/orthx/x.txt')))
y = unname(as.matrix(read.table('/home/st1864/model_selection/data/orthx/y_hsnr.txt')))
sigma = unname(as.numeric(read.table('/home/st1864/model_selection/data/orthx/sd_hsnr.txt')))
beta = as.numeric(unname(as.matrix(read.table('/home/st1864/model_selection/data/orthx/beta.txt'))))

# 15 minutes
result = run.cv(x, y, orthx=TRUE)
saveRDS(result$betahat, paste0('/scratch/st1864/boss/betahat/',type,'_n',n,'_p',p,'_',names(snr),'.rds')) # save in the tmp folder
saveRDS(result$i.min.cv, paste0('/scratch/st1864/boss/tmp/cv/',type,'_n',n,'_p',p,'_',names(snr),'.rds')) # save in the tmp folder
ptm = proc.time()
bdf_bs = bdf.bs.orthx(x, y)
proc.time() - ptm
saveRDS(bdf_bs, paste0('/scratch/st1864/boss/tmp/bdf/',type,'_n',n,'_p',p,'_',names(snr),'.rds')) # save in the tmp folder

betahat = readRDS(paste0('/scratch/st1864/boss/betahat/',type,'_n',n,'_p',p,'_',names(snr),'.rds')) # save in the tmp folder

calc.edf(lapply(betahat$bs, function(xx){x %*% xx}), y, sigma)

indir = paste('true_model/orthogonal/sparse/ex1/p0_6','/n_',200,'/p_',14, sep='')
df = readRDS(paste(base, '/code/as_gram/results/',indir,'/hsnr/df_bs.rds',sep=""))

n = c(200)
p = c(14)
snr = c(7)
names(snr) = c('hsnr')
type = c('Orth-Sparse-Ex1')
