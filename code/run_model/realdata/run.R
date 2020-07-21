### This file fits the models on real data examples
### Final results are stored in the directory '../results/realdata/'
library(pbapply)
library(ISLR)
library(MASS)
library(BOSSreg)
library(glmnet)
library(sparsenet)
library(leaps)

### Source functions --------
source(paste0(getwd(), '/utils.R'))

### Processing the data --------
dataset = list()
# Boston Housing data
tmp = Boston
tmp = na.omit(tmp)
tmp$chas = as.factor(tmp$chas)
dataset$boston$x = data.matrix(tmp[,!names(tmp) %in% 'medv'])
dataset$boston$y = tmp$medv

# MLB hitters salary
tmp = Hitters
tmp = na.omit(tmp)
tmp[,c('League', 'Division', 'NewLeague')] = 
  lapply(tmp[,c('League', 'Division', 'NewLeague')], as.factor)
dataset$hitters$x = data.matrix(tmp[,!(names(tmp) %in% c('Salary'))])
dataset$hitters$y = tmp$Salary

# College data
tmp = College
tmp$Private = as.factor(tmp$Private)
dataset$college$x = data.matrix(tmp[,!(names(tmp) %in% c('Outstate'))])
dataset$college$y = tmp$Outstate

# Auto data
tmp = Auto
dataset$auto$x = data.matrix(tmp[,!(names(tmp) %in% c('mpg','name','origin'))])
dataset$auto$y = tmp$mpg

### Functions to be used --------
## Root mean squared error
rmse <- function(y_hat, y){
  sqrt(sum( (y_hat - y)^2 / length(y)) )
}
## Leave-one-out errors, number of variables chosen, and running time
rdresult <- function(x, y, nrep, seed){
  p = dim(x)[2]
  
  allmethods = c('lasso', 'sparsenet', 'boss', 'fs', 'bs')
  error = numvar = time = replicate(length(allmethods), rep(NA,nrep), simplify=F)
  names(error) = names(numvar) = names(time) = allmethods
  
  set.seed(seed)
  for(i in 1:nrep){
    index = 1:nrow(x)
    index = index[-i]
    
    x.train = x[index, , drop=FALSE]
    y.train = y[index]
    x.test = x[-index, , drop=FALSE]
    x.test.withint = cbind(rep(1,nrow(x.test)), x.test)
    y.test = y[-index]
    
    # BS
    ptm = proc.time()
    bs_cv_model = cv.bs(x.train, y.train)
    time_tmp = proc.time() - ptm
    bs_pred = as.numeric( x.test.withint %*%  bs_cv_model$betahat[,bs_cv_model$i.min])
    error$bs[i] = rmse(bs_pred, y.test)
    numvar$bs[i] = sum(bs_cv_model$betahat[,bs_cv_model$i.min]!=0)
    time$bs[i] = time_tmp[3]
    
    # BOSS
    ptm = proc.time()
    boss_model = boss(x.train, y.train, intercept = TRUE)
    time_tmp = proc.time() - ptm
    boss_pred = as.numeric( predict(boss_model, newx=x.test) )
    error$boss[i] = rmse(boss_pred, y.test)
    numvar$boss[i] = sum(coef(boss_model)!=0)
    time$boss[i] = time_tmp[3]
    
    # FS
    ptm = proc.time()
    boss_cv_model = cv.boss(x.train, y.train)
    time_tmp = proc.time() - ptm
    fs_pred = as.numeric( predict(boss_cv_model, newx=x.test, method='fs') )
    error$fs[i] = rmse(fs_pred, y.test)
    numvar$fs[i] = sum(coef(boss_cv_model, method='fs')!=0)
    time$fs[i] = time_tmp[3]
    
    # LASSO
    ptm = proc.time()
    lasso_model = glmnet(x.train, y.train, intercept=TRUE)
    lasso_aicc = as.numeric(calc.ic(predict(lasso_model, newx=x.train), y.train, 
                                    ic='aicc', df=lasso_model$df+1))
    lasso_pred = predict(lasso_model, newx=x.test, s=lasso_model$lambda[which.min(lasso_aicc)])
    time_tmp = proc.time() - ptm
    error$lasso[i] = rmse(lasso_pred, y.test)
    numvar$lasso[i] = sum(coef(lasso_model, s=lasso_model$lambda[which.min(lasso_aicc)])!=0)
    time$lasso[i] = time_tmp[3]
    
    # SparseNet
    ptm = proc.time()
    sparsenet_cv_model = cv.sparsenet(x.train, y.train)
    time_tmp = proc.time() - ptm
    sparsenet_pred = predict(sparsenet_cv_model, newx=x.test, which='parms.min')
    error$sparsenet[i] = rmse(sparsenet_pred, y.test)
    numvar$sparsenet[i] = sum(coef(sparsenet_cv_model, which='parms.min')!=0)
    time$sparsenet[i] = time_tmp[3]
  }
  return(list(error=error, numvar=numvar, time=time))
}

### Fit the models --------
## It takes around 20 minutes to run on a single core of a local machine (i7 CPU and 16 GB RAM)
result = pblapply(dataset, function(xx){rdresult(xx$x, xx$y, nrow(xx$x), seed=66)})
saveRDS(result, paste0(getwd(), '/run_model/realdata/results/boston_hitters_college_auto.rds'))
