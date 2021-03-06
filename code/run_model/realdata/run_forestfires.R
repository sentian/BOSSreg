### This file fits the models on the forest fires dataset, available from https://archive.ics.uci.edu/ml/datasets/Forest+Fires
### This file was called by run_forestfires.sh, in order to be run on an HPC server parallelly
### Final results are stored in the directory '../results/realdata'

## Since the dimension p=55, we fit BS via MIO algorithm that was proposed by Bertsimas et. al. (2015)
## An R package for the implementation is 'bestsubset', available at https://github.com/ryantibs/best-subset

library(bestsubset)
library(BOSSreg)
library(glmnet)
library(sparsenet)

### Functions ---------
rmse <- function(y_hat, y){
  sqrt(sum( (y_hat - y)^2 / length(y)) )
}
## Fit BS via MIO
bs.mio <- function(x, y, k_stop){
  bs_model = bs(x, y, k=0:k_stop, time.limit = 180)
  # add intercept to beta
  bs_coef = rbind(bs_model$by - bs_model$bx %*% bs_model$beta, bs_model$beta)
  return(bs_coef)
}
## cross-validation for BS via MIO
cv.bs.mio <- function(x, y, n.folds=10, k_stop){
  n = dim(x)[1]
  p = dim(x)[2]

  fold.index = sample(rep(1:n.folds,length.out=n)) # randomly assign a fold to each observation
  cv_tmp = matrix(NA,nrow=n.folds,ncol=k_stop+1)
  
  for(fold in 1:n.folds){
    # split the training and testing sets
    test.index = which(fold.index==fold)
    x.test = x[test.index,, drop=FALSE]
    y.test = y[test.index]
    x.train = x[-test.index,, drop=FALSE]
    y.train = y[-test.index]
    # fit n training set
    betahat_train = bs.mio(x.train, y.train, k_stop) 
    # evaluate on the testing set
    cv_tmp[fold,] = colMeans(sweep(cbind(1, x.test)%*%betahat_train, y.test, 1, '-')^2)
  }
  cv = colMeans(cv_tmp)
  betahat = bs.mio(x, y, k_stop)
  return(list(betahat=betahat, i.min=which.min(cv)))
}

# i indicates which partition of LOO to run
rdresult <- function(x, y, i, seed){
  p = dim(x)[2]
  
  allmethods = c('lasso', 'sparsenet', 'boss', 'fs', 'bs')
  error = numvar = time = replicate(length(allmethods), NA, simplify=F)
  names(error) = names(numvar) = names(time) = allmethods
  
  set.seed(seed)
  
  index = 1:nrow(x)
  index = index[-i]
  
  x.train = x[index, , drop=FALSE]
  y.train = y[index]
  x.test = x[-index, , drop=FALSE]
  x.test.withint = cbind(1, x.test)
  y.test = y[-index]
  
  # BS
  ptm = proc.time()
  bs_cv_model = cv.bs.mio(x.train, y.train, 10)
  time_tmp = proc.time() - ptm
  bs_pred = as.numeric( x.test.withint %*%  bs_cv_model$betahat[,bs_cv_model$i.min])
  error$bs = rmse(bs_pred, y.test)
  numvar$bs = sum(bs_cv_model$betahat[,bs_cv_model$i.min]!=0)
  time$bs = time_tmp[3]
  
  # BOSS
  ptm = proc.time()
  boss_model = boss(x.train, y.train, intercept = TRUE)
  time_tmp = proc.time() - ptm
  boss_pred = as.numeric( predict(boss_model, newx=x.test) )
  error$boss = rmse(boss_pred, y.test)
  numvar$boss = sum(coef(boss_model)!=0)
  time$boss = time_tmp[3]
  
  # FS
  ptm = proc.time()
  boss_cv_model = cv.boss(x.train, y.train)
  time_tmp = proc.time() - ptm
  fs_pred = as.numeric( predict(boss_cv_model, newx=x.test, method='fs') )
  error$fs = rmse(fs_pred, y.test)
  numvar$fs = sum(coef(boss_cv_model, method='fs')!=0)
  time$fs = time_tmp[3]
  
  # LASSO
  ptm = proc.time()
  lasso_model = glmnet(x.train, y.train, intercept=TRUE)
  lasso_aicc = as.numeric(calc.ic(predict(lasso_model, newx=x.train), y.train, 
                                  ic='aicc', df=lasso_model$df+1))
  lasso_pred = predict(lasso_model, newx=x.test, s=lasso_model$lambda[which.min(lasso_aicc)])
  time_tmp = proc.time() - ptm
  error$lasso = rmse(lasso_pred, y.test)
  numvar$lasso = sum(coef(lasso_model, s=lasso_model$lambda[which.min(lasso_aicc)])!=0)
  time$lasso = time_tmp[3]
  
  # SparseNet
  ptm = proc.time()
  sparsenet_cv_model = cv.sparsenet(x.train, y.train)
  time_tmp = proc.time() - ptm
  sparsenet_pred = predict(sparsenet_cv_model, newx=x.test, which='parms.min')
  error$sparsenet = rmse(sparsenet_pred, y.test)
  numvar$sparsenet = sum(coef(sparsenet_cv_model, which='parms.min')!=0)
  time$sparsenet = time_tmp[3]
  
  return(list(error=error, numvar=numvar, time=time))
}


## Fit the model --------
base = '/scratch/st1864/boss' # the basis of the directory to store all the outputs
## Environment parameters: which partition of the leave-one-out to run
i = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

dir.create(paste0(base, '/tmp/forestfires'), showWarnings = FALSE, recursive = TRUE)

## Read and process the data
data = read.table(paste0(getwd(), '/forestfires.csv'), sep = ",", quote = "\"", header = T)
x = data[,!names(data) %in% c('area', 'month', 'day')]
x = model.matrix( ~.^2, x)[,-1] # add interaction terms
y = data$area

start = 1 + 50*(num-1)
end = 50*num
if(num == 10){
  end = nrow(x)[1]
}
result = lapply(start:end, function(i){rd.result(x, y, i, 66)})
saveRDS(result, paste0(base, '/tmp/forestfires/', i, '.rds'))

### Combine all the results --------
# torun = FALSE
# if(torun){
#   dir.create(paste0(base, '/results/realdata'), showWarnings = FALSE, recursive = TRUE)
#   
#   base = '/scratch/st1864/boss' # the basis of the directory to store all the outputs
#   tmp = lapply(1:517, function(i){readRDS(paste0(base, '/tmp/forestfires/', i, '.rds'))})
#   result = lapply(1:3, do.call(c, lapply(tmp, '[[', i)))
#   names(result) = c('error', 'numvar', 'time')
#   saveRDS(result, paste0(base, '/results/realdata/forestfires.rds'))
# }

