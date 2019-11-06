### This file fits the models and evaluate them
### This file was called by run_orthx.sh and run_generalx.sh, in order to be run on an HPC server parallelly for all the configurations of true models

###### Parameters for each model configuration, to be used as input files for fitting models on HPC server --------
para.forhpc <- function(orthx){
  count = 0
  n = c(200, 2000)
  p = c(14, 30, 60, 180)
  snr = c(0.2, 1.5, 7)
  if(orthx){
    type = c('Orth-Sparse-Ex1', 'Orth-Sparse-Ex2', 'Orth-Dense')
    for(i in 1:length(type)){
      for(j in 1:length(n)){
        for(k in 1:length(p)){
          for(l in 1:length(snr)){
            count = count + 1
            write.table(c(type[i], n[j], p[k], snr[l]),
                        file = paste(base, "/run_model/para_forhpc/orthx_",count,".txt",sep=''), sep = '', col.names = FALSE, row.names = FALSE)

          }
        }
      }
    }
  }else{
    type = c('Sparse-Ex1', 'Sparse-Ex2', 'Sparse-Ex3', 'Sparse-Ex4', 'Dense')
    rho = c(0, 0.5, 0.9)
    for(i in 1:length(type)){
      for(j in 1:length(n)){
        for(k in 1:length(p)){
          for(l in 1:length(snr)){
            for(m in 1:length(rho)){
              count = count + 1
              write.table(c(type[i], n[j], p[k], snr[l], rho[m]),
                          file = paste(base, "/run_model/para_forhpc/generalx_",count,".txt",sep=''), sep = '', col.names = FALSE, row.names = FALSE)
            }
          }
        }
      }
    }

  }


}
# base = '/Users/sentian/Dropbox/Sen/Research/Model_selection/rpackage/BOSSreg/code'
# para.forhpc(TRUE)
# para.forhpc(FALSE)

###### Fit the models on HPC server --------
base = '/scratch/st1864/boss'
# Read the parameters
args = (commandArgs(TRUE))
orthx = as.logical( as.character(args[1]) )
num = as.character(Sys.getenv('SLURM_ARRAY_TASK_ID'))

if(orthx){
  para = c(t(read.table(file = paste(base, "/para_forhpc/orthx_", num, ".txt",sep=''), sep = '', header=FALSE)))
  # create folders to store the results (both temporary and final results)
  dir.create(base, 'tmp/bdf_bs', recursive = TRUE)
  dir.create(base, 'tmp/result_cv/orthx', recursive = TRUE)
  dir.create(base, 'results/orthx', recursive = TRUE)
}else{
  para = c(t(read.table(file = paste(base, "/para_forhpc/generalx_", num, ".txt",sep=''), sep = '', header=FALSE)))
  rho = as.numeric(para[5])
  # create folders to store the results (both temporary and final results)
  dir.create(base, 'tmp/result_cv/generalx', recursive = TRUE)
  dir.create(base, 'results/generalx', recursive = TRUE)
}
type = para[1]
n = as.numeric(para[2])
p = as.numeric(para[3])
snr = as.numeric(para[4])

snr_all = c(0.2, 1.5, 7)
names(snr_all) = c('lsnr', 'msnr', 'hsnr')
snr_name = names(snr_all)[which(snr_all == snr)]

# Source all the required libraries and functions
source('/home/st1864/boss/code/utils.R')

# Generate simulated datasets
if(orthx){
  data = gen.data.orthx(n=n, p=p, snr=snr, type=type)
}else{
  data = gen.data.generalx(n=n, p=p, rho=rho, snr=snr, type=type)
}

# Bootstrapped df for BS, only when X is orthogonal
if(orthx){
  print('bdf')
  bdf_bs = bdf.bs.orthx(data$x, data$y)
  saveRDS(bdf_bs, paste0(base, '/tmp/bdf_bs/',type,'_n',n,'_p',p,'_',snr_name,'.rds'))
}else{
  bdf_bs = NULL
}

# Cross-valications
print('cv')
ptm = proc.time()
result_cv = run.cv(data$x, data$y, orthx=orthx)
proc.time() - ptm
if(orthx){
  saveRDS(result_cv, paste0(base, '/tmp/result_cv/orthx/',type,'_n',n,'_p',p,'_',snr_name,'.rds'))
}else{
  saveRDS(result_cv, paste0(base, '/tmp/result_cv/generalx/',type,'_n',n,'_p',p,'_',snr_name,'_rho',gsub("[.]","",as.character(rho)),'.rds'))
}

# Evaulate the selected subsets for all methods
print('evaluation')
result = eval.metrics(data$x, data$y, data$beta, data$sigma, result_cv, bdf_bs, orthx)
if(orthx){
  saveRDS(result, paste0(base, '/results/orthx/',type,'_n',n,'_p',p,'_',snr_name,'.rds'))
}else{
  saveRDS(result, paste0(base, '/results/generalx/',type,'_n',n,'_p',p,'_',snr_name,'_rho',gsub("[.]","",as.character(rho)),'.rds'))
}











n.folds = 10
n.rep = 1
intercept = FALSE

n = dim(x)[1]
p = dim(x)[2]
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
    maxstep = trunc(min(n - n/n.folds, p))
    # matrix to store the CV error
    cv_rep_boss = cv_rep_fs = matrix(NA, nrow=n.rep, ncol=maxstep+1)

    for(replication in 1:n.rep){
      fold.index = sample(rep(1:n.folds, length.out=n)) # randomly assign a fold to each observation
      cv_tmp_boss = cv_tmp_fs = matrix(NA, nrow=n.folds, ncol=maxstep+1)

      for(fold in 1:n.folds){
        # split the training and testing sets
        test.index = which(fold.index==fold)
        x.test = x[test.index, , drop=FALSE]
        y.test = y[test.index, rep]
        x.train = x[-test.index, , drop=FALSE]
        y.train = y[-test.index, rep]
        boss_result <- boss(x.train, y.train, intercept, hdf.ic.boss=FALSE)
        beta_fs= boss_result$beta_fs
        beta_boss = boss_result$beta_boss
        # if intercept
        if(intercept){
          x.test = cbind(rep(1,nrow(x.test)), x.test)
        }
        sweep(x.test%*%beta_fs, 1, y.test, '-')
        cv_tmp_fs[fold, ] = Matrix::colMeans(sweep(x.test%*%beta_fs, 1, y.test, '-')^2)
        cv_tmp_boss[fold, ] = Matrix::colMeans(sweep(x.test%*%beta_boss, 1, y.test, '-')^2)

      }
      cv_rep_fs[replication, ] = Matrix::colMeans(cv_tmp_fs)
      cv_rep_boss[replication, ] = Matrix::colMeans(cv_tmp_boss)
    }


    # i.min.cv$boss[[rep]] = boss_cv$i.min.boss
    # i.min.cv$fs[[rep]] = boss_cv$i.min.fs
    # betahat$boss[[rep]] = boss_cv$boss$beta_boss
    # betahat$fs[[rep]] = boss_cv$boss$beta_fs
    # ic_hdf_boss[[rep]] = boss_cv$boss$IC_boss
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
  relaxlasso_cv = cvrelaxo(x, y[,rep], K=10, keep.data=FALSE)
  relaxlasso_fullsample = relaxo(x, y[,rep], keep.data=FALSE, phi = seq(0, 1, length = 10))
  tmp = which(relaxlasso_fullsample$lambda == relaxlasso_cv$lambda & relaxlasso_fullsample$phi == relaxlasso_cv$phi)
  if(length(tmp) != 1){
    stop('relaxed LASSO CV wrong')
  }
  i.min.cv$relaxlasso[[rep]] = which(relaxlasso_fullsample$lambda == relaxlasso_cv$lambda & relaxlasso_fullsample$phi == relaxlasso_cv$phi)
  betahat$relaxlasso[[rep]] = Matrix(t(relaxlasso_fullsample$beta), sparse = TRUE)
}


saveRDS(list(x=x, y=y), '/Users/sentian/Downloads/data.rds')
data = readRDS('/Volumes/HDD/Dropbox/data.rds')
boss_result = boss(data$x, data$y, intercept=FALSE, hdf.ic.boss=FALSE)
x = data$x
y = data$y
test = fs(x, y, intercept = FALSE)
test2 = boss(x, y, intercept = FALSE)

