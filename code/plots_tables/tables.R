### This file reproduces all the tables in the paper
### The tables are based on the fitting results (stored in the '/run_model/results' directory), obtained by running 'run.R' 

library(xtable)

## Specify the directory to save the tables
base = getwd()
base_tables = paste0(base, '/plots_tables/tables')
## Specify the directory to read the results
base_results = paste0(base, '/run_model/results')
## Create the directory to save the tables
dir.create(base_tables, recursive = TRUE, showWarnings = FALSE)

### Functions to be used throughout this file --------
source(paste0(base, '/utils.R'))

### Table 4: LBS and BS, selection rule is Cp-edf --------

## Calculate the evaluation metrics
summary.result <- function(data, result, methods, rmse_oracle){
  # % worse than the oracle performance (best possible BS or BOSS)
  percentloss = lapply(methods, function(xx){round(100*(mean(result[[xx]][['rmse']]) / mean(rmse_oracle) - 1), 0)})
  
  # Sparsistency (extra variables)
  sparsistency.extravar <- function(sparsistency_method, extravariable_method){
    if(sum(data$beta == 0) == 0){
      # Dense model
      round(mean(sparsistency_method), 1)
    }else{
      # Sparse model
      paste(round(mean(sparsistency_method), 1), '(', round(mean(extravariable_method), 1), ')', sep='')
    }
  }
  numvar = lapply(methods, function(xx){sparsistency.extravar(result[[xx]][['sparsistency']], result[[xx]][['extravariable']])})
  
  # Relative efficiency
  mu = data$x %*% data$beta
  rmse_null = sqrt(sum(mu^2) / n) # null model
  betahat_full = ginv(data$x) %*% data$y # full OLS
  rmse_full = mean( sqrt(colSums(sweep(data$x%*%betahat_full, 1, mu,'-')^2)/n) )
  rmse_method = lapply(methods, function(xx){mean(result[[xx]][['rmse']])})
  rmse_min = min(c(unlist(rmse_method), rmse_null, rmse_full))
  efficiency = lapply(rmse_method, function(xx){round(rmse_min/xx, 2)})
  
  names(percentloss) = names(numvar) = names(efficiency) = methods
  return(list(percentloss = percentloss, 
              numvar = numvar,
              efficiency = efficiency))
}


n = c(200, 2000)
p = c(30, 180)
snr = c(7, 0.2)
nrep = 1000
type = c('Orth-Sparse-Ex1', 'Orth-Sparse-Ex2', 'Orth-Dense')
methods = c('bs', 'lbs')

count = 0
output_percentloss = output_efficiency = output_numvar = data.frame(matrix(nrow=0,ncol=6))

for(i in 1:length(n)){
  for(j in 1:length(snr)){
    for(k in 1:length(p)){
      for(l in 1:length(type)){
        # Generate the datasets
        data = gen.data.orthx(n[i], p[k], snr[j], type[l])
        mu = data$x %*% data$beta
        # Performance of LBS and BS, selection rule is Cp-edf
        result_cp = eval.metrics.bs.lbs.cp.orthx(data)
        # Best possible BS
        betahat_bs = lapply(1:nrep, function(rep){bs.orthx(data$x, data$y[,rep])})
        rmse_bs_allrep = lapply(betahat_bs, function(xx){sqrt(colSums(sweep(data$x%*%xx,1,mu,'-')^2)/n[i])})
        betahat_bs_best = do.call(cbind, Map(function(xx, yy){yy[,which.min(xx)]}, rmse_bs_allrep, betahat_bs))
        result_bs_best = rmse.sparsistency.extravariable(betahat_bs_best, data$x, data$beta)
        rmse_oracle = result_bs_best$rmse
        
        metrics = summary.result(data, result_cp, methods, rmse_oracle)
        
        
      }
    }
  }
}
  
  
