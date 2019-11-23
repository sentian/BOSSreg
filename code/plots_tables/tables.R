### This file reproduces all the tables in the paper
### The tables are based on the fitting results (stored in the '/run_model/results' directory), obtained by running 'run.R' 
### Resulting tables (in Tex format) are stored in '../paper/tables' directory
library(xtable)

## Specify the directory to save the tables
# The default directory on my machine is 'code/', change it to its parent directory 'BOSSreg/'
setwd("..")
base = getwd()
base_tables = paste0(base, '/paper/tables')
## Specify the directory to read the results
base_results = paste0(base, '/code/run_model/simulation/results')
## Create the directory to save the tables
dir.create(paste0(base_tables, '/supplement/boss'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(base_tables, '/supplement/bs_ic'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(base_tables, '/supplement/bs_regu'), recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(base_tables, '/supplement/boss'), recursive = TRUE, showWarnings = FALSE)

### Functions to be used throughout this file --------
source(paste0(base, '/code/utils.R'))
## Calculate the evaluation metrics: % worse than oracle method, relative efficiency, sparsistency (extra variables)
summary.result <- function(data, result, rmse.oracle, type.table='maintext.method'){
  n = dim(data$x)[1]
  # % worse than the oracle performance (best possible BS or BOSS)
  percentloss = unlist(lapply(result$rmse, function(xx){round(100*(xx / rmse.oracle - 1), 0)}))
  
  # Sparsistency (extra variables)
  sparsistency.extravar <- function(sparsistency_method, extravariable_method){
    if(sum(data$beta == 0) == 0){
      # Dense model
      as.character(round(sparsistency_method, 1))
    }else{
      # Sparse model
      paste(round(sparsistency_method, 1), '(', round(extravariable_method, 1), ')', sep='')
    }
  }
  numvar = unlist(Map(sparsistency.extravar, result$sparsistency, result$extravariable))
  
  # Relative efficiency
  mu = data$x %*% data$beta
  rmse_null = sqrt(sum(mu^2) / n) # null model
  betahat_full = ginv(data$x) %*% data$y # full OLS
  rmse_full = mean( sqrt(colSums(sweep(data$x%*%betahat_full, 1, mu,'-')^2)/n) )
  rmse_method = unlist(result$rmse)
  rmse_min = min(c(rmse_method, rmse_null, rmse_full))
  efficiency = round(rmse_min/rmse_method, 2)
  
  output = list(percentloss = percentloss, efficiency = efficiency, numvar = numvar)
  
  # Different styles of output for different types of tables
  if(type.table == 'maintext.method'){
    return(output)
  }else if(type.table == 'maintext.bsrule'){
    tmp_function <- function(xx){
      # We don't have BS results for p>30
      if(length(xx) == 12){
        xx = c(xx[1:12], '-')
      }
      c(xx[1], paste(xx[2:4],collapse='/'), xx[5], paste(xx[6:8],collapse='/'), xx[9], paste(xx[10:12],collapse='/'), xx[13])
    } 
    return(lapply(output, tmp_function))
  }else if(type.table == 'maintext.bossregu'){
    tmp_function <- function(xx){
      # We don't have BS results for p>30
      if(length(xx) == 5){
        xx = c(xx[1:2], '-', xx[3:5])
      }
      c(paste(xx[1:2],collapse='/'), xx[3:6])
    } 
    return(lapply(output, tmp_function))
  }else if(type.table == 'supplement.bsrule'){
    tmp_function <- function(xx){
      # We don't have BS results for p>30
      if(length(xx) == 16){
        xx = c(xx[1:16], '-')
      }
      c(xx[1], paste(xx[2:4],collapse='/'), xx[5], paste(xx[6:8],collapse='/'), xx[9], paste(xx[10:12],collapse='/'), 
        xx[13], paste(xx[14:16],collapse='/'), xx[17])
    } 
    return(lapply(output, tmp_function))
  }else if(type.table == 'supplement.bsregu'){
    tmp_function <- function(xx){
      c(xx[1], paste(xx[2:3],collapse='/'), paste(xx[4:5],collapse='/'), xx[6:7])
    }
    return(lapply(output, tmp_function))
  }else if(type.table == 'supplement.bossregu'){
    tmp_function <- function(xx){
      # We don't have BS results for p>30
      if(length(xx) == 10){
        xx = c(xx[1:3], '-', xx[4:10])
      }
      c(paste(xx[1:3],collapse='/'), xx[4:5], paste(xx[6:7],collapse='/'), paste(xx[8:9],collapse='/'), xx[10:11])
    } 
    return(lapply(output, tmp_function))
  }

}
## Evaluation metrics for all scenarios considered in a table
df.summary.result <- function(n, p, snr, type, rho=NULL, methods, type.table){
  count = 0
  output = list()
  if(!is.null(rho)){
    # X is general
    for(i in 1:length(n)){
      for(j in 1:length(snr)){
        for(m in 1:length(rho)){
          for(k in 1:length(p)){
            metrics = list()
            count = count + 1
            for(l in 1:length(type)){
              # Read data
              data = gen.data.generalx(n[i], p[k], rho[m], snr[j], type[l])
              # Read the simulation results
              filename = paste0('/generalx/', type[l], '_n', n[i], '_p', p[k], '_', names(snr)[j], '_rho', gsub("[.]","",as.character(rho[m])))
              result = readRDS(paste0(base_results, filename, '.rds'))
              if(!'bs' %in% names(result)){
                result$bs = list(NULL)
              }
              result_target = lapply(methods, function(xx){result[[xx]]})
              result_target = Filter(Negate(is.null), result_target) # Remove the CV results for BS when p>30
              # Best possible BOSS
              rmse_oracle = mean(result$boss$best$rmse)
              # The evaluation metrics for the specified methods
              result_tosummary = list()
              for(which.metric in c('rmse', 'sparsistency', 'extravariable')){
                result_tosummary[[which.metric]] = lapply(result_target, function(xx){mean(xx[[which.metric]])})
              }
              # Calculate the three metrics to be presented in the table
              metrics[[l]] = summary.result(data, result_tosummary, rmse_oracle, type.table)
            }
            output[[count]] = lapply(1:3, function(ii){do.call(c, lapply(metrics, '[[', ii))})
          }
        }
      }
    }
  }else{
    # X is orthogonal
    for(i in 1:length(n)){
      for(j in 1:length(snr)){
        for(k in 1:length(p)){
          metrics = list()
          count = count + 1
          for(l in 1:length(type)){
            # Read data
            data = gen.data.orthx(n[i], p[k], snr[j], type[l])
            # Read the simulation results
            filename = paste0('/orthx/', type[l], '_n', n[i], '_p', p[k], '_', names(snr)[j])
            result = readRDS(paste0(base_results, filename, '.rds'))
            if(!'bs' %in% names(result)){
              result$bs = list(NULL)
            }
            result_target = lapply(methods, function(xx){result[[xx]]})
            result_target = Filter(Negate(is.null), result_target) # Remove the CV results for BS when p>30
            # Best possible BS
            rmse_oracle = mean(result$bs$best$rmse)
            # Average evaluation metrics for the specified methods
            result_tosummary = list()
            for(which.metric in c('rmse', 'sparsistency', 'extravariable')){
              result_tosummary[[which.metric]] = lapply(result_target, function(xx){mean(xx[[which.metric]])})
            }
            # Calculate the three metrics to be presented in the table
            metrics[[l]] = summary.result(data, result_tosummary, rmse_oracle, type.table)
          }
          output[[count]] = lapply(1:3, function(ii){do.call(c, lapply(metrics, '[[', ii))})
        }
      }
    }
  }
  output = do.call(rbind, lapply(1:3, function(ii){do.call(rbind, lapply(output, '[[', ii))}))
  return(output) 
}

### Table 1: Selection rules for BS, Orth-Sparse-Ex1 is the true model --------
table.bs.selectrule <- function(type, title, label, scale.box){
  # Parameters
  n = c(200, 2000)
  p = c(30, 180)
  snr = c(7, 0.2)
  names(snr) = c('hsnr', 'lsnr')
  nrep = 1000

  output = df.summary.result(n, p, snr, type, 
                             methods = do.call(c, list(lapply(c('edf', 'ndf', 'hdf', 'bdf'), function(xx){c('bs', 'ic', 'cp', xx)}),
                                                       lapply(c('edf', 'ndf', 'hdf', 'bdf'), function(xx){c('bs', 'ic', 'aicc', xx)}),
                                                       lapply(c('edf', 'ndf', 'hdf', 'bdf'), function(xx){c('bs', 'ic', 'bic', xx)}),
                                                       list(c('bs', 'cv')))),
                             type.table = 'maintext.bsrule')
  
  output = cbind(rep(c("\\multirow{4}[4]{*}{n=200}", "", "\\cmidrule{2-10}", "", "\\midrule \\multirow{4}[4]{*}{n=2000}", "", "\\cmidrule{2-10}", ""), 3),
                 rep(c("\\multirow{2}[2]{*}{hsnr}", "", "\\multirow{2}[2]{*}{lsnr}", ""), 6),
                 rep(c("p=30", "p=180"), 12),
                 output)
  command = c(paste("\\toprule \n",
                    "\\multicolumn{1}{|r}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{2}{c|}{C$_p$} & \\multicolumn{2}{c|}{AICc} & \\multicolumn{2}{c|}{BIC} & \\multirow{2}[4]{*}{CV} \\\\\n",
                    "\\cmidrule{4-9}\\multicolumn{1}{|r}{} & \\multicolumn{1}{c}{} &       & edf   & ndf/hdf/bdf & edf   & ndf/hdf/bdf & edf   & ndf/hdf/bdf &     \\\\\n",
                    "\\midrule \n",
                    "\\multicolumn{1}{|r}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{7}{c|}{\\% worse than the best possible BS} \\\\\n",
                    "\\midrule \n"),
              paste("\\midrule \n",
                    "\\multicolumn{1}{|r}{} & \\multicolumn{1}{r}{} &       & \\multicolumn{7}{c|}{Relative efficiency} \\\\\n",
                    "\\midrule \n"),
              paste("\\midrule \n",
                    "\\multicolumn{1}{|r}{} & \\multicolumn{1}{r}{} &       & \\multicolumn{7}{c|}{Sparsistency (number of extra variables)} \\\\\n",
                    "\\midrule \n"),
              paste("\\bottomrule \n"
              )
  )
  print(xtable(output,
               align = "l|c|c|c|cc|cc|cc|c|",  # align and put a vertical line (first "l" again represents column of row numbers)
               label = label,
               caption = title),
        #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
        scalebox = scale.box,
        caption.placement = "top",
        include.rownames = FALSE, 
        include.colnames = FALSE, 
        hline.after = NULL, 
        floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
        sanitize.text.function = force, # Important to treat content of first column as latex function
        add.to.row = list(pos = list(-1,
                                     8,
                                     16,
                                     nrow(output)),
                          command = command
        ),
        file = paste0(base_tables, "/bs_selectrule_", type, ".tex"), compress = FALSE
  )
}
table.bs.selectrule(type = 'Orth-Sparse-Ex1', 
                    title = "The performance of AICc-hdf. The true model setup is Orth-Sparse-Ex1. 
                    The columns involving `edf' refer to infeasible selection rules since edf is estimated as if the true model is known, 
                    while other columns correspond to feasible rules.",
                    label = "tab:ic_df_orthx_sparseex1",
                    scale.box = 0.65)

### Table 2: Selection rules for BS, Orth-Dense is the true model --------
table.bs.selectrule(type = 'Orth-Dense', 
                    title = "The performance of AICc-hdf. The true model setup is Orth-Dense. Details of the columns can be referred to 
                    the caption in Table \\ref{tab:ic_df_orthx_sparseex1}",
                    label = "tab:ic_df_orthx_dense",
                    scale.box = 0.7)

### Table 3: Selection rules for BS, Orth-Dense is the true model --------
table.bs.regu <- function(title, scale.box){
  # Parameters
  n = c(200, 2000)
  p = c(30, 180)
  snr = c(7, 0.2)
  names(snr) = c('hsnr', 'lsnr')
  nrep = 1000
  type = c('Orth-Sparse-Ex1', 'Orth-Sparse-Ex2', 'Orth-Dense')
  
  output = df.summary.result(n, p, snr, type, 
                             methods = list(c('bs', 'ic', 'aicc', 'hdf'), c('lasso', 'ic', 'aicc'), c('sparsenet', 'cv')),
                             type.table = 'maintext.method')
  
  output = cbind(rep(c("\\multirow{4}[4]{*}{n=200}", "", "\\cmidrule{2-12}", "", "\\midrule \\multirow{4}[4]{*}{n=2000}", "", "\\cmidrule{2-12}", ""), 3),
                 rep(c("\\multirow{2}[2]{*}{hsnr}", "", "\\multirow{2}[2]{*}{lsnr}", ""), 6),
                 rep(c("p=30", "p=180"), 12),
                 output)
  command = c(paste("\\toprule \n",
                    "\\multicolumn{1}{|r}{} & \\multicolumn{1}{r}{} &       & \\multicolumn{3}{c|}{Orth-Sparse-Ex1}   & \\multicolumn{3}{c|}{Orth-Sparse-Ex2}             & \\multicolumn{3}{c|}{Orth-Dense} \\\\\n",
                    "\\cmidrule{4-12}\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & BS    & LASSO & SparseNet & BS    & LASSO & SparseNet & BS    & LASSO & SparseNet  \\\\\n",
                    "\\midrule \n",
                    "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{9}{c|}{\\% worse than the best possible BS}  \\\\\n",
                    "\\midrule \n"),
              paste("\\midrule \n",
                    "\\multicolumn{1}{|r}{} & \\multicolumn{1}{r}{} &       & \\multicolumn{9}{c|}{Relative efficiency} \\\\\n",
                    "\\midrule \n"),
              paste("\\midrule \n",
                    "\\multicolumn{1}{|r}{} & \\multicolumn{1}{r}{} &       & \\multicolumn{9}{c|}{Sparsistency (number of extra variables)} \\\\\n",
                    "\\midrule \n"),
              paste("\\bottomrule \n"
              )
  )
  print(xtable(output,
               align = "l|c|c|c|ccc|ccc|ccc|",  # align and put a vertical line (first "l" again represents column of row numbers)
               label = 'tab:bs_regularization',
               caption = title),
        #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
        scalebox = scale.box,
        caption.placement = "top",
        include.rownames = FALSE, 
        include.colnames = FALSE, 
        hline.after = NULL, 
        floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
        sanitize.text.function = force, # Important to treat content of first column as latex function
        add.to.row = list(pos = list(-1,
                                     8,
                                     16,
                                     nrow(output)),
                          command = command
        ),
        file = paste0(base_tables, "/bs_regu.tex"), compress = FALSE
  )
}
table.bs.regu(title = "The performance of BS compared to regularization methods. The selection rules are AICc-hdf for BS, 
		AICc for LASSO and 10-fold CV for SparseNet, respectively.", scale.box = 0.7)


### Table 4: LBS and BS, selection rule is Cp-edf --------
# output: bs_lbs.tex
table.bs.lbs <- function(title, scale.box=0.7){
  # Parameters
  n = c(200, 2000)
  p = c(30, 180)
  snr = c(7, 0.2)
  names(snr) = c('hsnr', 'lsnr')
  nrep = 1000
  type = c('Orth-Sparse-Ex1', 'Orth-Sparse-Ex2', 'Orth-Dense')
  
  output = df.summary.result(n, p, snr, type, 
                             methods = list(c('bs', 'ic', 'cp', 'edf'), c('lbs', 'ic', 'cp', 'edf')),
                             type.table = 'maintext.method')

  output = cbind(rep(c("\\multirow{4}[4]{*}{n=200}", "", "\\cmidrule{2-9}", "", "\\midrule \\multirow{4}[4]{*}{n=2000}", "", "\\cmidrule{2-9}", ""), 3),
                 rep(c("\\multirow{2}[2]{*}{hsnr}", "", "\\multirow{2}[2]{*}{lsnr}", ""), 6),
                 rep(c("p=30", "p=180"), 12),
                 output)
  command = c(paste("\\toprule \n",
                    "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{2}{c|}{Orth-Sparse-Ex1} & \\multicolumn{2}{c|}{Orth-Sparse-Ex2} & \\multicolumn{2}{c|}{Orth-Dense} \\\\\n",
                    "\\cmidrule{4-9}\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & BS    & LBS   & BS    & LBS   & BS    & LBS  \\\\\n",
                    "\\midrule \n",
                    "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{6}{c|}{\\% worse than the best possible BS} \\\\\n",
                    "\\midrule \n"),
              paste("\\midrule \n",
                    "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{6}{c|}{Relative efficiency} \\\\\n",
                    "\\midrule \n"),
              paste("\\midrule \n",
                    "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{6}{c|}{Sparsistency (number of extra variables)} \\\\\n",
                    "\\midrule \n"),
              paste("\\bottomrule \n"
              )
  )
  print(xtable(output,
               align = "l|c|c|c|cc|cc|cc|",  # align and put a vertical line (first "l" again represents column of row numbers)
               label = 'tab:lbs_bs',
               caption = title),
        #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
        scalebox = scale.box,
        caption.placement = "top",
        include.rownames = FALSE, 
        include.colnames = FALSE, 
        hline.after = NULL, 
        floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
        sanitize.text.function = force, # Important to treat content of first column as latex function
        add.to.row = list(pos = list(-1,
                                     8,
                                     16,
                                     nrow(output)),
                          command = command
        ),
        file = paste0(base_tables, "/lbs_bs.tex"), compress = FALSE
  )
}
table.bs.lbs("The performance of BS and LBS. The selection rule for both methods is C$_p$-edf. We assume knowledge of $\\mu$ and $\\sigma$.")

### Table 5: BOSS compared to other methods --------
# output: boss_regu.tex
table.boss.regu <- function(title, scale.box=0.43){
  # Parameters
  n = c(200, 2000)
  p = c(30, 180)
  snr = c(7, 0.2)
  names(snr) = c('hsnr', 'lsnr')
  rho = c(0.5, 0.9)
  nrep = 1000
  type = c('Sparse-Ex3', 'Sparse-Ex4', 'Dense')
  
  # BOSS: AICc-hdf/CV, BS, FS, LASSO, SparseNet
  output = df.summary.result(n, p, snr, type, rho, 
                             methods = list(c('boss', 'ic', 'aicc', 'hdf'), c('boss', 'cv'), c('bs', 'cv'), 
                                            c('fs', 'cv'), c('lasso', 'ic', 'aicc'), c('sparsenet', 'cv')),
                             type.table = 'maintext.bossregu')

  output = cbind(rep(c("\\multirow{8}[4]{*}{n=200}", "",  "",  "", "\\cmidrule{2-19}", "",  "",  "", "\\midrule \\multirow{8}[4]{*}{n=2000}", "",  "",  "", "\\cmidrule{2-19}", "",  "",  ""), 3),
                 rep(c("\\multirow{4}[2]{*}{hsnr}", "",  "",  "", "\\multirow{4}[2]{*}{lsnr}", "",  "",  ""), 6),
                 rep(c("\\multirow{2}[1]{*}{$\\rho=0.5$}", "", "\\multirow{2}[1]{*}{$\\rho=0.9$}", ""), 12),
                 rep(c("p=30", "p=180"), 24),
                 output)
  command = c(paste("\\toprule \n",
                    "\\multicolumn{1}{|r}{} & \\multicolumn{1}{r}{} & \\multicolumn{1}{r}{} &       & \\multicolumn{5}{c|}{Sprse-Ex3}        & \\multicolumn{5}{c|}{Sparse-Ex4}       & \\multicolumn{5}{c|}{Dense} \\\\\n",
                    "\\cmidrule{5-19}\\multicolumn{1}{|r}{} & \\multicolumn{1}{r}{} & \\multicolumn{1}{r}{} &       & BOSS  & BS    & FS    & LASSO & SparseNet & BOSS  & BS    & FS    & LASSO & SparseNet & BOSS  & BS    & FS    & LASSO & \\multicolumn{1}{c|}{SparseNet}  \\\\\n",
                    "\\cmidrule{5-19}\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{15}{c|}{\\% worse than the best possible BOSS}  \\\\\n",
                    "\\midrule \n"),
              paste("\\midrule \n",
                    "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{15}{c|}{Relative efficiency} \\\\\n",
                    "\\midrule \n"),
              paste("\\midrule \n",
                    "\\multicolumn{1}{|c}{} & \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{} &       & \\multicolumn{15}{c|}{Sparsistency (number of extra variables)} \\\\\n",
                    "\\midrule \n"),
              paste("\\bottomrule \n"
              )
  )
  print(xtable(output,
               align = "l|c|c|c|c|ccccc|ccccc|ccccc|",  # align and put a vertical line (first "l" again represents column of row numbers)
               label = 'tab:boss_regu',
               caption = title),
        #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
        scalebox = scale.box,
        caption.placement = "top",
        include.rownames = FALSE, 
        include.colnames = FALSE, 
        hline.after = NULL, 
        floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
        sanitize.text.function = force, # Important to treat content of first column as latex function
        add.to.row = list(pos = list(-1,
                                     16,
                                     32,
                                     nrow(output)),
                          command = command
        ),
        file = paste0(base_tables, "/boss_regu.tex"), compress = FALSE
  )
}
table.boss.regu("The performance of BOSS compared to other methods. Selection rules are for 'AICc-hdf/CV' for BOSS, 
                AICc for LASSO and CV for the remaining methods in the table, respectively.")

### Table 6: Real data examples ---------
table.realdata <- function(title, scale.box){
  ## Read the results
  result = readRDS(paste0(base, '/code/run_model/realdata/results/boston_hitters_college_auto.rds'))
  result_ff = readRDS(paste0(base, '/code/run_model/realdata/results/forestfire.rds'))
  result = c(result[c('boston', 'hitters', 'auto', 'college')], result_ff)
  
  # Merge them into a data frame
  output = round(do.call(rbind, lapply(result, function(xx){do.call(rbind, lapply(xx, function(yy){do.call(c, lapply(yy, mean))}) )})), 3)
  output = output[,c('boss', 'bs', 'fs', 'lasso', 'sparsenet')]
  
  # Make the minimum value in each row bold
  tmp_function <- function(xx){
    ind = which(xx == min(xx))
    xx[ind] = paste0('\\textbf{', xx[ind], '}')
    xx
  }
  output = t(apply(output, 1, tmp_function))
  
  output = cbind(c("\\midrule\\multirow{3}[2]{*}{Housing (506, 13)}", "",  "", "\\midrule\\multirow{3}[2]{*}{Hitters (263, 19)}", "",  "", "\\midrule \\multirow{3}[2]{*}{Auto (392, 6)}", "",  "", "\\midrule\\multirow{3}[2]{*}{College (777, 17)}", "",  "",  "\\midrule\\multirow{3}[2]{*}{Forest Fires (517, 55)}", "",  ""),
                 rep(c("RMSE", "\\# predictors",  "running time (s)"), 5),
                 output)
  command = c(paste("\\toprule \n",
                    "Dataset (n, p) & Metrics & BOSS  & BS    & FS    & LASSO & SparseNet  \\\\\n"),
              paste("\\bottomrule \n"
              )
  )
  print(xtable(output,
               align = "l|c|c|ccc|c|c|",  # align and put a vertical line (first "l" again represents column of row numbers)
               label = 'tab:realdata',
               caption = title),
        #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
        scalebox = scale.box,
        caption.placement = "top",
        include.rownames = FALSE, 
        include.colnames = FALSE, 
        hline.after = NULL, 
        floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
        sanitize.text.function = force, # Important to treat content of first column as latex function
        add.to.row = list(pos = list(-1,
                                     nrow(output)),
                          command = command
        ),
        file = paste0(base_tables, "/realdata.tex"), compress = FALSE
  )
}
table.realdata(title = "Performance of subset selection methods on real datasets. 
               The results are averages of leave-one-out (LOO) testing. 
               The selection rules are AICc-hdf for BOSS, AICc for LASSO and 10-fold CV for the rest, respectively. 
               The intercept term is always fitted and is not counted in the number of predictors. 
               Minimal values for the metrics for each dataset are given in bold face."
              , scale.box = 0.9)
### Supplemental material: BS, selection rules --------
table.bs.selectrule.supplement <- function(scale.box = 0.55){
  # Parameters
  p = c(14, 30, 60, 180)
  snr = c(7, 1.5, 0.2)
  names(snr) = c('hsnr', 'msnr', 'lsnr')
  nrep = 1000
  
  for(type in c(paste0('Orth-Sparse-Ex', 1:2), 'Orth-Dense')){
    for(n in c(200, 2000)){
      title = paste0('The performance of BS using different selection rules, ', type, ', n=', n)
      filename = paste0(type, '_n', n)
      
      output = df.summary.result(n, p, snr, type, 
                                 methods = do.call(c, list(lapply(c('edf', 'ndf', 'hdf', 'bdf'), function(xx){c('bs', 'ic', 'cp', xx)}),
                                                           lapply(c('edf', 'ndf', 'hdf', 'bdf'), function(xx){c('bs', 'ic', 'aicc', xx)}),
                                                           lapply(c('edf', 'ndf', 'hdf', 'bdf'), function(xx){c('bs', 'ic', 'bic', xx)}),
                                                           lapply(c('edf', 'ndf', 'hdf', 'bdf'), function(xx){c('bs', 'ic', 'gcv', xx)}),
                                                           list(c('bs', 'cv')))),
                                 type.table = 'supplement.bsrule')
      
      output = cbind(rep(c("\\midrule\\multirow{4}[2]{*}{hsnr}", "", "", "", "\\midrule\\multirow{4}[2]{*}{msnr}", "", "", "", "\\midrule\\multirow{4}[2]{*}{lsnr}", "", "", ""), 3),
                     rep(paste0('p=', c(14,30,60,180)), 9),
                     output)
      command = c(paste("\\toprule \n",
                        "\\multicolumn{1}{|c}{} &       & \\multicolumn{2}{c|}{C$_p$} & \\multicolumn{2}{c|}{AICc} & \\multicolumn{2}{c|}{BIC} & \\multicolumn{2}{c|}{GCV} & \\multirow{2}[2]{*}{CV} \\\\\n",
                        "\\multicolumn{1}{|c}{} &       & edf   & ndf/hdf/bdf & edf   & ndf/hdf/bdf & edf   & ndf/hdf/bdf & edf   & ndf/hdf/bdf &       \\\\\n",
                        "\\cmidrule{3-11}\\multicolumn{1}{|c}{} &       & \\multicolumn{9}{c|}{\\% worse than the best possible BS} \\\\\n"),
                  paste("\\midrule \n",
                        "\\multicolumn{1}{|c}{} &       & \\multicolumn{9}{c|}{Relative efficiency} \\\\\n"),
                  paste("\\midrule \n",
                        "\\multicolumn{1}{|c}{} &       & \\multicolumn{9}{c|}{Sparsistency (number of extra variables)} \\\\\n"),
                  paste("\\bottomrule \n"
                  )
      )
      print(xtable(output,
                   align = "l|c|c|cc|cc|cc|cc|c|",  # align and put a vertical line (first "l" again represents column of row numbers)
                   label = NULL,
                   caption = title),
            #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
            scalebox = scale.box,
            caption.placement = "top",
            include.rownames = FALSE, 
            include.colnames = FALSE, 
            hline.after = NULL, 
            floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
            sanitize.text.function = force, # Important to treat content of first column as latex function
            add.to.row = list(pos = list(-1,
                                         12,
                                         24,
                                         nrow(output)),
                              command = command
            ),
            file = paste0(base_tables, "/supplement/bs_ic/", filename, ".tex"), compress = FALSE
      )
      
      
    }
  }
  
  

  
  
}
table.bs.selectrule.supplement()

### Supplemental material: BS and regularization methods --------
table.bs.regu.supplement <- function(scale.box = 0.75){
  # Parameters
  p = c(14, 30, 60, 180)
  snr = c(7, 1.5, 0.2)
  names(snr) = c('hsnr', 'msnr', 'lsnr')
  nrep = 1000
  
  for(type in c(paste0('Orth-Sparse-Ex', 1:2), 'Orth-Dense')){
    for(n in c(200, 2000)){
      title = paste0('The performance of BS compared to regularization methods, ', type, ', n=', n)
      filename = paste0(type, '_n', n)
      
      output = df.summary.result(n, p, snr, type, 
                                 methods = list(c('bs', 'ic', 'aicc', 'hdf'), c('lasso', 'ic', 'aicc'), c('lasso', 'cv'), 
                                                c('gamlr', 'ic', 'aicc'), c('gamlr', 'cv'), c('sparsenet', 'cv'), c('relaxlasso', 'cv')),
                                 type.table = 'supplement.bsregu')
      
      output = cbind(rep(c("\\midrule\\multirow{4}[2]{*}{hsnr}", "", "", "", "\\midrule\\multirow{4}[2]{*}{msnr}", "", "", "", "\\midrule\\multirow{4}[2]{*}{lsnr}", "", "", ""), 3),
                     rep(paste0('p=', c(14,30,60,180)), 9),
                     output)
      command = c(paste("\\toprule \n",
                        "\\multicolumn{1}{|c}{} &       & BS    & LASSO & Gamma LASSO & SparseNet & \\multicolumn{1}{c|}{rLASSO}  \\\\\n",
                        "\\multicolumn{1}{|c}{} &       & AICc-hdf & AICc/CV & AICc/CV & CV    & \\multicolumn{1}{c|}{CV}       \\\\\n",
                        "\\midrule\\multicolumn{1}{|c}{} &       & \\multicolumn{5}{c|}{\\% worse than the best possible BS} \\\\\n"),
                  paste("\\midrule \n",
                        "\\multicolumn{1}{|c}{} &       & \\multicolumn{5}{c|}{Relative efficiency} \\\\\n"),
                  paste("\\midrule \n",
                        "\\multicolumn{1}{|c}{} &       & \\multicolumn{5}{c|}{Sparsistency (number of extra variables)} \\\\\n"),
                  paste("\\bottomrule \n"
                  )
      )
      print(xtable(output,
                   align = "l|c|c|ccccc|",  # align and put a vertical line (first "l" again represents column of row numbers)
                   label = NULL,
                   caption = title),
            #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
            scalebox = scale.box,
            caption.placement = "top",
            include.rownames = FALSE, 
            include.colnames = FALSE, 
            hline.after = NULL, 
            floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
            sanitize.text.function = force, # Important to treat content of first column as latex function
            add.to.row = list(pos = list(-1,
                                         12,
                                         24,
                                         nrow(output)),
                              command = command
            ),
            file = paste0(base_tables, "/supplement/bs_regu/", filename, ".tex"), compress = FALSE
      )
      
      
    }
  }
  
  
}
table.bs.regu.supplement()

### Supplemental material: BOSS, FS, BS and regularization methods --------
# output: supplement/boss/*.tex
table.boss.regu.supplement <- function(scale.box=0.7){
  
  # Parameters
  p = c(14, 30, 60, 180)
  snr = c(7, 1.5, 0.2)
  names(snr) = c('hsnr', 'msnr', 'lsnr')
  nrep = 1000
  
  count = 0
  for(type in c(paste0('Sparse-Ex', 1:4), 'Dense')){
    for(rho in c(0, 0.5, 0.9)){
     for(n in c(200, 2000)){
       count = count + 1
       if(count %% 5 == 0){
         print(paste0(count,'/',30,' are done'))
       }
       title = paste0('The performance of BOSS compared to other methods, ', type, ', $\\rho$=', rho, ', n=', n)
       filename = paste0(type, '_n', n, '_rho_', gsub("[.]","",as.character(rho)))
       
       output = df.summary.result(n, p, snr, type, rho, 
                                  methods = list(c('boss', 'ic', 'cp', 'hdf'), c('boss', 'ic', 'aicc', 'hdf'), c('boss', 'cv'), 
                                                 c('bs', 'cv'), c('fs', 'cv'), c('lasso', 'ic', 'aicc'), c('lasso', 'cv'), 
                                                 c('gamlr', 'ic', 'aicc'), c('gamlr', 'cv'), c('sparsenet', 'cv'), c('relaxlasso', 'cv')),
                                  type.table = 'supplement.bossregu')
       
       output = cbind(rep(c("\\midrule\\multirow{4}[2]{*}{hsnr}", "",  "",  "", "\\midrule\\multirow{4}[2]{*}{msnr}", "",  "",  "", "\\midrule\\multirow{4}[2]{*}{lsnr}", "",  "",  ""), 3),
                      rep(paste0('p=',p), 9),
                      output)
       command = c(paste("\\toprule \n",
                         "\\multicolumn{1}{|c}{} &       & BOSS  & BS    & FS    & LASSO & Gamma LASSO & SparseNet & \\multicolumn{1}{c|}{rLASSO} \\\\\n",
                         "\\multicolumn{1}{|c}{} &       & C$_p$-hdf/AICc-hdf/CV & CV    & CV    & AICc/CV & AICc/CV & CV    & \\multicolumn{1}{c|}{CV}  \\\\\n",
                         "\\cmidrule{3-9}\\multicolumn{1}{|c}{} &       & \\multicolumn{7}{c|}{\\% worse than the best possible BOSS}  \\\\\n"
                         ),
                   paste("\\midrule \n",
                         "\\multicolumn{1}{|c}{} &       & \\multicolumn{7}{c|}{Relative efficiency} \\\\\n"),
                   paste("\\midrule \n",
                         "\\multicolumn{1}{|c}{} &       & \\multicolumn{7}{c|}{Sparsistency (number of extra variables)} \\\\\n"),
                   paste("\\bottomrule \n"
                   )
       )
       print(xtable(output,
                    align = "l|c|c|ccccccc|",  # align and put a vertical line (first "l" again represents column of row numbers)
                    label = NULL,
                    caption = title),
             #size = size, #Change size; useful for bigger tables "normalsize" "footnotesize"
             scalebox = scale.box,
             caption.placement = "top",
             include.rownames = FALSE, 
             include.colnames = FALSE, 
             hline.after = NULL, 
             floating = TRUE, # whether \begin{Table} should be created (TRUE) or not (FALSE)
             sanitize.text.function = force, # Important to treat content of first column as latex function
             add.to.row = list(pos = list(-1,
                                          12,
                                          24,
                                          nrow(output)),
                               command = command
             ),
             file = paste0(base_tables, "/supplement/boss/", filename, ".tex"), compress = FALSE
       )
       
     }
    }
  }
}
table.boss.regu.supplement()


