### This file fits and evaluates various subset selection models, where the dimension of X is larger than the sample size
### This file was called by run_generalx_highdim.sh, in order to be run on an HPC server parallelly for all the configurations of true models
### Final results are combined with the main results and are stored in the directory '../results'

para.forhpc.highdim <- function(){
  count = 360
  n = rep(200, 4)
  p = c(550, 1000, 4000, 10000)
  snr = c(0.2, 1.5, 7)
  type = c('Sparse-Ex1', 'Sparse-Ex2', 'Sparse-Ex3', 'Sparse-Ex4', 'Dense')
  rho = c(0, 0.5, 0.9)
  
  for(i in 1:length(type)){
    for(j in 1:length(n)){
      for(l in 1:length(snr)){
        for(m in 1:length(rho)){
          count = count + 1
          write.table(c(type[i], n[j], p[j], snr[l], rho[m]),
                      file = paste(base, "/run_model/simulation/para_forhpc_highdim/generalx_",count,".txt",sep=''), sep = '', col.names = FALSE, row.names = FALSE)
        }
      }
    }
  }
}

# base = '/Volumes/HDD/Dropbox/Sen/Research/Model_selection/BOSSreg/code'
# para.forhpc.highdim()

###### Fit the models on HPC server --------
base = '/scratch/st1864/boss' # the basis of the directory to store all the outputs
## Environment parameters 
num = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

pathname = 'generalx'

## create folders to store the results (both temporary and final results)
# folders for temporary results in case the code is terminated somewhere
dir.create(paste0(base, '/tmp/betahat/', pathname), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(base, '/tmp/result_intermediate/', pathname), showWarnings = FALSE, recursive = TRUE)
# folder for the final results
dir.create(paste0(base, '/results/', pathname), showWarnings = FALSE, recursive = TRUE)


## Read the parameters
para = c(t(read.table(file = paste(base, "/para_forhpc/", pathname, "_", num, ".txt",sep=''), sep = '', header=FALSE)))
type = para[1]
n = as.numeric(para[2])
p = as.numeric(para[3])
snr = as.numeric(para[4])
rho = as.numeric(para[5])

snr_all = c(0.2, 1.5, 7)
names(snr_all) = c('lsnr', 'msnr', 'hsnr')
snr_name = names(snr_all)[which(snr_all == snr)]

## Filename
filename = paste0(type, '_n', n, '_p', p, '_', snr_name, '_rho', gsub("[.]","",as.character(rho)))

## Source all the required libraries and functions
## The code and results are stored at different file systems on the server
source(paste0(getwd(), '/utils.R'))

# ## Generate simulated datasets
# data = gen.data.generalx(n=n, p=p, rho=rho, snr=snr, type=type)

# ## Cross-valications
# print('cv')
# ptm = proc.time()
# result_cv = run.cv.simplified(data$x, data$y)
# proc.time() - ptm
# saveRDS(result_cv, paste0(base, '/tmp/result_cv/', filename, '.rds'))
# 
# ## Evaulate the selected subsets for all methods
# print('evaluation')
# result = eval.metrics.simplified(data$x, data$y, data$beta, data$sigma, result_cv)

result = run.all.simplified(n, p, rho, snr, type, nrep=1000, p0=54, write.tmp.to.file = TRUE)
## Save the results
saveRDS(result, paste0(base, '/results/', pathname, '/', filename, '.rds'))








