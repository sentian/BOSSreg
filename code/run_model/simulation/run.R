### This file fits the models and evaluate them
### This file was called by run_orthx.sh and run_generalx.sh, in order to be run on an HPC server parallelly for all the configurations of true models
### Final results are stored in the directory '../results'

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
                          file = paste(base, "/run_model/simulation/para_forhpc/generalx_",count,".txt",sep=''), sep = '', col.names = FALSE, row.names = FALSE)
            }
          }
        }
      }
    }

  }


}

# base = '/Volumes/HDD/Dropbox/Sen/Research/Model_selection/BOSSreg/code'
# para.forhpc(TRUE)
# para.forhpc(FALSE)

###### Fit the models on HPC server --------
base = '/scratch/st1864/boss' # the basis of the directory to store all the outputs
## Environment parameters 
args = (commandArgs(TRUE))
orthx = as.logical( as.character(args[1]) )
num = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

if(orthx){
  pathname = 'orthx'
}else{
  pathname = 'generalx'
}

## create folders to store the results (both temporary and final results)
# folders for temporary results in case the code is terminated somewhere
dir.create(paste0(base, '/tmp/bdf_bs/', pathname), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(base, '/tmp/result_cv/', pathname), showWarnings = FALSE, recursive = TRUE)
# folder for the final results
dir.create(paste0(base, '/results/', pathname), showWarnings = FALSE, recursive = TRUE)


## Read the parameters
para = c(t(read.table(file = paste(base, "/para_forhpc/", pathname, "_", num, ".txt",sep=''), sep = '', header=FALSE)))
type = para[1]
n = as.numeric(para[2])
p = as.numeric(para[3])
snr = as.numeric(para[4])
if(!orthx){
  rho = as.numeric(para[5])
}
snr_all = c(0.2, 1.5, 7)
names(snr_all) = c('lsnr', 'msnr', 'hsnr')
snr_name = names(snr_all)[which(snr_all == snr)]

## Filename
filename = paste0(pathname, '/', type, '_n', n, '_p', p, '_', snr_name)
if(!orthx){
  filename = paste0(filename, '_rho', gsub("[.]","",as.character(rho)))
}

## Source all the required libraries and functions
## The code and results are stored at different file systems on the server
source(paste0(getwd(), '/utils.R'))

## Generate simulated datasets
if(orthx){
  data = gen.data.orthx(n=n, p=p, snr=snr, type=type)
}else{
  data = gen.data.generalx(n=n, p=p, rho=rho, snr=snr, type=type)
}

## Bootstrapped df for BS, only when X is orthogonal
if(orthx){
  print('bdf')
  bdf_bs = bdf.bs.orthx(data$x, data$y)
  saveRDS(bdf_bs, paste0(base, '/tmp/bdf_bs/', filename, '.rds'))
}else{
  bdf_bs = NULL
}

## Cross-valications
print('cv')
ptm = proc.time()
result_cv = run.cv(data$x, data$y, orthx=orthx)
proc.time() - ptm
saveRDS(result_cv, paste0(base, '/tmp/result_cv/', filename, '.rds'))


## Evaulate the selected subsets for all methods
print('evaluation')
result = eval.metrics(data$x, data$y, data$beta, data$sigma, result_cv, bdf_bs, orthx)

## When X is orthogonal, run LBS (Lagrangian BS) with Cp-edf
if(orthx){
  print('LBS')
  result_lbs = eval.metrics.lbs.cp.orthx(data$x, data$y, data$beta, data$sigma)
  result$lbs$ic$cp$edf = result_lbs
}

## Save the results
saveRDS(result, paste0(base, '/results/', filename, '.rds'))

