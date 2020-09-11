### This file fits and evaluates the simplifed relaxed lasso from Hastie (2017)
### This file was called by run_orthx.sh and run_generalx.sh, in order to be run on an HPC server parallelly for all the configurations of true models
### Final results are combined with the main results and are stored in the directory '../results'

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
dir.create(paste0(base, '/tmp/result_cv/srlasso/', pathname), showWarnings = FALSE, recursive = TRUE)
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
# source(paste0('/home/st1864/boss/code/utils.R'))

## Generate simulated datasets
if(orthx){
  data = gen.data.orthx(n=n, p=p, snr=snr, type=type)
}else{
  data = gen.data.generalx(n=n, p=p, rho=rho, snr=snr, type=type)
}

## Cross-validations
print('cv')
ptm = proc.time()
result_cv = run.cv.srlasso(data$x, data$y)
proc.time() - ptm
saveRDS(result_cv, paste0(base, '/tmp/result_cv/srlasso/', filename, '.rds'))


## Evaluate the selected subsets for all methods
print('evaluation')
result = eval.metrics.srlasso(data$x, data$y, data$beta, result_cv)
result_all = readRDS(paste0(base, '/results/', filename, '.rds'))
result_all = c(result_all, result)
## Save the results
saveRDS(result_all, paste0(base, '/results/', filename, '.rds'))

