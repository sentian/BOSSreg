library(JuliaCall)
library(BOSSreg)
library(bestsubset)
library(MASS)

# julia = julia_setup() 
# julia = julia_setup(JULIA_HOME="/share/apps/julia/1.5.3/bin")
julia = julia_setup(JULIA_HOME='/Applications/Julia-1.5.app/Contents/Resources/julia/bin')

# julia_install_package_if_needed("Distributions") # if you don't already have the package installed
julia_library("SubsetSelection")
julia_library("SubsetSelectionCIO")
julia_library("StatsBase")

# n = 8
# julia_assign("n", n)
# julia_command("sqrt(n)")
# (a = julia_eval("sqrt(n)"))
# (a = julia_call("sqrt", n))
# (a = n %>J% sqrt())

gen.data.generalx <- function(n, p, rho, snr, type, nrep=1000, seed=66, print.r2=FALSE, center.y=TRUE){
  set.seed(seed)
  
  # Covariance matrix and true beta
  if(type == 'Sparse-Ex1'){
    covmatrix = rho^abs(outer(1:p, 1:p, '-'))
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


n = 200
p = 180
p0 = 6
data = gen.data.generalx(n, p, rho=0.5, snr=7, type='Sparse-Ex1', nrep=1, seed=66, print.r2=FALSE, center.y=TRUE)

indices = betahat = list()
ptm = proc.time()
for(k in 1:180){
	Sparse_Regressor = julia_call("subsetSelection", julia_eval("OLS()"), julia_call("Constraint", k), c(data$y), data$x)
	indices[[k]] = Sparse_Regressor$indices
	betahat[[k]] = Sparse_Regressor$w
}
proc.time() - ptm


indices = betahat = list()
ptm = proc.time()
for(k in 1:10){
  Sparse_Regressor = julia_call("oa_formulation", julia_eval("SubsetSelection.OLS()"), c(data$y), data$x, k, )
  indices[[k]] = Sparse_Regressor$indices
  betahat[[k]] = Sparse_Regressor$w
}
proc.time() - ptm


ptm = proc.time()
boss_model = boss(data$x, data$y, intercept=FALSE )
proc.time() - ptm

ptm = proc.time()
bs_model = bs(x, y, k=0:n, verbose=TRUE)
proc.time() - ptm

