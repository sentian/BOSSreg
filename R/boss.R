#' Best orthogonalized subset selection (BOSS).
#'
#'\itemize{
#'  \item Compute the solution path of BOSS and forward stepwise selection (FS).
#'  \item Compute various information criteria based on a heuristic degrees of freedom
#'   that can serve as the selection rule to choose the optimal subset given by BOSS.
#'   Only work when n>p.
#'}
#' @param x A matrix of predictors, with \code{nrow(x)=length(y)=n} observations and
#'   \code{ncol(x)=p} predictors. Intercept shall not be included.
#' @param y A vector of response variable, with \code{length(y)=n}.
#' @param intercept Whether to include an intercept term.
#' @param fs.only Whether to ignore BOSS and perform FS only.
#' @param hdf.ic.boss Whether to calculate the heuristic degrees of freedom (hdf)
#'   and information criteria (IC) for BOSS. IC includes AIC, BIC, AICc, BICc, GCV,
#'   Cp. Note that if the option fs.only=TRUE or n<=p, \code{hdf.ic.boss=FALSE} no matter what.
#'
#' @return
#' \itemize{
#'   \item beta_fs: A matrix of regression coefficients for each step performed by FS,
#'   from a null model until stop, with \code{nrow=p} and \code{ncol=min(n,p)+1}, where min(n,p) is
#'   the maximum number of steps performed.
#'   \item beta_boss: A matrix of regression coefficients for each step performed by
#'   BOSS, with \code{nrow=p} and \code{ncol=min(n,p)+1}. Note that unlike beta_fs and due to the nature of BOSS,
#'   the number of non-zero components in columns of beta_boss may not be unique, i.e.
#'   there maybe multiple columns corresponding to the same size of subset. \code{beta_boss=NULL}
#'   if the option \code{fs.only=TRUE}.
#'   \item steps_fs: A vector of numbers representing which predictor joins at each step,
#'   with \code{length(steps_fs)=min(n,p)}.
#'   \item hdf_boss: A vector of heuristic degrees of freedom (hdf) for BOSS, with
#'   \code{length(hdf_boss)=p+1}. Note that \code{hdf_boss=NULL} if n<=p or \code{hdf.ic.boss=FALSE}.
#'   \item IC_boss: A list of information criteria (IC) for BOSS, where each element
#'   in the list is a vector representing values of a given IC for each candidate subset
#'   of BOSS (or each column in beta_boss). The output IC includes AIC, BIC, AICc, BICc,
#'   GCV and Mallows' Cp. Note that each IC is calculated by plugging in hdf_boss.
#'
#' }
#'
#' @details This function computes the full solution path given by FS and (or) BOSS on a given
#'   dataset (x,y) with n observations and p predictors. Meanwhile, in the case where n>p, it calculates
#'   the heuristic degrees of freedom for BOSS, and various information criteria, which can further
#'   be used to select the optimal candidate along the path. Please refer to the example section below
#'   for implementation details and Tian et al. (2018) for methodology details.
#'
#' @author Sen Tian
#' @references Tian, Hurvich and Simonoff (2018), On the use of information criterion
#'   in least squares based subset selection problems. (Link to be added)
#' @example R/example/eg.boss.R
#' @export
boss <- function(x, y, intercept=FALSE, fs.only=FALSE, hdf.ic.boss=TRUE){
  n = dim(x)[1]
  p = dim(x)[2]

  maxstep = min(n, p)

  # standardize x (mean 0 and norm 1) and y (mean 0)
  std_result = std(x, y)
  x = std_result$x_std
  colnames(x) = seq(1, p)
  y = std_result$y_std
  mean_x = std_result$mean_x
  mean_y = std_result$mean_y
  sd_demanedx = std_result$sd_demeanedx

  steps = rep(NA, maxstep) # the step in order of the variables
  # determine the first variable to step in
  steps[1] = which.max(abs(t(x) %*% y))
  #x_active = x[,steps[1]]
  x_remain = x[, -steps[1]]
  qr_result = qr(x[, steps[1]])
  Q = qr.Q(qr_result, complete=T)
  Ql = Q[, 1, drop=FALSE]
  Qr = Q[, -1, drop=FALSE]
  R = qr.R(qr_result, complete=F)

  resid_tmp = x_remain # residuals of regressing x_remain on Ql
  # following steps
  for(i in 2:maxstep){
    # determine which variable to step in, based on partial correlation
    if(i < p){
      resid_tmp = resid_tmp - Ql[, dim(Ql)[2]] %*% t(Ql[, dim(Ql)[2]]) %*% x_remain
      j_remain = which.max(abs(t(y) %*% scale(resid_tmp, center=F, scale=sqrt(colSums(resid_tmp^2)))))
    }else{
      j_remain = 1
    }

    steps[i] = as.numeric(colnames(x_remain)[j_remain])
    # update QR
    updateQR_result = updateQR(Ql, Qr, R, x[, steps[i]])
    Ql = updateQR_result$Ql
    Qr = updateQR_result$Qr
    R = updateQR_result$R

    # update others
    if(i < p){
      resid_tmp = resid_tmp[, -j_remain,drop=FALSE]
      x_remain = x_remain[, -j_remain,drop=FALSE]
    }
  }

  # coefficients
  z = t(Ql) %*% y

  # fs
  beta_q = matrix(rep(z, maxstep), nrow=maxstep, byrow=F)
  beta_q = beta_q * upper.tri(beta_q, diag=T)
  beta_q = cbind(0, beta_q)
  if(n<p){
    steps_expand = c(steps, setdiff(seq(1, p),steps))
    beta_fs = rbind(backsolve(R, beta_q)[order(steps),], matrix(0, nrow=p-n, ncol=maxstep+1))[order(steps_expand), ]
  }else{
    beta_fs = backsolve(R, beta_q)[order(steps), ]
  }
  # scale back
  beta_fs = diag(1/sd_demanedx) %*% beta_fs
  #beta_fs = cbind(0,beta_fs)
  if(intercept){
    beta_fs = rbind((mean_y - mean_x %*% beta_fs), beta_fs)
  }
  beta_fs = Matrix::Matrix(beta_fs, sparse=T)

  beta_boss = hdf_result = IC_result = NULL

  if(!fs.only){
    # boss
    beta_q = matrix(0, nrow=p, ncol=maxstep+1)
    for(j in 1:maxstep){
      beta_q[order(-z^2)[1:j], (j+1)] = z[order(-z^2)[1:j]]
    }
    # project back and change the order
    if(n<p){
      beta_boss = rbind(backsolve(R, beta_q), matrix(0, nrow=p-n, ncol=maxstep+1))[order(steps_expand), ]
    }else{
      beta_boss = backsolve(R, beta_q)[order(steps), ]
    }
    # scale back
    beta_boss = diag(1/sd_demanedx) %*% beta_boss
    #beta_boss = cbind(0,beta_boss)
    if(intercept){
      beta_boss = rbind((mean_y - mean_x %*% beta_boss), beta_boss)
    }
    beta_boss = Matrix::Matrix(beta_boss, sparse=T)

    # hdf and IC
    if(n<=p & hdf.ic.boss){
      warning('hdf not available when n<=p')
    }else if(hdf.ic.boss){
      hdf_result = calc.hdf(Ql, y)
      IC_result = calc.ic.all(beta_q, Ql, y, hdf_result$hdf, hdf_result$sigma)
    }

  }

  return(list(beta_fs=beta_fs, beta_boss=beta_boss, steps_fs=steps, hdf_boss=hdf_result$hdf, IC_boss=IC_result))
}





