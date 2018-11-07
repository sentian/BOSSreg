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
#' @param mu True mean vector, used in the calculation of hdf. Default is NULL, and is estimated via full OLS.
#' @param sigma True standard deviation of the error, used in the calculation of hdf. Default is NULL,
#'   and is estimated via full OLS.
#' @param use.lasso.sigma Whether to use LASSO (10 fold CV) based estimates of sigma, rather than full OLS,
#'   to be plugged into hdf. Details about LASSO based sigma can be found in Reid et al. (2016).
#' @param ... Extra parameters to allow flexibility. Currently none argument allows or requires, just for
#'   the convinience of call from other parent functions like cv.boss.
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
#' @references Reid, Tibshirani and Friedman (2016), A study of error variance estimation
#'   in LASSO regression.
#' @example R/example/eg.boss.R
#' @export
boss <- function(x, y, intercept=TRUE, fs.only=FALSE, hdf.ic.boss=TRUE, mu=NULL, sigma=NULL, use.lasso.sigma=FALSE, ...){
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
    if(intercept){
      beta_boss = rbind((mean_y - mean_x %*% beta_boss), beta_boss)
    }
    beta_boss = Matrix::Matrix(beta_boss, sparse=T)

    # hdf and IC
    if(n<=p & hdf.ic.boss){
      warning('hdf not available when n<=p')
    }else if(hdf.ic.boss){
      hdf_result = calc.hdf(Ql, y, mu, sigma)
      if(intercept){
        hdf_result$hdf = hdf_result$hdf + 1
      }
      if(is.null(sigma)){
        IC_result = calc.ic.all(beta_q, Ql, y, hdf_result$hdf, hdf_result$sigma)
      }else{
        IC_result = calc.ic.all(beta_q, Ql, y, hdf_result$hdf, sigma)
      }

    }

  }
  # output
  out = list(beta_fs=beta_fs,
             beta_boss=beta_boss,
             steps_fs=steps,
             hdf_boss=hdf_result$hdf,
             IC_boss=IC_result,
             call=list(intercept=intercept, fs.only=fs.only))
  class(out) = 'boss'
  invisible(out)
}





#' Select coefficient vector(s) from boss object.
#'
#' This function returns coefficient vector(s) for given step(s) of FS and BOSS.
#' For BOSS, it can also return the optimal coefficient vector selected by AICc
#' (by default) or other IC.
#'
#' @param object The boss object, returned from calling 'boss' function.
#' @param select.fs Given step(s) of FS, corresponding to columns in the coefficient
#' matrix. For example, the first column in beta_fs corresponds to step 1, the second
#' column corresponds to step 2, etc. We enforce the first column, that has all 0 entries,
#' to be step 1, just for convinicence of usage.
#' @param select.boss Given columns(s) in beta_boss, similar to select.fs.
#' @param method.boss Which IC is used to select the optimal coefficient vector for BOSS.
#' @param ... Extra arguments (unused for now)
#'
#' @return
#' \itemize{
#'   \item fs: The chosen coefficient vector(s) for FS.
#'   \item boss: The chosen coefficient vector(s) for FS.
#' }
#' @details If \code{select.fs} or \code{select.boss} is specified, the function returns
#' corresponding column(s) in the coefficient matrix.
#'
#' If \code{select.fs} is unspecified,
#' we return the entire coefficient matrix (\code{beta_fs}) for FS, since we currently do
#' not have a hands-on way of calculating IC for FS. But user can use their own rules and
#' specify e.g. \code{select.fs=which.min(rule)} in order to pick the optimal coefficient
#' vector.
#'
#' If \code{select.boss} is unspecified, we will try to return the optimal coefficient
#' vector selected by AICc-hdf (other choice of IC can be specified in \code{method.boss}).
#' The only exception is when n>=p, where hdf is not well defined, and we will return the
#' entire coefficient matrix.
#'
#' @examples
#' # See the example in the section of \code{boss}. Or type ?boss in R.
#'
#' @importFrom stats coef
#' @export
coef.boss <- function(object, select.fs=NULL, select.boss=NULL, method.boss=c('aicc','bicc','aic','bic','gcv','cp'), ...){
  # for fs, return the full coef matrix if not specified the columns
  if(is.null(select.fs)){
    select.fs = 1:ncol(object$beta_fs)
  }else if(select.fs == 0){
    select.fs = 1:ncol(object$beta_fs)
  }
  select.fs[select.fs > ncol(object$beta_fs)] = ncol(object$beta_fs)
  beta_fs_opt = object$beta_fs[, select.fs, drop=FALSE]

  # for boss, the default is to return coef selected by AICc
  # if fs only
  if(object$call$fs.only){
    beta_boss_opt = NULL
  }else{
    if(is.null(select.boss)){
      if(is.null(object$IC_boss)){
        # if we are in the case where n>=p
        if(dim(object$beta_boss)[1]+1 >= dim(object$beta_boss)[2]){
          warning('hdf does not work when n>=p, full coef matrix returned')
        }else{
          # this is where hdf.ic.boss is flagged FALSE
          warning("rerun boss with argument 'hdf.ic.boss=TRUE', and call coef.boss again")
        }
        select.boss = 1:ncol(object$beta_boss)
      }else{
        method.boss = match.arg(method.boss)
        select.boss = which.min(object$IC_boss[[method.boss]])
      }
    }else if(select.boss == 0){
      select.boss = 1:ncol(object$select.boss)
    }
    select.boss[select.boss > ncol(object$beta_boss)] = ncol(object$beta_boss)
    beta_boss_opt = object$beta_boss[, select.boss, drop=FALSE]
  }
  return(list(fs=beta_fs_opt, boss=beta_boss_opt))
}


#' Prediction given new data entries.
#'
#' This function returns the prediction(s) given new observation(s), for FS and BOSS,
#' where the optimal coefficient vector is chosen via certain selection rule.
#'
#' @param object The boss object, returned from calling 'boss' function.
#' @param newx A new data entry or several entries. It can be a vector, or a matrix with
#' \code{nrow(newx)} being the number of new entries and \code{ncol(newx)=p} being the
#' number of predictors. The function takes care of the intercept, NO need to add \code{1}
#' to \code{newx}.
#' @param ... Extra arguments to be plugged into \code{coef}, such as \code{select.boss},
#' see the description of \code{coef} for more details.
#'
#' @return
#' \itemize{
#'   \item fs: The prediction(s) for FS.
#'   \item boss: The prediction(s) for BOSS.
#' }
#'
#' @details The function basically calculates \eqn{x * coef}, where \code{coef}
#' is a coefficient vector chosen by a selection rule. See more details about the default
#' and available choices of the selection rule in the description of \code{coef.boss}.
#'
#' @examples
#' #See the example in the section of \code{boss}. Or type ?boss in R.
#'
#' @importFrom stats predict
#' @export
predict.boss <- function(object, newx, ...){
  # coefficients
  coef_result = coef(object, ...)
  beta_fs_opt = coef_result$fs
  beta_boss_opt = coef_result$boss

  # make newx a matrix
  # if newx is an array or a column vector, make it a row vector
  if(is.null(dim(newx))){
    newx = matrix(newx, nrow=1)
  }else if(dim(newx)[2] == 1){
    newx = t(newx)
  }
  # if intercept, add 1 to newx
  if(object$call$intercept){
    newx = cbind(rep(1,nrow(newx)), newx)
  }

  # check the dimension
  if(ncol(newx) != nrow(beta_fs_opt)){
    stop('Mismatch dimension of newx and coef for FS. Note do NOT add 1 to newx when intercept=TRUE')
  }else{
    mu_fs_opt = newx %*% beta_fs_opt
  }
  if(is.null(beta_boss_opt)){
    mu_boss_opt = NULL
  }else{
    if(ncol(newx) != nrow(beta_boss_opt)){
      stop('Mismatch dimension of newx and coef for BOSS. Note do NOT add 1 to newx when intercept=TRUE')
    }else{
      mu_boss_opt = newx %*% beta_boss_opt
    }
  }
  return(list(fs=mu_fs_opt, boss=mu_boss_opt))

}








