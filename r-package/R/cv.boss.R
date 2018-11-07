#' Cross validation for BOSS.
#'
#' Cross validation for BOSS and FS.
#'
#' @usage cv.boss(x, y, n.folds=10, n.rep=1, ...)
#'
#' @param x A matrix of predictors, see \code{boss}.
#' @param y A vector of response variable, see \code{boss}.
#' @param n.folds The number of cross validation folds.
#' @param n.rep The number of replications of cross validation.
#' @param ... Arguments to \code{boss}.
#'
#' @return
#' \itemize{
#'   \item boss: An object \code{boss} that fits on the full dataset.
#'   \item n.folds: The number of cross validation folds.
#'   \item cvm.fs: Mean OOS deviance for each candidate given by FS.
#'   \item cvm.boss: Mean OSS deviance for each candidate given by BOSS.
#'   \item i.min.fs: The index of minimum cvm.fs.
#'   \item i.min.boss: The index of minimum cvm.boss.
#' }
#'
#' @details This function fits BOSS and FS (\code{boss}) on the full dataset, and performs \code{n.folds}
#'   cross validation. The cross validation process can be repeated \code{n.rep} times to evaluate the
#'   out-of-sample (OOS) performance for the candidate subsets given by both methods.
#'
#' @author Sen Tian
#' @references Tian, Hurvich and Simonoff (2018), On the use of information criterion
#'   in least squares based subset selection problems. (Link to be added)
#' @example R/example/eg.cv.boss.R
#' @export
cv.boss <- function(x, y, n.folds=10, n.rep=1, ...){
  # # arguments
  argu = list(...)
  # argu_boss = c('intercept', 'fs.only', 'hdf.ic.boss') # arguments that boss accepts
  # # arguments that user specify but unused
  # argu_unused = setdiff(names(argu), argu_boss)
  # if(length(argu_unused) > 0){
  #   warning(paste(argu_unused, ' are not valid arguments for boss, check spelling maybe?', sep=''))
  # }

  if(is.null(argu$fs.only)) argu$fs.only = FALSE

  # overide hdf.ic.boss option in '...', to be used in CV
  boss.nohdf <- function(x, y, ..., hdf.ic.boss) boss(x, y, ..., hdf.ic.boss=FALSE)

  # start the CV process
  n = dim(x)[1]
  p = dim(x)[2]
  maxstep = trunc(min(n - n/n.folds, p))
  if(maxstep < p){
    warning('the number of observations in each fold, does not allow evaluating the full path, some large sets of variables are ignored')
  }

  # matrix to store the CV error
  cv_rep_boss = NULL
  if(!argu$fs.only){
    cv_rep_boss = cv_rep_fs = matrix(NA, nrow=n.rep, ncol=maxstep+1)
  }

  for(replication in 1:n.rep){
    fold.index = sample(rep(1:n.folds, length.out=n)) # randomly assign a fold to each observation
    cv_tmp_boss = NULL
    if(!argu$fs.only){
      cv_tmp_boss = cv_tmp_fs = matrix(NA, nrow=n.folds, ncol=maxstep+1)
    }
    for(fold in 1:n.folds){
      # split the training and testing sets
      test.index = which(fold.index==fold)
      x.test = x[test.index, , drop=FALSE]
      y.test = y[test.index]
      x.train = x[-test.index, , drop=FALSE]
      y.train = y[-test.index]
      boss_result <- boss.nohdf(x.train, y.train, ...)
      beta_fs= boss_result$beta_fs
      beta_boss = boss_result$beta_boss
      # if intercept
      if(dim(beta_fs)[1] == p+1){
        x.test = cbind(rep(1,nrow(x.test)), x.test)
      }
      cv_tmp_fs[fold, ] = Matrix::colMeans((matrix(rep(y.test,each=maxstep+1), ncol=maxstep+1, byrow=T) - x.test%*%beta_fs)^2)
      if(!argu$fs.only){
        cv_tmp_boss[fold, ] = Matrix::colMeans((matrix(rep(y.test,each=maxstep+1), ncol=maxstep+1, byrow=T) - x.test%*%beta_boss)^2)
      }
    }
    cv_rep_fs[replication, ] = Matrix::colMeans(cv_tmp_fs)
    if(!argu$fs.only){
      cv_rep_boss[replication, ] = Matrix::colMeans(cv_tmp_boss)
    }
  }

  cv_boss=NULL
  cv_fs = Matrix::colMeans(cv_rep_fs)
  if(!argu$fs.only){
    cv_boss = Matrix::colMeans(cv_rep_boss)
  }

  # fit on the full sample
  boss_result <- boss(x, y, ...)

  # output
  out = list(boss=boss_result,
             n.folds=n.folds,
             cvm.fs=cv_fs,
             cvm.boss=cv_boss,
             i.min.fs=which.min(cv_fs),
             i.min.boss=which.min(cv_boss))
  class(out) = 'cv.boss'
  invisible(out)
}


#' Select coefficient vector based on cross validation (CV).
#'
#' This function returns coefficient vector that minimizes out-of-sample (OOS) cross
#' validation score.
#'
#' @param object The cv.boss object, returned from calling 'cv.boss' function.
#' @param ... Extra arguments (unused for now).
#'
#' @return
#' \itemize{
#'   \item fs: The chosen coefficient vector for FS.
#'   \item boss: The chosen coefficient vector for FS.
#' }
#'
#' @examples
#' # See the example in the section of \code{cv.boss}. Or type ?cv.boss in R.
#'
#' @importFrom stats coef
#' @export
coef.cv.boss <- function(object, ...){
  coef_result = coef(object$boss, select.fs=object$i.min.fs, select.boss=object$i.min.boss)
  beta_fs_opt = coef_result$fs
  beta_boss_opt = coef_result$boss
  return(list(fs=beta_fs_opt, boss=beta_boss_opt))
}

#' Prediction given new data entries.
#'
#' This function returns the prediction(s) given new observation(s), for FS and BOSS,
#' where the optimal coefficient vector is chosen via cross-validation.
#'
#' @param object The cv.boss object, returned from calling 'cv.boss' function.
#' @param newx A new data entry or several entries. It can be a vector, or a matrix with
#' \code{nrow(newx)} being the number of new entries and \code{ncol(newx)=p} being the
#' number of predictors. The function takes care of the intercept, NO need to add \code{1}
#' to \code{newx}.
#' @param ... Extra arguments (unused for now).
#'
#' @return
#' \itemize{
#'   \item fs: The prediction for FS.
#'   \item boss: The prediction for BOSS.
#' }
#'
#' @examples
#' # See the example in the section of \code{cv.boss}. Or type ?cv.boss in R.
#'
#' @importFrom stats predict
#' @export
predict.cv.boss <- function(object, newx, ...){
  predict_result = predict(object$boss, newx, select.fs=object$i.min.fs, select.boss=object$i.min.boss)
  mu_fs_opt = predict_result$fs
  mu_boss_opt = predict_result$boss
  return(list(fs=mu_fs_opt, boss=mu_boss_opt))
}
