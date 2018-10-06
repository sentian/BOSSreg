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
cv.boss <- function(x, y, intercept=FALSE, fs.only=FALSE, hdf.ic.boss=TRUE, n.folds=10, n.rep=1){
  n = dim(x)[1]
  p = dim(x)[2]
  maxstep = trunc(min(n - n/n.folds, p))
  if(maxstep < p){
    warning('the number of observations in each fold, does not allow evaluating the full path, some large sets of variables are ignored')
  }

  # matrix to store the CV error
  cv_rep_boss = NULL
  if(!fs.only){
    cv_rep_boss = cv_rep_fs = matrix(NA, nrow=n.rep, ncol=maxstep+1)
  }

  for(replication in 1:n.rep){
    fold.index = sample(rep(1:n.folds, length.out=n)) # randomly assign a fold to each observation
    cv_tmp_boss = NULL
    if(!fs.only){
      cv_tmp_boss = cv_tmp_fs = matrix(NA, nrow=n.folds, ncol=maxstep+1)
    }
    for(fold in 1:n.folds){
      # split the training and testing sets
      test.index = which(fold.index==fold)
      x.test = x[test.index, ]
      y.test = y[test.index]
      x.train = x[-test.index, ]
      y.train = y[-test.index]
      boss_result <- boss(x.train, y.train, intercept, fs.only, hdf.ic.boss=FALSE)
      beta_fs= boss_result$beta_fs
      beta_boss = boss_result$beta_boss
      cv_tmp_fs[fold, ] = Matrix::colMeans((matrix(rep(y.test,each=maxstep+1), ncol=maxstep+1, byrow=T) - x.test%*%beta_fs)^2)
      if(!fs.only){
        cv_tmp_boss[fold, ] = Matrix::colMeans((matrix(rep(y.test,each=maxstep+1), ncol=maxstep+1, byrow=T) - x.test%*%beta_boss)^2)
      }
    }
    cv_rep_fs[replication, ] = Matrix::colMeans(cv_tmp_fs)
    if(!fs.only){
      cv_rep_boss[replication, ] = Matrix::colMeans(cv_tmp_boss)
    }
  }

  cv_boss=NULL
  cv_fs = Matrix::colMeans(cv_rep_fs)
  if(!fs.only){
    cv_boss = Matrix::colMeans(cv_rep_boss)
  }

  # fit on the full sample
  boss_result <- boss(x, y, intercept, fs.only, hdf.ic.boss)

  return(list(boss=boss_result, n.folds=n.folds, cvm.fs=cv_fs, cvm.boss=cv_boss, i.min.fs=which.min(cv_fs), i.min.boss=which.min(cv_boss)))
}
