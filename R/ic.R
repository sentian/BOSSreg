#' Calculate information criterion.
#'
#' Calculate a specified information criterion (IC) for an estimate or a group of estimates.
#' Such IC includes AIC, BIC, AICc, BICc, GCV and Mallows' Cp.
#'
#' @param y_hat A vector of fitted values with \code{length(y_hat)=length(y)=n}, or
#'   a matrix, with \code{nrow(coef)=length(y)=n} and \code{ncol(y_hat)=m}, containing m different fits.
#' @param y A vector of response variable, with \code{length(y)=n}.
#' @param method A specified IC to calculate. Default is AICc ('aicc'). Other choices include AIC ('aic'),
#'   BIC ('bic'), BICc ('bicc'), GCV ('gcv') and Mallows' Cp ('cp').
#' @param df A number if y_hat is a vector, or a vector with \code{length(df)=ncol(y_hat)=m} if y_hat is
#'   a matrix. df represents the degrees of freedom for each fit.
#' @param sigma Standard deviation of the error term. It only needs to be specified if \code{method='cp'}.
#'
#' @return The value(s) of the specified IC for each fit.
#'
#' @details This function enables the computation of various common IC for model fits, which can
#'   further be used to choose the optimal fit. This allows user comparing the effect of different IC.
#'   In order to calculate IC, df needs to be specified. To be more specific, here are the formulas used
#'   to calculate each IC:
#'
#'   \deqn{AIC = \log(\frac{RSS}{n}) + 2\frac{df}{n}}{AIC = log(RSS/n) + 2*df/n}
#'   \deqn{BIC = \log(\frac{RSS}{n}) + \log(n)\frac{df}{n}}{BIC = log(RSS/n) + log(n)*df/n}
#'   \deqn{AICc = \log(\frac{RSS}{n}) + 2\frac{df+1}{n-df-2}}{AICc = log(RSS/n) + 2*(df+1)/(n-df-2)}
#'   \deqn{BICc = \log(\frac{RSS}{n}) + \log(n)\frac{df+1}{n-df-2}}{BICc = log(RSS/n) + log(n)*(df+1)/(n-df-2)}
#'   \deqn{GCV = \frac{RSS}{(n-df)^2}}{GCV = RSS/(n-df)^2}
#'   \deqn{Mallows' Cp = RSS + 2\times \sigma^2 \times df}{AIC = RSS + 2*\sigma^2*df}
#'
#' @author Sen Tian
#' @example R/example/eg.ic.R
#' @export
calc.ic <- function(y_hat, y, method='aicc', df, sigma=NULL){
  if(! method %in% c('aic','bic','aicc','bicc','gcv','cp')){
    stop("the IC specified is not among ('aic','bic','aicc','bicc','gcv','cp')")
  }
  # unify dimensions
  y = matrix(y, ncol=1)
  df = matrix(df, nrow=1)
  if(is.null(dim(y_hat))){
    y_hat = matrix(y_hat, ncol=1)
  }else if(dim(y_hat)[1]==1){
    y_hat = matrix(y_hat, ncol=1)
  }
  # sanity check
  if(ncol(y_hat) != ncol(df)){
    stop('the number of fits does not match the number of df')
  }
  if(method=='cp' & is.null(sigma)){
    stop("need to specify sigma for Mallow's Cp")
  }

  # for AICc and BICc df larger than n-2 will cause trouble, round it
  if(method=='aicc' | method=='bicc'){
    df[which(df>=n-2)]=n-3
  }

  n <- nrow(y)
  nfit <- ncol(y_hat)
  if(nfit>1){
    y = matrix(rep(y,each=nfit),ncol=nfit,byrow=TRUE)
  }
  fit_wellness = Matrix::colSums((y - y_hat)^2)

  if(method=='aic'){return(log(fit_wellness/n)  + 2*df/n)}
  else if(method=='bic'){return(log(fit_wellness/n)  + log(n)*df/n)}
  else if(method=='aicc'){return(log(fit_wellness/n) + 2*(df+1)/(n-df-2))}
  else if(method=='bicc'){return(log(fit_wellness/n) + log(n)*(df+1)/(n-df-2))}
  else if(method=='gcv'){return(fit_wellness / (n-df)^2)}
  else if(method=='cp'){return(fit_wellness + 2*sigma^2*df)}
}
