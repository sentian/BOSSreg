% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boss.R
\name{coef.boss}
\alias{coef.boss}
\title{Select coefficient vector(s) for BOSS.}
\usage{
\method{coef}{boss}(object, ic = c("aicc", "bicc", "aic", "bic", "gcv",
  "cp"), select.boss = NULL, ...)
}
\arguments{
\item{object}{The boss object, returned from calling the \code{boss} function.}

\item{ic}{Which information criterion is used to select the optimal coefficient vector for BOSS.
The default is AICc-hdf.}

\item{select.boss}{The index (or indicies) of columns in the coefficient matrix for which
one wants to select. By default (NULL) it's selected by the information criterion specified in
'ic'.}

\item{...}{Extra arguments (unused for now)}
}
\value{
The chosen coefficient vector(s) for BOSS.
}
\description{
This function returns the optimal coefficient vector of BOSS selected by AICc
(by default) or other types of information criterion.
}
\details{
If \code{select.boss} is specified, the function returns
corresponding column(s) in the coefficient matrix.

If \code{select.boss} is unspecified, the function returns the optimal coefficient
vector selected by AICc-hdf (other choice of IC can be specified in \code{ic}).
The only exception is when n>=p, where hdf is not well defined, and the entire coefficient matrix
is returned.
}
\examples{
# See the example in the section of \\code{boss}. Or type ?boss in R.

}