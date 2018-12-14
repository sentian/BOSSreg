// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

// givens rotation
void givens(double a, double b, double &c, double &s){
  if(b == 0){
    c = 1;
    s = 0;
  }else{
    if(fabs(b) >= fabs(a)){
      double t = -a / b;
      s = 1 / sqrt(1 + pow(t,2));
      c = s * t;
    }else{
      double t = -b / a;
      c = 1 / sqrt(1 + pow(t,2));
      s = c * t;
    }
  }
}

// update Qr and v
void updateQv(arma::vec &v, arma::mat &Qr){
  int r = Qr.n_cols;
  double c, s;
  for (int i=r-1; i>=1; i--){
    givens(v(i-1), v(i), c, s);
    // update v
    v(i-1) = c * v(i-1) - s * v(i);
    // update Q
    arma::mat rotate = { {c, s},
                         {-s, c} };
    Qr.cols(i-1,i) =  Qr.cols(i-1,i) * rotate;
  }
}

// update QR decomposition
void updateQR(arma::mat &Ql, arma::mat &Qr, arma::mat &R, arma::vec newcol){
  int n = Ql.n_rows;
  int l = Ql.n_cols;
  arma::mat Q = join_rows(Ql,Qr);
  arma::vec v = Q.t() * newcol;
  arma::vec v_toupdate = v.subvec(l,n-1);
  updateQv(v_toupdate, Qr);
  // update v
  v.subvec(l,n-1) = v_toupdate;
  // update Q
  Ql = join_rows(Ql, Qr.col(0));
  Qr.shed_col(0);
  // update R
  R = join_cols(R, arma::zeros(l).t());
  R = join_rows(R, v.subvec(0,l));
}

// guided QR decomposition, based on semi-partial correlation
// [[Rcpp::export]]
List guideQR(arma::mat x, arma::vec y, int maxstep){
  // int n = x.n_rows;
  int p = x.n_cols;
  // which_remain stores the indices of remaining preditors
  arma::vec which_remain = arma::vec(maxstep);
  std::iota (std::begin(which_remain), std::end(which_remain), 0);

  arma::mat Q, R, x_remain = x;
  arma::vec steps = arma::vec(maxstep);

  // first variable to step in
  arma::vec cor_tmp = arma::abs(x.t() * y);
  arma::uword i = cor_tmp.index_max();
  steps(0) = i;
  which_remain.shed_row(i);
  x_remain.shed_col(i);
  arma::qr(Q, R, x.col(i));
  arma::mat Ql = Q, Qr = Q;
  Ql = Ql.col(0);
  Qr.shed_col(0);
  R = R.row(0);

  arma::mat resid_tmp = x_remain;
  int j_remain, j_tmp;
  // step 2 and beyond
  //int i_step = 1;
  for (int i_step=1; i_step<=maxstep-1; i_step++){
    // calculate the semi-partial correlations for each remaining predictor
    if(i_step < p-1){
      j_tmp = Ql.n_cols - 1; // -1 to match index in C
      resid_tmp = resid_tmp - Ql.col(j_tmp) * Ql.col(j_tmp).t() * x_remain;
      j_remain = arma::abs(y.t() * arma::normalise(resid_tmp)).index_max();
    }else{
      j_remain = 0;
    }

    steps(i_step) = which_remain(j_remain);
    // update QR
    updateQR(Ql, Qr, R, x_remain.col(j_remain));

    // update others
    if(i_step < p-1){
      resid_tmp.shed_col(j_remain);
      x_remain.shed_col(j_remain);
      which_remain.shed_row(j_remain);
    }
  }

  // add 1 to match the index in R
  steps = steps + 1;
  return List::create(Named("Q") = Ql,
                      Named("R") = R,
                      Named("steps") = steps);

}


/* // Old versions

// givens rotation
void givens(double a, double b, double *c, double *s){
  if(b == 0){
    *c = 1;
    *s = 0;
  }else{
    if(fabs(b) >= fabs(a)){
    double t = -a / b;
    *s = 1 / sqrt(1 + pow(t,2));
    *c = *s * t;
    }else{
      double t = -b / a;
      *c = 1 / sqrt(1 + pow(t,2));
      *s = *c * t;
    }
  }
}

List updateQv(arma::vec v, arma::mat Qr){
    int r = Qr.n_cols;
    double c, s;
    for (int i=r-1; i>=1; i--){
        givens(v[i-1], v[i], &c, &s);
        // update v
        v[i-1] = c * v[i-1] - s * v[i];
        // update Q
        arma::mat rotate = { {c, s},
                       {-s, c} };
        Qr.cols(i-1,i) =  Qr.cols(i-1,i) * rotate;
    }
    return List::create(Named("v") = v,
                        Named("Qr") = Qr);
}
*/


/*
# corresponding R functions
# Givens rotation
# algorithm 1.1 in the reference http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.142.2571&rep=rep1&type=pdf
givens <- function(a, b){
  if(b == 0){
    c = 1
    s = 0
  }else{
    if(abs(b) >= abs(a)){
      t = -a/b
      s = 1/sqrt(1+t^2)
      c = s*t
    }else{
      t = -b/a
      c = 1/sqrt(1+t^2)
      s = c*t
    }
  }
  return(list(c=c,s=s))
}

# add a column to X_before, and update QR
# updateQv takkes only the minimum required input and do the update for efficiency
updateQv <- function(v, Qr){
  n = dim(Qr)[1]
  r = dim(Qr)[2]

  for(i in r:2){
    givens_result = givens(v[i-1], v[i])
    c =  givens_result$c
    s =  givens_result$s
    # update v
    v[i-1] = c * v[i-1] - s * v[i]
    # update Q
    Qr[,(i-1):i] =  Qr[,(i-1):i] %*% cbind(c(c,-s),c(s,c))
  }
  return(list(v=v, Qr=Qr))
}

## update the QR decomposition when a column is added to X -------------------------
# here's a reference for the algorithms
# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.142.2571&rep=rep1&type=pdf

# X_before = Q_before R_before, X_before1:n*l, Q_before:n*n, R_before:n*l
# Q_before = [Ql_before Qr_before], where Ql_before: n*l and Qr_before: n*(n-l)
# add a column u to X_active_update: X_after=[X_before u], and update Q_before and R_before, get X_after=Q_after R_after
# X_after: n*(r+1), Q_after = [Ql_after Qr_after] where Ql_after = [Ql_before Ql_newcol]
# algorithm 2.19 in the reference
updateQR <- function(Ql, Qr, R, newcol){
  n = dim(Ql)[1]
  l = dim(Ql)[2]
  v = t(cbind(Ql, Qr)) %*% newcol
  updateQv_result = updateQv(v[(l+1):n], Qr)
  # update v
  v[(l+1):n] = updateQv_result$v
  # update Q
  Ql = cbind(Ql, updateQv_result$Qr[, 1])
  Qr = updateQv_result$Qr[, -1]
  # update R
  R = rbind(R, rep(0, l))
  R = cbind(R, v[1:(l+1)] )
  return(list(Ql=Ql, Qr=Qr, R=R))
}


guideQR <- function(x,y,maxstep){
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
}
*/

