// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;

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

// update Qr and v

// [[Rcpp::export]]
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

*/
