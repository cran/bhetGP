
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
NumericVector diag_quad_mat(NumericMatrix A, NumericMatrix B) {
  // Returns diag(A %*% B %*% t(A)) for s2 calculations
  int Arow = A.nrow();
  int Brow = B.nrow();
  NumericVector s2(Arow);
  for (int i = 0; i < Arow; i++) {
    s2(i) = 0.0;
    for (int j = 0; j < Brow; j++) {
      double temp_sum = 0.0;
      for (int n = 0; n < Brow; n++)
        temp_sum += A(i, n) * B(n, j);
      s2(i) += temp_sum * A(i, j);
    }
  }
  return s2;
}

// [[Rcpp::export]]
arma::mat Exp2(arma::mat distmat, const double tau2, const double theta,
                  const double g) { 
  // distmat = matrix of SQUARED distances
  int n1 = distmat.n_rows;
  int n2 = distmat.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      r = distmat(i, j) / theta;
      covmat(i, j) = tau2 * exp(-r);
    }
  }
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g; 
    //covmat(i, i) += tau2 * g * Ainv?
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat Exp2vec(arma::mat distmat, const double tau2, const double theta,
               const arma::vec g) { 
  // distmat = matrix of SQUARED distances
  int n1 = distmat.n_rows;
  int n2 = distmat.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      r = distmat(i, j) / theta;
      covmat(i, j) = tau2 * exp(-r);
    }
  }
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g(i); 
    //covmat(i, i) += tau2 * g * Ainv?
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat Exp2Sep(arma::mat x1, arma::mat x2, const double tau2, 
                     const arma::vec theta, const double g) { 
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (theta.n_elem != x1.n_cols)
    stop("length of theta does not match dimension of x");
  arma::mat covmat(n1, n2);
  double r;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      r = 0.0;
      for (int k = 0; k < d; k++)
        r += (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
      covmat(i, j) = tau2 * exp(-r);
    }
  }
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g;
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat Exp2SepVec(arma::mat x1, arma::mat x2, const double tau2, 
                  const arma::vec theta, const arma::vec g) { 
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (theta.n_elem != x1.n_cols)
    stop("length of theta does not match dimension of x");
  arma::mat covmat(n1, n2);
  double r;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      r = 0.0;
      for (int k = 0; k < d; k++)
        r += (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
      covmat(i, j) = tau2 * exp(-r);
    }
  }
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g(i);
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat Matern(arma::mat distmat, const double tau2, const double theta,
                    const double g, const double v) { 
  // distmat = matrix of SQUARED distances
  int n1 = distmat.n_rows;
  int n2 = distmat.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  if (v == 0.5) { 
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = sqrt(distmat(i, j) / theta);
        covmat(i, j) = tau2 * exp(-r);
      }
    }
  } else if (v == 1.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = sqrt(3 * distmat(i, j) / theta);
        covmat(i, j) = tau2 * (1 + r) * exp(-r);
      }
    }
  } else if (v == 2.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = sqrt(5 * distmat(i, j) / theta);
        covmat(i, j) = tau2 * (1 + r + pow(r, 2) / 3) * exp(-r);
      }
    }
  } 
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g;
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat MaternSep(arma::mat x1, arma::mat x2, const double tau2, const arma::vec theta,
                 const double g, const double v) { 
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (theta.n_elem != x1.n_cols)
    stop("length of theta does not match dimension of x");
  arma::mat covmat(n1, n2);
  double r;
  if (v == 0.5) { 
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        for (int k = 0; k < d; k++)
          r += (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
        covmat(i, j) = tau2 * exp(-sqrt(r));
      }
    }
  } else if (v == 1.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        for (int k = 0; k < d; k++)
          r += 3 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
        covmat(i, j) = tau2 * (1 + sqrt(r)) * exp(-sqrt(r));
      }
    }
  } else if (v == 2.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        for (int k = 0; k < d; k++)
          r += 5 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
        covmat(i, j) = tau2 * (1 + sqrt(r) + r / 3) * exp(-sqrt(r));
      }
    }
  } 
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g;
  }
  return covmat;
}
  
// [[Rcpp::export]]
arma::mat MaternVec(arma::mat distmat, const double tau2, const double theta,
                   const arma::vec g, const double v) { 
  // distmat = matrix of SQUARED distances
  int n1 = distmat.n_rows;
  int n2 = distmat.n_cols;
  arma::mat covmat(n1, n2);
  double r;
  if (v == 0.5) { 
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = sqrt(distmat(i, j) / theta);
        covmat(i, j) = tau2 * exp(-r);
      }
    }
  }else if (v == 1.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = sqrt(3 * distmat(i, j) / theta);
        covmat(i, j) = tau2 * (1 + r) * exp(-r);
      }
    }
  }else if (v == 2.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = sqrt(5 * distmat(i, j) / theta);
        covmat(i, j) = tau2 * (1 + r + pow(r, 2) / 3) * exp(-r);
      }
    }
  } 
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g(i);
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat MaternSepVec(arma::mat x1, arma::mat x2, const double tau2, const arma::vec theta,
                    const arma::vec g, const double v) { 
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (theta.n_elem != x1.n_cols)
    stop("length of theta does not match dimension of x");
  arma::mat covmat(n1, n2);
  double r;
  if (v == 0.5) { 
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        for (int k = 0; k < d; k++)
          r += (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
        covmat(i, j) = tau2 * exp(-sqrt(r));
      }
    }
  } else if (v == 1.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        for (int k = 0; k < d; k++)
          r += 3 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
        covmat(i, j) = tau2 * (1 + sqrt(r)) * exp(-sqrt(r));
      }
    }
  } else if (v == 2.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        for (int k = 0; k < d; k++)
          r += 5 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
        covmat(i, j) = tau2 * (1 + sqrt(r) + r / 3) * exp(-sqrt(r));
      }
    }
  } 
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g(i);
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat MaternProdSep(arma::mat x1, arma::mat x2, const double tau2, const arma::vec theta,
                        const double g, const double v) { 
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  double tmp;
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (theta.n_elem != x1.n_cols)
    stop("length of theta does not match dimension of x");
  arma::mat covmat(n1, n2);
  double r;
  double rexp;
  if (v == 1.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        rexp = 0.0;
        tmp = 1.0;
        for (int k = 0; k < d; k++){
          r = 3 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
          rexp += sqrt(r);
          tmp *= 1 + sqrt(r);
        }
        covmat(i, j) = tau2 * (tmp) * exp(-rexp);
      }
    }
  } else if (v == 2.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        rexp = 0.0;
        tmp = 1.0;
        for (int k = 0; k < d; k++){
          r = 5 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k); 
          rexp += sqrt(r);           
          tmp *= (1 + sqrt(r) + r/3);
        }
        covmat(i, j) = tau2 * (tmp) * exp(-rexp);
      }
    }
  } 
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g;
  }
  return covmat;
}

// [[Rcpp::export]]
arma::mat MaternProdSepVec(arma::mat x1, arma::mat x2, const double tau2, const arma::vec theta,
                           const arma::vec g, const double v) { 
  int n1 = x1.n_rows;
  int n2 = x2.n_rows;
  int d = x1.n_cols;
  if (x1.n_cols != x2.n_cols) 
    stop("dimension of x1 and x2 do not match");
  if (theta.n_elem != x1.n_cols)
    stop("length of theta does not match dimension of x");
  arma::mat covmat(n1, n2);
  double r;
  double rexp;
  double tmp;
  if (v == 1.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        rexp = 0.0;
        tmp = 1;
        for (int k = 0; k < d; k++){
          r = 3 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k);
          rexp += sqrt(r);
          tmp *= 1 + sqrt(r);
        }
        covmat(i, j) = tau2 * (tmp) * exp(-rexp);
      }
    }
  } else if (v == 2.5) {
    for (int i = 0; i < n1; i++) {
      for (int j = 0; j < n2; j++) {
        r = 0.0;
        rexp = 0.0;
        tmp = 1;
        for (int k = 0; k < d; k++){
          r = 5 * (x1(i, k) - x2(j, k)) * (x1(i, k) - x2(j, k)) / theta(k); 
          rexp += sqrt(r);           
          tmp *= (1 + sqrt(r) + r/3);
        }
        covmat(i, j) = tau2 * (tmp) * exp(-rexp);   
      }
    }
  } 
  if (n1 == n2) {
    for (int i = 0; i < n1; i++) 
      covmat(i, i) += tau2 * g(i);
  }
  return covmat;
}


