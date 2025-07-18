
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

/*
 * Code derived from GPvecchia package (Katzfuss et al.)
 */

#include "cov.h"
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
NumericVector forward_solve_raw(NumericMatrix U, NumericVector z,
                            NumericMatrix NNarray) {
  // Solves U * y = z for y
  // Uses raw form of U (in create_U use raw_form = TRUE)
  
  int n = U.nrow();
  NumericVector y(n);
  int mp1 = NNarray.ncol(); // m plus 1

  y(0) = z(0) / U(0, 0);
  
  for (int i = 1; i < n; i++) {
    int B = min(i + 1, mp1);
    y(i) = z(i);
    for (int j = 1; j < B; j++) { // continue comes back up here. Ignores update
      if (NNarray(i, j) == 0) continue; // this avoids passing NA altogether & ignores 0
      if (NumericVector::is_na(NNarray(i, j))) continue; //this ignores NA (in case passed)
      y(i) -= U(i, j) * y(NNarray(i, j) - 1); // U(i , j) = 0 if NN(i, j) = NA
      // if NNarray(i, j) was NA, you would just have y(i) = y(i) - 0 --> y(i)
    } 
    y(i) = y(i) / U(i, 0); //scaling to OG
  }
  return y;
}

// [[Rcpp::export]]
arma::mat rev_matrix(arma::mat x) {
  return reverse(x, 1);
}

double d2_vector(arma::rowvec x1, arma::rowvec x2) { 
  int n = x1.size();
  double d2 = 0.0;
  for(int k = 0; k < n; k++) {
    d2 += (x1[k] - x2[k]) * (x1[k] - x2[k]);
  }
  return d2;
}

arma::mat d2_matrix(arma::mat x) {
  int outrows = x.n_rows;
  int outcols = x.n_rows;
  arma::mat d2(outrows, outcols);
  for (int i = 0 ; i < outrows ; i++) {
    for (int j = 0 ; j < outcols ; j++) {
      d2(i, j) = d2_vector(x.row(i), x.row(j));
    }
  }
  return d2;
}

// [[Rcpp::export]]
arma::mat U_entries_Hom (const int Ncores, const arma::mat& x, const arma::umat& revNNarray,
                     const double tau2, const double theta, const double g, const double v){
  
  const int m = revNNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  arma::mat covmat;
  
  #ifdef _OPENMP
    
  #pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
    for (int k = 0; k < n; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds > 0)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat dist = d2_matrix(x.rows(inds00));
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2(dist, tau2, theta, g);
      } else {
        covmat = Matern(dist, tau2, theta, g, v);
      }
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(k, span(0, n0 - 1)) = M.t();
    }
    
  #else
    
    for (int k = 0; k < n; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds > 0)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat dist = d2_matrix(x.rows(inds00));
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2(dist, tau2, theta, g);
      } else {
        covmat = Matern(dist, tau2, theta, g, v);
      }
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(k, span(0, n0 - 1)) = M.t();
    }
    
  #endif
  
  return Lentries;
}

// [[Rcpp::export]]
arma::mat U_entries(const int Ncores, const arma::mat& x, const arma::umat& revNNarray,
                     const double tau2, const double theta, const arma::vec g, const double v){
  
  const int m = revNNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  arma::mat covmat;
  
#ifdef _OPENMP

// # pragma omp parallel for num_threads(threads) shared(Lentries) schedule(static)
# pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
  for (int k = 0; k < n; k++) {
    arma::uvec inds = revNNarray.row(k).t();
    arma::uvec inds00 = inds.elem(find(inds > 0)) - 1;
    uword n0 = inds00.n_elem;
    arma::mat dist = d2_matrix(x.rows(inds00));
    arma::mat covmat(n0, n0);
    arma::vec g00 = g(inds00);
    if (v == 999) {
      covmat = Exp2vec(dist, tau2, theta, g00); //think i need to pass only the specific g
    }else {
      covmat = MaternVec(dist, tau2, theta, g00, v);
    }
    arma::vec onevec = zeros(n0);
    onevec[n0 - 1] = 1;
    arma::vec M = solve(chol(covmat, "upper"), onevec);
    Lentries(k, span(0, n0 - 1)) = M.t();
  }
  
#else
  
  for (int k = 0; k < n; k++) {
    arma::uvec inds = revNNarray.row(k).t();
    arma::uvec inds00 = inds.elem(find(inds > 0)) - 1;
    uword n0 = inds00.n_elem;
    arma::mat dist = d2_matrix(x.rows(inds00));
    arma::mat covmat(n0, n0);
    arma::vec g00 = g(inds00);
    if (v == 999) {
      covmat = Exp2vec(dist, tau2, theta, g00); 
    } else {
      covmat = MaternVec(dist, tau2, theta, g00, v);
    }
    arma::vec onevec = zeros(n0); // initialize a vector onevec
    onevec[n0 - 1] = 1; // last element = 1
    arma::vec M = solve(chol(covmat, "upper"), onevec); // cholesky decomposition of covmat
    // solve system of equations // diagonals get used
    Lentries(k, span(0, n0 - 1)) = M.t(); //stores solution?
  }
  
#endif
  
  return Lentries;
}

// [[Rcpp::export]]
arma::mat U_entries_sep_Hom (const int Ncores, const arma::mat& x, const arma::umat& revNNarray, 
                     const double tau2, const arma::vec theta, const double g, const double v){
  
  const int m = revNNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  arma::mat covmat;
  
  #ifdef _OPENMP
    
  #pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
    for (int k = 0; k < n; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds > 0)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2Sep(x.rows(inds00), x.rows(inds00), tau2, theta, g);
      } else if( v > 1000){
        covmat = MaternProdSep(x.rows(inds00), x.rows(inds00), tau2, theta, g, (v-1000));
      } else {
        covmat = MaternSep(x.rows(inds00), x.rows(inds00), tau2, theta, g, v);
      }
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(k, span(0, n0 - 1)) = M.t();
    }
    
  #else
    
    for (int k = 0; k < n; k++) {
      arma::uvec inds = revNNarray.row(k).t();
      arma::uvec inds00 = inds.elem(find(inds > 0)) - 1;
      uword n0 = inds00.n_elem;
      arma::mat covmat(n0, n0);
      if (v == 999) {
        covmat = Exp2Sep(x.rows(inds00), x.rows(inds00), tau2, theta, g);
      } else if( v > 1000){
        covmat = MaternProdSep(x.rows(inds00), x.rows(inds00), tau2, theta, g, (v-1000));
      } else {
        covmat = MaternSep(x.rows(inds00), x.rows(inds00), tau2, theta, g, v);
      }
      arma::vec onevec = zeros(n0);
      onevec[n0 - 1] = 1;
      arma::vec M = solve(chol(covmat, "upper"), onevec);
      Lentries(k, span(0, n0 - 1)) = M.t();
    }
    
  #endif
  
  return Lentries;
}

// [[Rcpp::export]]
arma::mat U_entries_sep (const int Ncores, const arma::mat& x, const arma::umat& revNNarray, 
                         const double tau2, const arma::vec theta, const arma::vec g, const double v){
  
  const int m = revNNarray.n_cols - 1;
  const int n = x.n_rows;
  arma::mat Lentries = zeros(n, m + 1);
  arma::mat covmat;
  
#ifdef _OPENMP
  
#pragma omp parallel for num_threads(Ncores) shared(Lentries) schedule(static)
  for (int k = 0; k < n; k++) {
    arma::uvec inds = revNNarray.row(k).t();
    arma::uvec inds00 = inds.elem(find(inds > 0)) - 1;
    uword n0 = inds00.n_elem;
    arma::mat covmat(n0, n0);
    if (v == 999) {
      covmat = Exp2SepVec(x.rows(inds00), x.rows(inds00), tau2, theta, g(inds00));
    } else if( v > 1000){
      covmat = MaternProdSepVec(x.rows(inds00), x.rows(inds00), tau2, theta, g(inds00), (v-1000));
    } else {
      covmat = MaternSepVec(x.rows(inds00), x.rows(inds00), tau2, theta, g(inds00), v);
    }
    arma::vec onevec = zeros(n0);
    onevec[n0 - 1] = 1;
    arma::vec M = solve(chol(covmat, "upper"), onevec);
    Lentries(k, span(0, n0 - 1)) = M.t();
  }
  
#else
  
  for (int k = 0; k < n; k++) {
    arma::uvec inds = revNNarray.row(k).t();
    arma::uvec inds00 = inds.elem(find(inds > 0)) - 1;
    uword n0 = inds00.n_elem;
    arma::mat covmat(n0, n0);
    if (v == 999) {
      covmat = Exp2SepVec(x.rows(inds00), x.rows(inds00), tau2, theta, g(inds00));
    }  else if( v > 1000){
      covmat = MaternProdSepVec(x.rows(inds00), x.rows(inds00), tau2, theta, g(inds00), (v-1000));
    } else {
      covmat = MaternSepVec(x.rows(inds00), x.rows(inds00), tau2, theta, g(inds00), v);
    }
    arma::vec onevec = zeros(n0);
    onevec[n0 - 1] = 1;
    arma::vec M = solve(chol(covmat, "upper"), onevec);
    Lentries(k, span(0, n0 - 1)) = M.t();
  }
  
#endif
  
  return Lentries;
}

// [[Rcpp::export]]
void check_omp () {
  #ifdef _OPENMP
    // DO NOTHING
  #else 
    Rcout << "NOTE: OpenMP install suggested for best results; see ?bhetGP for details \n";
  #endif
}

// [[Rcpp::export]]
int check_cores (const int Ncores) {
int threads = Ncores;
#ifdef _OPENMP
   int max_threads = omp_get_max_threads(); // this gives max cores on comp
   if (threads > max_threads) {
   // Rcpp::Rcout << "Adjusting Ncores: requested " << threads
   //            << ", using max available (" << limit_threads << ").\n"; //sanity
   threads = max_threads;
   }
#else 
  threads = 0;
#endif
return(threads);
}
  
// [[Rcpp::export]]
arma::mat row_col_pointers(const arma::umat& NNarray) {
    
  const int m = NNarray.n_cols- 1; // removing one col. dim
  const int n = NNarray.n_rows; 
  int start, col_count;
    
  int length = (n - m) * (m + 1);
  for (int i = 1; i <= m; i ++)
    length += i;
    
  arma::mat pointers = zeros(length, 2); // 2 columns only
  
  // To get rid of NAs: first pass them as zeros; 
  // Since they are always in first m columns starting from last index, they 
  // are never accessed because i <= m ---> for first m, only does first i elements
  
  start = 0;
  for (int i = 1; i <= n; i++) { // starts from index 1
    if (i <= m) { // accesses only i of m elements
      col_count = i - 1; // how many elements will the row have; gives rhs start point
      for (int j = start; j < start + i; j++) { // start with 0 based index
        pointers(j, 0) = i; // row index = column first element
        pointers(j, 1) = NNarray(i - 1, col_count); // From NN array use (i - 1), col neighbor
        col_count -= 1; // move left of columns
      }
      start += i; // move pointer
    } else { // reached max cols
      col_count = m;
      for (int j = start; j < start + m + 1; j++) {
        pointers(j, 0) = i;
        pointers(j, 1) = NNarray(i - 1, col_count);
        col_count -= 1;
      }
      start += m + 1; // move by by m + 1 positions (next row start)
    }
  }
  return pointers; // returns row <-- neighbors starting from RHS
}

