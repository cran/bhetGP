
# Function Contents (Annie) ---------------------------------------------------
# Internal:
#   rand_mvn_vec: samples from MVN Gaussian using vecchia approximation
#   create_U: uses C++ to create sparse Matrix U
#   create_approx: creates vecchia approximation object
#   update_obs_in_approx: updates x_ord inside approx object
#   add_pred_to_approx: incorporates predictive locations inside vecchia approx

# Random MVN ------------------------------------------------------------------

rand_mvn_vec <- function(approx, theta, v, mean = 0,
                         g = eps, scale = 1, sep) {
  
  # unique x_ord
  z <- rnorm(nrow(approx$x_ord))
  # Add sep = T for seperable inner layer?
  
  # create U based on x_approx
  Uraw <- create_U(approx, g, theta, v, sep = sep, raw_form = TRUE)/sqrt(scale)
  
  #forward solve to obtain lam prior
  sample <- forward_solve_raw(Uraw, z, approx$NNarray)
  
  # arrange lam prior for original x order (not x_approx ord)
  sample <- sample[approx$rev_ord_obs] + mean
  return(sample)
}

# Create U --------------------------------------------------------------------

create_U <- function(approx, g, theta, v, sep = FALSE,
                     raw_form = FALSE) {
  
  revNN <- approx$revNN
  revNN[is.na(revNN)] <- 0
  
  n <- nrow(approx$x_ord)
  if(length(g) == 1){
    if (sep) {
      U <- U_entries_sep_Hom(approx$n_cores, approx$x_ord, revNN, 
                         1, theta, g, v)
    } else U <- U_entries_Hom(approx$n_cores, approx$x_ord, revNN, 
                          1, theta, g, v)
  }else{
    if (sep) {
      U <- U_entries_sep(approx$n_cores, approx$x_ord, revNN, 
                         1, theta, g, v)
    } else U <- U_entries(approx$n_cores, approx$x_ord, revNN, 
                          1, theta, g, v)
    
  }
  
  if (raw_form) {
    m <- ncol(approx$NNarray) - 1
    U <- rev_matrix(U)
    for (i in 1:m) {
      zeros <- (U[i, ] == 0)
      U[i, ] <- c(U[i, !zeros], rep(0, times = sum(zeros))) 
    }
    return(U)
  } else {
    U <- c(t(U))[approx$notNA]
    U <- Matrix::sparseMatrix(i = approx$pointers[, 2], 
                              j = approx$pointers[, 1], x = U, 
                              dims = c(n, n))
    return(U)
  } 
}


# Create Approximation --------------------------------------------------------

create_approx <- function(x, m, ordering = NULL, cores = NULL) {
  
  n <- nrow(x)
  if (is.null(ordering)) 
    ordering <- sample(1:n, n, replace = FALSE)
  rev_ord_obs <- order(ordering)
  x_ord <- x[ordering, , drop = FALSE]
 
  NNarray <- GpGp::find_ordered_nn(x_ord, m)
  revNN <- rev_matrix(NNarray)
  notNA <- as.vector(t(!is.na(NNarray)))
  pointers <- row_col_pointers(NNarray)

  if(is.null(cores))
    n_cores <- min(parallel::detectCores(all.tests = FALSE, logical = TRUE), 2)
  else 
    n_cores <- cores
  
  out <- list(m = m, ord = ordering, NNarray = NNarray, revNN = revNN, 
              notNA = notNA, pointers = pointers,
              n_cores = n_cores, rev_ord_obs = rev_ord_obs, x_ord = x_ord)
  return(out)
}


# Update Observation in Approx ------------------------------------------------

update_obs_in_approx <- function(approx, x_new, col_index = NULL) {
  
  if (is.null(col_index)) {
    approx$x_ord <- x_new[approx$ord, , drop = FALSE]
  } else {
    approx$x_ord[, col_index] <- x_new[approx$ord]
  }
  
  return(approx)
}


# Add Pred to Approx ----------------------------------------------------------

# --- combines approx$x_ord with x_pred (which is also ordered in a new way)
# --- gets nn for new X_comb in same format as "approx"

add_pred_to_approx <- function(approx, x_pred, m, ordering_new = NULL) {
  
  n <- nrow(approx$x_ord)
  n_pred <- nrow(x_pred)
  if (is.null(ordering_new)) 
    ordering_new <- sample(1:n_pred, n_pred, replace = FALSE)
  rev_ord_pred <- order(ordering_new)
  
  ord <- c(approx$ord, ordering_new + n) # observed data FIRST
  x_ord <- rbind(approx$x_ord, x_pred[ordering_new, , drop = FALSE])
  observed <- c(rep(TRUE, n), rep(FALSE, n_pred))
  
  NNarray <- GpGp::find_ordered_nn(x_ord, m)
  revNN <- rev_matrix(NNarray)
  notNA <- as.vector(t(!is.na(NNarray)))
  pointers <- row_col_pointers(NNarray)

  out <- list(m = m, ord = ord, NNarray = NNarray, revNN = revNN, 
              notNA = notNA, pointers = pointers,
              n_cores = approx$n_cores, rev_ord_obs = approx$rev_ord_obs,
              rev_ord_pred = rev_ord_pred, observed = observed, x_ord = x_ord)
  return(out)
}

