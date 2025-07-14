# Function Contents (Annie) ---------------------------------------------------
#   invdet: calculates inverse and log determinant of a matrix using C
#           (credit given to the "laGP package, R.B. Gramacy & Furong Sun)
#   sq_dist
#   score
#   ifel: returns second or third element depending on first argumentifel

eps <- sqrt(.Machine$double.eps)

# Matrix Inverse --------------------------------------------------------------

# this is to use C function from within R
invdet <- function(M) {

  n <- nrow(M)
  out <- .C("inv_det_R",
            n = as.integer(n),
            M = as.double(M),
            Mi = as.double(diag(n)),
            ldet = double(1),
            PACKAGE = "bhetGP")

  return(list(Mi = matrix(out$Mi, ncol=n), ldet = out$ldet))
}


sq_dist <- function(X1, X2 = NULL) {

  X1 <- as.matrix(X1)
  n1 <- nrow(X1)
  m <- ncol(X1)

  if(is.null(X2)) {
    outD <- .C("distance_symm_R",
               X = as.double(t(X1)),
               n = as.integer(n1),
               m = as.integer(m),
               D = double(n1 * n1),
               PACKAGE = "bhetGP")
    return(matrix(outD$D, ncol = n1, byrow = TRUE))
  } else {
    X2 <- as.matrix(X2)
    n2 <- nrow(X2)
    if(ncol(X1) != ncol(X2)) stop("dimension of X1 & X2 do not match")
    outD <- .C("distance_R",
               X1 = as.double(t(X1)),
               n1 = as.integer(n1),
               X2 = as.double(t(X2)),
               n2 = as.integer(n2),
               m = as.integer(m),
               D = double(n1 * n2),
               PACKAGE = "bhetGP")
    return(matrix(outD$D, ncol = n2, byrow = TRUE))
  }
}

crps <- function(y, mu, s2) {
  sigma <- sqrt(s2)
  z <- (y - mu) / sigma
  return(mean(sigma * (-(1 / sqrt(pi)) + 2 * dnorm(z) + z * (2 * pnorm(z) - 1))))
}

# If else ---------------------------------------------------------------------

ifel <- function(logical, yes, no) {
  if (logical) {
    return(yes)
  } else return(no)
}

## Mappings - use with vdims

# map <- function(reps_vdims, y){
#   map <- list()
#   map$xv <- reps_vdims$X0
#   map$yv <- sapply(reps_vdims$Zlist, function(i) mean(y[i]))
#   map$yvs2 <- sapply(reps_vdims$Zlist, function(i) sum((y[i] - mean(y[i]))^2))
#   map$orderF <- reps_vdims$Z
#   map$multF <- reps_vdims$mult
#   return(map)
# }

