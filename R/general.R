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

score <- function(ytrue, mu, s2, mult = 1){
  
  mF <- rep(mu, mult)
  s2F <- rep(s2, mult)
  seF <- (ytrue - mF)^2
  sc <- - seF/s2F - log(s2F)
  score <- mean(sc)
  
  return(score)
}

rmse <- function(ytrue, mu, mult = 1){
  mF <- rep(mu, mult)
  seF <- (ytrue - mF)^2
  rmse <- sqrt(mean(seF))
  
  return(rmse)
}

# Noise simulator

rn <- function(Xlam, mean = 0, theta.lam = 0.2, tau2.lam = 1, 
               px = c(0.2, 0.8), py = c(-2, 2)){
  
  if(!is.matrix(Xlam)) Xlam <- as.matrix(Xlam)
  px <- matrix(px, ncol = 1)
  
  condMu <- Exp2Sep(x1= Xlam ,x2= px, 1, theta =theta.lam, 1e-8) %*%
    solve(Exp2Sep(x1 = px, x2 = px, tau2 = 1, theta =theta.lam, g = 1e-8)) %*% py
  
  condSig <- Exp2Sep(x1=Xlam, x2 = Xlam, tau2 = 1, theta =theta.lam, g = 1e-8) -
    Exp2Sep(x1 = Xlam, x2= px, tau2 = 1, theta =theta.lam, g = 1e-8) %*%
    solve(Exp2Sep(x1=px, x2 = px, tau2 = 1, theta =theta.lam, g = 1e-8)) %*%
    Exp2Sep(x1 = px, x2 = Xlam, tau2 = 1, theta =theta.lam, g = 1e-8)
  
  llam <- sqrt(tau2.lam) * (drop(rmvnorm(1, mean = condMu, sigma= condSig)) + mean)
}

