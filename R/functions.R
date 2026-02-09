# ----------- Log Likelihood Functions -----------------------------------------
# loglw: Woodbury log lik for separable cases
# logl: Complete log lik for separable cases
# loglw_iso: Woodbury log lik for isotropic cases
# logl_iso: Complete log lik for isotropic cases
# logl_vec: Complete vecchia log likelihood for sep and iso
# loglw_vec: woodbury vecchia log likelihood for sep and iso
# mvn_prior: Lambda prior for non vecchia cases
# ------------------------------------------------------------------------------

# woodbury llik; outer = TRUE for prof llik; 

loglw <- function(ys2, yn, xn, nugs, A, v, theta, outer = TRUE, calc_tau2 = TRUE, 
                  mean = 0, scale = 1, a , b){
  
  # Calculate the nugget vector: A_inv * Lam
  lam_adj <- nugs/A

  # Build Covariance (Vec = vector)
  if(v == 999){
    Kn <- Exp2SepVec(xn, xn, tau2 = scale, theta = theta, g = lam_adj)
  }else if(v > 1000){
    Kn <- MaternProdSepVec(xn, xn, tau2 = scale, theta = theta, g = lam_adj, v= (v- 1000))
  }
  else{
    Kn <- MaternSepVec(xn, xn, tau2 = scale, theta = theta, g = lam_adj, v= v)
  }
  
  # Kinv <- solve(Kn)
  Kinv <- invdet(Kn)
  # if(mean != 0) mean <- drop(rowSums(Kinv) %*% yn/(sum(Kinv))) # 1t Kn(-1) Delta (1t Kn(-1) 1) ^ (-1)
  # ys2 <- unlist(lapply(yc, function(x) sum((x - mean(x))^2)))
  yn <- yn - mean
  
  N <- sum(A)
  q <- sum(ys2/nugs)
  
  quadterm <-  q + t(yn) %*% Kinv$Mi %*% yn # solve(Kn,yn) # woodbury adjusted quad 
  
  logdet <- (-0.5) * Kinv$ldet # determinant(Kn, logarithm=TRUE)$modulus # log (det(Kn)) 
  
  # Adjused linear term
  S <- (-0.5) * sum((A-1) * log(nugs) + log(A))
  
  # outer = prof likelihood (always true for Y layer)
  if (outer) {
    logl <- logdet - ((N + a) * 0.5) * log(quadterm + b) + S
  } else {
    logl <- logdet - 0.5 * quadterm + S
  }
  
  # return tau2 if profll
  if (calc_tau2) {
    tau2 <- c(quadterm + b) / (N + a)
  } else tau2 <- NULL
  
  return(list(llik = logl, tau2 = tau2))
}

# complete llik (for lambda layer & full N version)
# outer = FALSE for lamba layer; TRUE for full Y version
logl <- function(y, x, nugs = NULL, theta, v, outer = TRUE, calc_tau2 = TRUE, 
                 mean = 0, scale = 1, a , b){
  
  if(v == 999){
    if(length(nugs) == 1) KN <- Exp2Sep(x, x, tau2 = scale, theta = theta, g = nugs) # homGP
    else KN <- Exp2SepVec(x, x, tau2 = scale, theta = theta, g = nugs) # hetGP
  }else if(v > 1000){
    if(length(nugs) == 1) KN <- MaternProdSep(x, x, tau2 = scale, theta = theta, g = nugs, v= (v - 1000))
    else KN <- MaternProdSepVec(x, x, tau2 = scale, theta = theta, g = nugs, v= (v - 1000))
  }
  else{
    if(length(nugs) == 1) KN <- MaternSep(x, x, tau2 = scale, theta = theta, g = nugs, v= v)
    else KN <- MaternSepVec(x, x, tau2 = scale, theta = theta, g = nugs, v= v)
  }
  
  # Kinv <- solve(KN)
  Kinv <- invdet(KN)
  N <- length(y)
  
  # if(mean != 0) mean <- drop(rowSums(Kinv) %*% y/(sum(Kinv)))
  y <- y - mean

  quadterm <- (t(y) %*% Kinv$Mi %*% y)
  
  # outer for Y layer, False for lambda layer
  if(outer) {
    llik <- (-0.5 * (N + a)) * log (quadterm + b) - 0.5 * Kinv$ldet
      # (0.5 * determinant(KN, logarithm=TRUE)$modulus)
  }else {
    llik <- (-0.5) * Kinv$ldet - (0.5) * (quadterm)
      # ((-0.5) * determinant(KN, logarithm=TRUE)$modulus) - (0.5) * (quadterm)
  } 
  
  if(calc_tau2){
    tau2 <- c(quadterm + b) / (N + a)
  } else tau2 <- NULL
  
  return(list(llik = llik, tau2 = tau2))
}

# Product Kernel does not matter for iso functions
# woodbury for anisotropic case
loglw_iso <- function(ys2, yn, dx_n, nugs, A, v, theta, outer = TRUE, calc_tau2 = TRUE, 
                      mean = 0, scale = 1, a , b){
  
  # Lam_N * A_inv
  lam_adj <- nugs/A
  
  # if(length(nugs) == 1) lamN <- rep(nugs, sum(A)) # homGP
  # else lamN <- rep(nugs, A) # hetGP
  
  if(v == 999){
    Kn <- Exp2vec(dx_n, tau2 = scale, theta = theta, g = lam_adj)
  }else{
    if(v > 1000) v = v - 1000
    Kn <- MaternVec(dx_n, tau2 = scale, theta = theta, g = lam_adj, v= v)
  }

  N <- sum(A)
  n <- length(yn)
  
  # Kinv <- solve(Kn)
  Kinv <- invdet(Kn)
  # if(mean != 0) mean <- drop(rowSums(Kinv) %*% yn/(sum(Kinv))) # 1t Kn(-1) Delta (1t Kn(-1) 1) ^ (-1)
  yn <- yn - mean
  q <- sum(ys2/nugs)
  
  quadterm <-  q + t(yn) %*% Kinv$Mi %*% yn#solve(Kn,yn) 
  
  logdet <- (-0.5) * Kinv$ldet # determinant(Kn, logarithm=TRUE)$modulus ## log (det(Kn)) 
  
  S <- (-0.5) * sum((A-1) * log(nugs) + log(A))
  
  if (outer) { # prof llik
    logl <- logdet - ((N + a) * 0.5) * log(quadterm + b) + S
  } else {
    logl <- logdet - 0.5 * quadterm + S
  }
  
  if(calc_tau2) { # if prof llik and wish to infer tau2
    tau2 <- c(quadterm + b) / (N + a)
  } else tau2 <- NULL
  
  return(list(llik = logl, tau2 = tau2))
}

# complete llik for anisotropic case
logl_iso <- function(y, dx, nugs = NULL, theta, v, outer = TRUE, calc_tau2 = TRUE,
                     mean = 0, scale = 1, a , b){

  if(v == 999){
    if(length(nugs) == 1) KN <- Exp2(dx, tau2 = scale, theta = theta, g = nugs) # homGP
    else KN <- Exp2vec(dx, tau2 = scale, theta = theta, g = nugs) # hetGP
  }else{
    if(v > 1000) v = v - 1000
    if(length(nugs) == 1) KN <- Matern(dx, tau2 = scale, theta = theta, g = nugs, v= v) # homGP
    else KN <- MaternVec(dx, tau2 = scale, theta = theta, g = nugs, v= v) # hetGP
  }

  # Kinv <- solve(KN)
  Kinv <- invdet(KN)
  # if(mean != 0) mean <- drop(rowSums(Kinv) %*% y/(sum(Kinv)))
  y <- y - mean
  
  N <- length(y)
  quadterm <- (t(y) %*% Kinv$Mi %*% y)
  
  # outer for Y layer, false for Lam layer
  if(outer) {
    llik <- (- 0.5 *(N + a)) * log (quadterm + b) - (0.5) * Kinv$ldet
      # (0.5 * determinant(KN, logarithm=TRUE)$modulus)
  }else{
    # llik <- ((-0.5) * determinant(KN, logarithm=TRUE)$modulus) - (0.5) * (quadterm)
    llik <- (-0.5) * Kinv$ldet - (0.5) * (quadterm)
  }
  
  if(calc_tau2) {
    tau2 <- c(quadterm + b) / (N + a)
  } else tau2 <- NULL
  
  return(list(llik = llik, tau2 = tau2))
}

# ------------------ Vecchia Likelihoods ---------------------------------------

# Vecchia complete logl (use for full N and lambda layer)
# Lam are unordered, Y is ordered 
logl_vec <- function(out_vec, approx, nugs, theta, outer = TRUE, v, calc_tau2 = TRUE,
                     sep = FALSE, mean, scale = 1, latent = FALSE, a , b) {
  
  # a <- b <- 1.5
  n <- nrow(approx$x_ord)
  out_vec <- out_vec - mean
  
  if(length(nugs) == 1){
    # homGP
    if(latent == TRUE){
      # For lambda layer, out_vec = lambda which are unordered
      out_vec_ord <- out_vec[approx$ord] # order lambdas
      U_mat <- create_U(approx, nugs, theta, v = v, sep = sep) / sqrt(scale) # does either isotropic or separable
    }else{
      out_vec_ord <- out_vec # out_vec = Y --> no need to order. 
      U_mat <- create_U(approx, nugs, theta, v = v, sep = sep) / sqrt(scale) # does either isotropic or separable
    }
  } else{ # hetGP
    G <- nugs[approx$ord] # Always order nuggets
    out_vec_ord <- out_vec # out_vec -> Y -> already ordered
    U_mat <- create_U(approx, G, theta, v = v, sep = sep) / sqrt(scale) # does either isotropic or separable
  }

  Uty <- Matrix::crossprod(U_mat, out_vec_ord) # Ut * Y
  
  ytUUty <- sum(Uty^2) # Y(t) UUt Y
  logdet <- sum(log(Matrix::diag(U_mat))) # sum log UUt
  
  if (outer) {
    # logl <- logdet - (n * 0.5) * log(ytUUty) # prof llik
    logl <- logdet - ((n + a) * 0.5) * log((ytUUty + b)) # prof llik w prior IG(a/2, b/2)
  } else {
    logl <- logdet - 0.5 * ytUUty
  }
  
  if(calc_tau2) { # currently using a jeffrey's prior
    tau2 <- (c(ytUUty) + b) / (n + a)
  } else tau2 <- NULL
  
  return(list(llik = logl, tau2 = tau2))
}


# Works for homoskedastic GP with replicates
# Take in pre ordered Y's
loglw_vec <- function(out_vec_s2, out_vec, approx, A, nugs, theta, outer = TRUE, 
                      v, calc_tau2 = TRUE, sep = FALSE, mean = 0, scale = 1, a , b) {
  
  N <- sum(A)
  n <- length(out_vec)

  out_vec <- out_vec - mean
  
  # order according to X's
  A_ord <- A[approx$ord]
  if(length(nugs) ==1) nugs <- rep(nugs, length(out_vec)) # homGP -> repeat same nugs
  g_ord <- nugs[approx$ord] # order lambdas
  G <- g_ord/A_ord # calculate adjusted Lam (n) or g(n) * A_inv 
  
  U_mat <- create_U(approx, G, theta, v, sep = sep) / sqrt(scale) # does either isotropic or separable
  Ln.inv <- 1/g_ord # Lam(n) ^ (-1)
  
  Uty <- Matrix::crossprod(U_mat, out_vec) 
  ytUUty <- sum(Uty^2)
  
  term <- sum(out_vec_s2 * Ln.inv)
  quadterm <- term + ytUUty
  
  logdet <- sum(log(Matrix::diag(U_mat)))
  
  # linear term
  S <- (-0.5) * sum((A_ord -1) * log(g_ord) + log(A_ord))
  
  # adjust prof likelihood
  if (outer) {
    # logl <- logdet - ((N * 0.5) * log(quadterm) + S
    logl <- logdet - ((N + a) * 0.5) * log((quadterm + b)) + S
  } else {
    logl <- logdet - 0.5 * quadterm + S
  }
  
  if(calc_tau2) {
    tau2 <- c(quadterm + b) / (N + a)
  } else tau2 <- NULL
  
  return(list(llik = logl, tau2 = tau2))
}

# ------------------ Lambda Prior ----------------------------------------------

# Draw a prior for lambda's
mvn_prior <- function(x, v, theta0, g0, dx = NULL, mean0, scale0 = 1, sep = FALSE){
  
  if(length(g0) == 1) g0 <- rep(g0, nrow(x))
  if(length(mean0) == 1) mean0 <- rep(mean0, nrow(x))
  
  if(v == 999){
    if(sep) {Kxlam <- Exp2SepVec(x, x, tau2 = scale0, theta = theta0, g = g0)
    } else Kxlam <- Exp2vec(dx, tau2 = scale0, theta = theta0, g = g0) 
  }else{
    if(sep) {
      if(v > 1000) Kxlam <- MaternProdSepVec(x, x, tau2= scale0, theta = theta0, g = g0, v = (v- 1000))
      else Kxlam <- MaternSepVec(x, x, tau2= scale0, theta = theta0, g = g0, v = v)
    }else{
      if(v > 1000) v = v - 1000
      Kxlam <- MaternVec(dx, tau2 = scale0, theta = theta0, g = g0, v = v)
    } 
  }
  
  llam_nv <- drop(rmvnorm(1, mean = mean0, sigma = Kxlam))
  return(llam_nv)
}


# --------- Old and working (check to match if needed) ------------------------

loglw_N <- function(y, yn, xn, nugs, A, v, theta, outer = TRUE, calc_tau2 = TRUE, 
                    mean = 0, scale = 1, a, b){
  
  # Calculate the nugget vector: A_inv * Lam
  lam_adj <- nugs/A
  
  # Create Lam(N) for hom (len(nugs) == 1) and het(!that)
  if(length(nugs) == 1) lamN <- rep(nugs, sum(A))
  else lamN <- rep(nugs, A) 
  
  # Build Covariance (Vec = vector)
  if(v == 999){
    Kn <- Exp2SepVec(xn, xn, tau2 = scale, theta = theta, g = lam_adj)
  }else if(v > 1000){
    Kn <- MaternProdSepVec(xn, xn, tau2 = scale, theta = theta, g = lam_adj, v= (v- 1000))
  }
  else{
    Kn <- MaternSepVec(xn, xn, tau2 = scale, theta = theta, g = lam_adj, v= v)
  }
  
  # Kinv <- solve(Kn)
  Kinv <- invdet(Kn)
  # if(mean != 0) mean <- drop(rowSums(Kinv) %*% yn/(sum(Kinv))) # 1t Kn(-1) Delta (1t Kn(-1) 1) ^ (-1)
  
  y <- y - mean
  yn <- yn - mean
  
  N <- length(y)
  q1 <- sum(y^2 /lamN) # Y(N)t * Lam(N) ^ (-1) * Y(N)
  q2 <- sum(yn^2 * A/nugs) # Y(n)t * Lam(n) ^ (-1) * Y(n)
  
  quadterm <-  q1 - q2 + t(yn) %*% Kinv$Mi %*% yn #solve(Kn, yn) # woodbury adjusted quad 
   
  logdet <- (-0.5) * Kinv$ldet #determinant(Kn, logarithm=TRUE)$modulus # log (det(Kn)) 
  
  # Adjused linear term
  S <- (-0.5) * sum((A-1) * log(nugs) + log(A))
  
  # outer = prof likelihood (always true for Y layer)
  if (outer) {
    logl <- logdet - (N * 0.5) * log(quadterm) + S
  } else {
    logl <- logdet - 0.5 * quadterm + S
  }
  
  # return tau2 if profll
  if(calc_tau2) {
    tau2 <- c(quadterm + b) / (N + a)
  } else tau2 <- NULL
  
  return(list(llik = logl, tau2 = tau2))
}

loglw_vec_N <- function(out_vec_N, out_vec, approx, A, nugs, theta, outer = TRUE, 
                        v, calc_tau2 = TRUE, sep = FALSE, mean = 0, scale = 1, a , b) {
  
  # a <- b <- 1.5
  
  N <- length(out_vec_N)
  n <- length(out_vec)
  
  out_vec_N <- out_vec_N - mean
  out_vec <- out_vec - mean
  
  # order according to X's
  A_ord <- A[approx$ord]
  if(length(nugs) ==1) nugs <- rep(nugs, length(out_vec)) # homGP -> repeat same nugs
  g_ord <- nugs[approx$ord] # order lambdas
  G <- g_ord/A_ord # calculate adjusted Lam (n) or g(n) * A_inv 
  
  U_mat <- create_U(approx, G, theta, v, sep = sep) / sqrt(scale) # does either isotropic or separable
  Ln.inv <- 1/g_ord # Lam(n) ^ (-1)
  
  LN.inv <- rep(Ln.inv, A_ord)
  
  Uty <- Matrix::crossprod(U_mat, out_vec) 
  ytUUty <- sum(Uty^2)
  
  term1 <- sum((out_vec_N^2) * LN.inv)
  term2 <- sum((out_vec^2) * (A_ord * Ln.inv))
  
  quadterm <- term1 - term2 + ytUUty
  
  logdet <- sum(log(Matrix::diag(U_mat)))
  
  # linear term
  S <- (-0.5) * sum((A_ord -1) * log(g_ord) + log(A_ord))
  
  # adjust prof likelihood
  if (outer) {
    # logl <- logdet - ((N * 0.5) * log(quadterm) + S
    logl <- logdet - ((N + a) * 0.5) * log((quadterm + b)) + S
  } else {
    logl <- logdet - 0.5 * quadterm + S
  }
  
  if(calc_tau2) {
    tau2 <- c(quadterm + b) / (N + a)
  } else tau2 <- NULL
  
  return(list(llik = logl, tau2 = tau2))
}

loglw_iso_N <- function(y, yn, dx_n, nugs, A, v, theta, outer = TRUE, calc_tau2 = TRUE, 
                        mean = 0, scale = 1, a, b){
  
  # Lam_N * A_inv
  lam_adj <- nugs/A
  
  if(length(nugs) == 1) lamN <- rep(nugs, sum(A)) # homGP
  else lamN <- rep(nugs, A) # hetGP
  
  if(v == 999){
    Kn <- Exp2vec(dx_n, tau2 = scale, theta = theta, g = lam_adj)
  }else{
    if(v > 1000) v = v - 1000
    Kn <- MaternVec(dx_n, tau2 = scale, theta = theta, g = lam_adj, v= v)
  }
  
  N <- length(y)
  n <- length(yn)
  
  # Kinv <- solve(Kn)
  Kinv <- invdet(Kn)
  # if(mean != 0) mean <- drop(rowSums(Kinv) %*% yn/(sum(Kinv))) # 1t Kn(-1) Delta (1t Kn(-1) 1) ^ (-1)
  y <- y - mean
  yn <- yn - mean
  
  q1 <- sum(y^2 /lamN)  # Y(N)t * Lam(N) ^ (-1) * Y(N)
  q2 <- sum(yn^2 * A/nugs) # Y(n)t * Lam(n) ^ (-1) * Y(n)
  
  quadterm <-  q1 - q2 + t(yn) %*% Kinv$Mi %*% yn#solve(Kn,yn) 
  
  logdet <- (-0.5) * Kinv$ldet#determinant(Kn, logarithm=TRUE)$modulus ## log (det(Kn)) 
  
  S <- (-0.5) * sum((A-1) * log(nugs) + log(A))
  
  if (outer) { # prof llik
    logl <- logdet - (N * 0.5) * log(quadterm) + S
  } else {
    logl <- logdet - 0.5 * quadterm + S
  }
  
  if(calc_tau2) { # if prof llik and wish to infer tau2
    tau2 <- c(quadterm + b) / (N + a)
  } else tau2 <- NULL
  
  return(list(llik = logl, tau2 = tau2))
}
