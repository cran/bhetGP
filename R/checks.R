# Priors for all the hyper-parameters
check_settings <- function(settings, gp_type, x, y) {
  
  dab <- laGP::darg(NULL, x)$ab
  # gab <- laGP::garg(NULL, y)$ab
  
  # priors for proposal distribution for Length-scales
  if (is.null(settings$l)) settings$l <- 1
  if (is.null(settings$u)) settings$u <- 2
  
  # Priors for ls - Latent layer (Lambdas)
  if(gp_type == "hetgp"){
    if (is.null(settings$alpha$theta_lam)) settings$alpha$theta_lam <- dab[1]
    if (is.null(settings$beta$theta_lam)) settings$beta$theta_lam <- dab[2]
  }
  #else if(gp_type == "homgp"){ # initialize nugget hypers 
  if (is.null(settings$alpha$g)) settings$alpha$g <- 1.5
  if (is.null(settings$beta$g)) settings$beta$g <- 4
  #}
  
  # Priors for LS - y
  if (is.null(settings$alpha$theta_y)) settings$alpha$theta_y <- dab[1]
  if (is.null(settings$beta$theta_y)) settings$beta$theta_y <- dab[2]
  
  if (is.null(settings$a$tau2_y)) settings$a$tau2_y <- 10/2
  if (is.null(settings$b$tau2_y)) settings$b$tau2_y <- 4/2
  
  if (is.null(settings$a$tau2_lam)) settings$a$tau2_lam <- 10/2
  if (is.null(settings$b$tau2_lam)) settings$b$tau2_lam <- 4/2
  
  return(settings)
}

# for hetGP
check_initialization <- function(reps, initial, n, D, sep, v = NULL, vec = FALSE, verb = TRUE, 
                                 stratergy) {
  
  # first run inits() and then assign accordingly if any thing not specified?
  if(is.null(initial$prof_ll_lam)){
    initial$prof_ll_lam <- TRUE
    initial$inner <- TRUE # for Lam layer
    initial$inner_tau2 <- TRUE
    initial$scale_lam <- 1
  } 
  
  if(initial$prof_ll_lam){
    if(initial$scale_lam != 1){
      warning("cannot specify scale and use prof ll; will not infer tau2 lam")
      initial$inner <- initial$inner_tau2 <- FALSE # for Lam layer
    }
  }else
    initial$inner <- initial$inner_tau2  <- FALSE # for Lam layer
  
  if(verb) print("obtaining initialization")
  
  if(stratergy == "default"){
    # if(ncol(reps$X0) > 1 || nrow(reps$X0) <= 1000)
      init <- inits_bhetgp(reps, v = v, vec = vec)
    # else{
    #   warning("vec = T but dim(x) = 1. SVec needs dim(x) > 2. Setting stratergy = flat")
    #   stratergy = "flat"
    # }  
  }
  if (stratergy == "flat")
  
  if(verb) print("checking specifications")
  if(is.null(initial$mean_y)) initial$mean_y <- mean(reps$Z)
  if(is.null(initial$mean_lam)) initial$mean_lam <- init$mean0
  
  # if(is.null(initial$scale_y)) initial$scale_y <- init$scale_y
  if(!initial$prof_ll_lam) initial$scale_lam <- init$scale_lam
  
  initial$outer <- TRUE # for Y layer
  initial$tau2 <- TRUE
  
  if(is.null(initial$noise)) initial$noise = FALSE

  # nugget for lam layer numeric stability]
  if (is.null(initial$g)) initial$g <- 1e-5 
  
  if(is.null(initial$theta_lam)) {
    initial$theta_lam <- ifel(sep, init$tg, mean(init$tg))
  }
  if(is.null(initial$theta_y)) {
    initial$theta_y <- ifel(sep, init$ty, mean(init$ty))
  }
  
  if(sep){
    if (length(initial$theta_lam) == 1) 
      initial$theta_lam <- rep(initial$theta_lam, D)
    else if(length(initial$theta_lam) != D)
      stop("bad theta_lam; does not match dimensions")
    
    if(length(initial$theta_y) == 1)
      initial$theta_y <- rep(initial$theta_y, D)
    else if(length(initial$theta_y) != D)
      stop("bad theta_y; does not match dimensions")
  }else{
    if(length(initial$theta_lam) != 1) stop("enter scalar")
    if(length(initial$theta_y) != 1) stop("enter scalar")
  }
  
  if(is.null(initial$theta_check)) initial$theta_check <- FALSE
  else if(!is.logical(initial$theta_check)) stop("must be TRUE/FALSE")
  
  llam <- initial$llam
  if (!is.null(llam) & length(llam) != n) stop("wrong number of nuggets")
  else if (is.null(llam))  initial$llam <- init$llam_w
  # initialize lam values
  
  if (!is.matrix(initial$llam)) 
    initial$llam <- as.matrix(initial$llam)
  # ensure n lambdas
  
  if(nrow(initial$llam) != n){
    initial$llam <- matrix(rep(initial$llam, times = init$reps$mult), ncol = 1)
  }
  if (ncol(initial$llam) != 1) 
    stop("latent layer must have only one dimension")
  
  return(initial)
}

check_initialization_hom <- function(reps, initial, n, D, sep, v, vec = FALSE, verb = TRUE,
                                     stratergy) {
  
  if(verb) print("obtaining initialization")
  if(stratergy == "default"){
    # if(ncol(reps$X0) > 1 || nrow(reps$X0) <= 1000)
      init <- inits_bhomgp(reps, v = v, vec = vec)
    # else{
    #   warning("vec = T but dim(x) = 1. SVec needs dim(x) > 2. Setting stratergy = flat")
    #   stratergy = "flat"
    # } 
  }
  if(stratergy == "flat")
    init <- list(ty = rep(0.1, ncol(reps$X0)), mu_y= 0, scale = 1, g = var(reps$Z)*0.1)
  
  if(verb) print("checking specifications")
  if(is.null(initial$mean_y)) initial$mean_y <- mean(reps$Z)
  if(is.null(initial$scale_y)) initial$scale_y <- 1
  
  # can just directly set true in gibbs
  if(initial$prof_ll == TRUE || is.null(initial$prof_ll) == TRUE){
    initial$outer <- TRUE
    initial$tau2 <- TRUE
    initial$scale_y <- 1
  }else{
    initial$outer <- FALSE
    initial$tau2 <- FALSE
    if(is.null(initial$scale_y)){
      warning("using MLE estimate for scale")
      initial$scale_y <- init$scale 
    } 
  }
  
  if(is.null(initial$noise)) initial$noise <- TRUE
  
  # initialize nugget
  if(initial$noise){ # sample noise
    if(is.null(initial$g)) initial$g <- init$g
  }else{# no noise in model
      warning("very small noise for numeric stability")
      if(is.null(initial$g)) initial$g <- 1e-8 
    }
  
  if (is.null(initial$theta_y)) initial$theta_y <- init$ty
  
  if(sep){
    if(length(initial$theta_y) == 1)
      initial$theta_y <- rep(initial$theta_y, D)
    else if(length(initial$theta_y) != D)
      stop("bad theta_y; does not match dimensions")
  }
  
  return(initial)
}

check_inputs <- function(x, y, tau2) {
  
  if (!is.vector(y)) 
    stop("y must be a vector")
  if (nrow(x) != length(y)) 
    stop("dimensions of x and y do not match")
  if (min(x) < -5 | min(x) > 5 | max(x) < -4 | max(x) > 6) 
    warning("this function is designed for x over the range [0, 1]")

  if(tau2 > 100)
    warning("standardize response with variance 1 for better results")
  # if (max(y) > 3 * sqrt(tau2) | min(y) < (-3) * sqrt(tau2))
  #   warning("standardize response for mean zero and variance 1 for better results")
  
  return(NULL)
}

check_ordering <- function(ordering, n) {
  
  if (length(ordering) != n) 
    stop("length(ordering) must match dimension of x for fitting or x_new for predicting")
  if (min(ordering) != 1)
    stop("ordering must begin at index 1")
  if (max(ordering) != n)
    stop("ordering must end at index nrow(x) for fitting or nrow(x_new) for predicting")
  if (sum(duplicated(ordering)) > 0)
    stop("ordering must not have any duplicates")
  
  return(NULL)
}

check_scale <- function(x, v, vec = FALSE, priors, init){
  
  N <- length(init$llam)
  count <- 0
  init$tau2_lam <- 1000 # arbitrary
  while(init$tau2_lam > 10){
    init$theta_lam <- init$theta_lam/10
    obj <- precalc(x, out = NULL, outs2 = NULL, A = NULL, v, vec, priors, init, 
                   latent = TRUE)
    init$tau2_lam <- obj$tau2
    init$llik_lam <- obj$llik
    count <- count + 1
    if(count == 10) init$g <- init$g * 10
  }
  
  if(init$theta_check == TRUE & any(init$theta_lam < init$theta_y)) init$theta_y = init$theta_lam
  
  return(init)
}

check_inits_vdims <- function(reps, initial, n, D, sep, vdims = NULL, 
                              nlam, reps_vdims, v = NULL, stratergy, 
                              vec, verb = TRUE) {
  
  if(verb) print("obtaining initialization")
  if(stratergy == "default"){
    init <- inits_vdims(reps, reps_vdims, v = v, vec = vec)
  }
  if (stratergy == "flat"){
    out <- list(reps = reps, ty = rep(1, D), tg = rep(2, length(vdims)), 
                mean0 = 0, scale0 = 1, llam_w = rep(log(var(reps$Z) * 0.1), nlam), 
                mu_y = 0, scale = 1, g = 1e-5)
  }
    
  
  if(verb) print("checking specifications")
  if(is.null(initial$mean_y)) initial$mean_y <- mean(reps$Z)
  if(is.null(initial$mean_lam)) initial$mean_lam <- init$mean0
  
  # if(is.null(initial$scale_y)) initial$scale_y <- init$scale_y
  if(is.null(initial$scale_lam)) initial$scale_lam <- init$scale_lam
  
  if(is.null(initial$prof_ll_lam)){
    initial$prof_ll_lam <- TRUE
    initial$scale_lam <- 1
  } 
  
  if(initial$prof_ll_lam){
    initial$inner <- TRUE # for Lam layer
    initial$inner_tau2 <- TRUE
    # if(initial$scale_lam != 1) stop("cannot specify scale and use prof ll")
  }else{
    initial$inner <- FALSE # for Lam layer
    initial$inner_tau2 <- FALSE
  }
  
  initial$outer <- TRUE # for Y layer
  initial$tau2 <- TRUE
  
  if(is.null(initial$noise)) initial$noise = FALSE
  
  # nugget for lam layer numeric stability]
  if (is.null(initial$g)) initial$g <- 1e-5 
  
  if(is.null(initial$theta_lam)) initial$theta_lam <- init$tg[vdims]
  if(is.null(initial$theta_y)) initial$theta_y <- init$ty
  
  if(sep){
    if(length(initial$theta_lam) != length(vdims)) 
      stop("bad theta_lam; match length of vdims")
    if(length(initial$theta_y) == 1)
      initial$theta_y <- rep(initial$theta_y, D)
    else if(length(initial$theta_y) != D)
      stop("bad theta_y; does not match dimensions")
  }else{
    if(length(initial$theta_lam) != 1) stop("enter scalar")
    if(length(initial$theta_y) != 1) stop("enter scalar")
  }
  
  if(is.null(initial$theta_check)) initial$theta_check <- FALSE
  else if(!is.logical(initial$theta_check)) stop("must be TRUE/FALSE")
  
  llam <- initial$llam
  # initialize lam values
  if (is.null(llam)){
      llam <- llam_nv <- rep(0.1, nlam)
      initial$llam_nv <- llam_nv
      initial$lam_nv <- exp(llam_nv)
      if(!is.null(reps_vdims$mult)) llam <- rep(llam, reps_vdims$mult)[reps_vdims$mult$Z] # make it len of input space and order
    }
  
  # initial$lam_0 <- lam
  initial$llam <- llam
  
  if (!is.matrix(initial$llam)) 
    initial$llam <- as.matrix(initial$llam)
  # ensure n lambdas
  if(nrow(initial$llam) != n){
    warning("wrong number of nuggets") 
    initial$llam <- matrix(rep(initial$llam, times = init$reps$mult), ncol = 1)
  }
  if (ncol(initial$llam) != 1) 
    stop("latent layer must have only one dimension")
  
  return(initial)
}

check_scale_vdims <- function(x, v, init, vec = FALSE, sep){
  
  N <- length(init$llam)
  tau2 <- 1000
  init$theta_lam <- init$theta_lam * 10
  count <- 0
  while(tau2 > 10){
    # print(tau2)
    init$theta_lam <- init$theta_lam/10
    if(!vec){
      if(sep){
        if(v == 999) KN <- Exp2Sep(x, x, tau2 = init$scale, theta = init$theta_lam, g = init$g_0)
        else if(v > 1000) KN <- MaternProdSep(x, x, tau2 = init$scale, theta = init$theta_lam, g = init$g_0, v= (v - 1000))
        else KN <- MaternSep(x, x, tau2 = init$scale, theta = init$theta_lam, g = init$g_0, v = v)
      } 
      else{
        if(v > 1000) v = v - 1000
        if(v == 999) KN <- Exp2(sq_dist(x), tau2 = init$scale, theta = init$theta_lam, g = init$g_0)
        else KN <- Matern(sq_dist(x), tau2 = init$scale, theta = init$theta_lam, g = init$g_0, v= v)
      } 
      Kinv <- solve(KN)
      quadterm <- (t(init$llam_nv) %*% Kinv %*% init$llam_nv)
      tau2 <- c(quadterm)/N
    }else{
      llam_ord <- init$llam_nv[x$ord]
      U_mat <- create_U(x, init$g_0, init$theta_lam, v = v, sep = sep) / sqrt(init$scale)
      Uty <- Matrix::crossprod(U_mat, llam_ord) # Ut * Y
      ytUUty <- sum(Uty^2) 
      tau2 <- c(ytUUty)/N
    }
    count <- count + 1
    if(count == 10) init$g_0 <- init$g_0 * 10
  }
  
  init$scale_lam <- tau2
  if(init$theta_check == TRUE & any(init$theta_lam < init$theta_y)) init$theta_y = init$theta_lam
  
  return(init)
}

# add something for isotropic GP
precalc <- function(x, out = NULL, outs2 = NULL, A = NULL, v, vec, priors, init, 
                    latent = TRUE){

  out_vec <- ifel(latent, init$llam, out)
  theta <- ifel(latent, init$theta_lam, init$theta_y)
  outer <- ifel(latent, init$inner, init$outer)
  tau2 <- ifel(latent, init$inner_tau2, TRUE)
  mean <- ifel(latent, init$mean_lam, init$mean_y)
  a <- ifel(latent, priors$a$tau2_lam, priors$a$tau2_y)
  b <- ifel(latent, priors$b$tau2_lam, priors$b$tau2_y)
  D <- ifel(vec, ncol(x$x_ord), ncol(x))
  
  if(length(theta) == 1) theta <- rep(theta, D)
  
  if(!is.null(init$llam))
    nugs <- ifel(latent, init$g, exp(init$llam))
  else
    nugs <- init$g
  
  if(!vec){
    if(is.null(outs2)) calc <- logl(out_vec, x, nugs, theta, v, outer, tau2, mean, scale = 1, a , b)
    else calc <- loglw(outs2, out_vec, x, nugs, A, v, theta, outer, tau2, mean, scale = 1, a , b)
  }else{
    if(is.null(outs2)) calc <- logl_vec(out_vec, x, nugs, theta, outer, v, tau2,
                                        sep = TRUE, mean, scale = 1, latent, a , b)
    else calc <- loglw_vec(outs2, out_vec, x, A, nugs, theta, outer, v, tau2, 
                           sep = TRUE, mean, scale = 1, a , b)
        
  }
  
  return(calc)
}


