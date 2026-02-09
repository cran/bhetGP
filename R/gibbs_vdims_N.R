# Variance changes only in some dimensions of the input space in this version
# Full N i.e. no reps case

# ---------- Functions --------------------------------------------------------
# gibbs_sep_vdims_N: gibbs sampler for seperable theta
# gibbs_iso_vdims_N: gibbs for isotropic GP 
# ess_sample_N : ess sampling for full N version with vdims
# gibbs_sep_vdims_N_vec: gibbs sampler for seperable theta w/ vecchia
# gibbs_iso_vdims_N_vec: gibbs for isotropic GP w/ vecchia
# ess_sample_N_vec: ess sampling for full N version with vdims w/ vecchia
# -----------------------------------------------------------------------------

gibbs_sep_vdims_N <- function(Yn, xm, mcmc, initial, priors, v, vdims = NULL, verb = TRUE){
  
  # mean process
  nm <- nrow(xm)
  D <- ncol(xm)
  
  # variance process input space
  xv <- xm[ ,vdims , drop = FALSE]
  nv <- nrow(xv)
  d <- length(vdims)

  N <- length(Yn)
  
  # initialize
  llam_N <- matrix(nrow = mcmc, ncol = nm)
  theta_y <- matrix(nrow = mcmc, ncol = D)
  theta_lam <- matrix(nrow = mcmc, ncol = d)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_N[1, ] <- initial$llam
  
  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- g <- NULL
  llik_y_store <- llik_lam_store <- tau2 <- tau2_lam <- tau2_lam_store <- tau2_store <- NULL
  
  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y != 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g

  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # theta lam sampling
    for(i in 1:d){
      obj_lam <- mcmc_theta_l_sep(xv, llam_N[t-1, ], index = i, theta_cur_lam, llik_prev = llik_lam, 
                                  tau2_prev = tau2_lam, v = v, ls_check = initial$theta_check, 
                                  theta_y = theta_cur_y[vdims[i]], g0 = g[t - 1], alpha = priors$alpha$theta_lam, 
                                  beta = priors$beta$theta_lam, l = priors$l , u = priors$u,
                                  inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                                  mean0 = initial$mean_lam, scale0 = initial$scale_lam,
                                  a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam)
      
      theta_lam[t, ] <- obj_lam$theta
      theta_cur_lam <- obj_lam$theta
      llik_lam <- obj_lam$llik
      tau2_lam <- obj_lam$tau2
    }
    
    if(initial$noise){
      obj_g <- mcmc_g_sep_inner(llam_N[t - 1, ], xv, g[t - 1], llik_prev = llik_lam, tau2_prev = tau2_lam, 
                                theta = theta_lam[t, ], v = v, outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                                mean = initial$mean_lam, scale = initial$scale_lam, alpha = priors$alpha$g, 
                                beta = priors$beta$g + t, l = priors$l , u = priors$u,
                                a = priors$a$tau2_lam, b = priors$b$tau2_lam)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
      
    }else g[t] <- initial$g

    # theta y sampling
    for(i in 1:D){
      
      obj_y <- mcmc_theta_y_sep(YNs2 = NULL, yn= Yn, xn= xm, A = NULL, index = i, llam_N[t-1, ], 
                                theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, v = v,
                                ls_check = initial$theta_check, 
                                theta_lam = ifel(any(i == vdims), theta_cur_lam[which(i == vdims)], NULL),
                                outer = initial$outer, calc_tau2 = initial$tau2,
                                mean = initial$mean_y, scale = initial$scale_y, alpha = priors$alpha$theta_y, 
                                a = priors$a$tau2_y, b = priors$b$tau2_y,
                                beta = priors$beta$theta_y, l = priors$l, u = priors$u)
      
      theta_y[t, ] <- obj_y$theta
      theta_cur_y <- obj_y$theta
      llik_y <- obj_y$llik
      tau2 <- obj_y$tau2
    }
    
    # ess sampling for N lambdas
    llam_draw <- ess_sample_N(yn= Yn, xn= xm, xv, llam_N[t - 1, ], llik_prev = llik_y, 
                              theta_y[t, ], theta_lam = theta_lam[t, ], 
                              sep = TRUE, v = v, outer = initial$outer, calc_tau2 = initial$tau2, 
                              dx_n = NULL, dx_v = NULL, mean = initial$mean_y, scale = initial$scale_y, 
                              a = priors$a$tau2_y, b = priors$b$tau2_y,
                              mean0 = initial$mean_lam, 
                              scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam),
                              g0 = g[t])
    
    llam_N[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2

    obj_lam <- logl(llam_N[t, ], xv, nugs = g[t], theta_lam[t, ], v= v, outer = initial$inner, 
                    calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, scale = initial$scale_lam, 
                    a = priors$a$tau2_lam, b = priors$b$tau2_lam)
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2
    
    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)     
    tau2_store[t] <- tau2
  }
  
  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_N, tau2 = tau2_store, 
              tau2_lam = tau2_lam_store, g= g, llik_y = llik_y_store, llik_lam = llik_lam_store))
}

gibbs_iso_vdims_N <- function(Yn, xm, mcmc, initial, priors, v, vdims = NULL, verb = TRUE){
  
  # mean and variance process
  nm <- nrow(xm)
  D <- ncol(xm)
  
  xv <- xm[ ,vdims , drop = FALSE]
  nv <- nrow(xv)

  N <- length(Yn)
  
  # initialize
  llam_N <- matrix(nrow = mcmc, ncol = nm)
  theta_y <- matrix(nrow = mcmc, ncol = 1)
  theta_lam <- matrix(nrow = mcmc, ncol = 1)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_N[1, ] <- initial$llam
  
  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- g <- NULL
  llik_y_store <- llik_lam_store <- tau2 <- tau2_lam <- tau2_lam_store <- tau2_store <- NULL
  
  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y != 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g
  
  # pre calc distances
  dxv <- sq_dist(xv)
  dxm <- sq_dist(xm)
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # theta lam sampling
    obj_lam <- mcmc_theta_l(dxv, llam_N[t - 1, ], theta_cur_lam, llik_prev = llik_lam, tau2_prev = tau2_lam, 
                            v = v, ls_check = initial$theta_check, theta_y = theta_cur_y, g0 = g[t - 1], 
                            alpha = priors$alpha$theta_lam, beta = priors$beta$theta_lam, l = priors$l,
                            u = priors$u, inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                            mean0 = initial$mean_lam, scale0 = initial$scale_lam,
                            a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam)
    
    theta_lam[t] <- obj_lam$theta
    theta_cur_lam <- obj_lam$theta
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2
    
    if(initial$noise){
      obj_g <- mcmc_g_inner(llam_N[t-1, ], dxv, g[t - 1], llik_prev = llik_lam, tau2_prev = tau2_lam, 
                            theta = theta_lam[t], v = v, outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                            mean = initial$mean_lam, scale = initial$scale_lam, alpha = priors$alpha$g, 
                            beta = priors$beta$g + t, l = priors$l , u = priors$u,
                            a = priors$a$tau2_y, b = priors$b$tau2_y)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
    }else g[t] <- initial$g

    # theta_y sampling
    obj_y <- mcmc_theta_y(YNs2 = NULL, Yn, dxm, A = NULL, llam_N[t - 1, ], theta_cur_y, llik_prev = llik_y, 
                          tau2_prev = tau2, v = v, ls_check = initial$theta_check, theta_lam = theta_cur_lam,
                          outer = initial$outer, calc_tau2 = initial$tau2,
                          mean = initial$mean_y, scale = initial$scale_y,  alpha = priors$alpha$theta_y, 
                          beta = priors$beta$theta_y, l = priors$l, u = priors$u,
                          a = priors$a$tau2_y, b = priors$b$tau2_y)
    theta_y[t] <- obj_y$theta
    theta_cur_y <- obj_y$theta
    llik_y <- obj_y$llik
    tau2 <- obj_y$llik
    
    # Ess sampling for lambdas
    llam_draw <- ess_sample_N(yn= Yn, xn= xm, xv, llam_N[t - 1, ], llik_prev = llik_y, theta_y[t, ],
                              theta_lam = theta_lam[t, ], sep = FALSE, v = v, outer = initial$outer, 
                              calc_tau2 = initial$tau2, a = priors$a$tau2_y, b = priors$b$tau2_y,
                              dx_n = dxm, dx_v = dxv, mean = initial$mean_y, scale = initial$scale_y, 
                              mean0 = initial$mean_lam, 
                              scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam),
                              g0 = g[t])
    
    llam_N[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2

    obj_lam <- logl_iso(llam_N[t, ], dxv, nugs = g[t - 1], theta_lam[t, ], v = v, outer = initial$inner, 
                        calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, scale = initial$scale_lam,
                        a = priors$a$tau2_lam, b = priors$b$tau2_lam)
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2
    
    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)     
    tau2_store[t] <- tau2
  }
  
  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_N, tau2 = tau2_store,
              tau2_lam = tau2_lam_store, g = g, llik_y = llik_y_store, llik_lam = llik_lam_store))
}

ess_sample_N <- function(yn, xn, xv, llam_prev, llik_prev = NULL, theta_y, theta_lam, sep = TRUE, 
                         v, outer, calc_tau2, mean, scale, dx_n = NULL, dx_v = NULL, 
                         mean0, scale0, g0, a, b) { 
  
  # prev llik
  if(is.null(llik_prev)){
    if(sep) llik_prev <- logl(yn, xn, exp(llam_prev), theta_y, v = v, outer = outer, calc_tau2 = calc_tau2, 
                              mean = mean, scale = scale, a= a, b = b)$llik
    else llik_prev <- logl_iso(yn, dx_n, exp(llam_prev), theta_y, v = v, outer = outer, calc_tau2 = calc_tau2, 
                               mean = mean, scale = scale, a= a, b = b)$llik 
  }
  
  # llam_prior 
  llam_prior <- mvn_prior(x = xv, v, theta_lam, g0, dx = dx_v, mean0, scale0, sep = sep)  
  # angle
  gam <- runif(1,0,2*pi)
  
  # bounds
  gam_min <- gam - 2*pi
  gam_max <- gam
  
  while (1){
    
    llam_new <- llam_prev * cos (gam) + llam_prior * sin (gam) # proposal
    
    # calculate llik
    if(sep)
      obj_prop <- logl(yn, xn, nugs = exp(llam_new), theta_y, v = v, outer = outer, calc_tau2 = calc_tau2, 
                       mean = mean, scale = scale, a= a, b = b)
    else 
      obj_prop <- logl_iso(yn, dx_n, nugs = exp(llam_new), theta_y, v = v, outer = outer, calc_tau2 = calc_tau2, 
                           mean = mean, scale = scale, a= a, b = b)
    llik_prop <- obj_prop$llik
    tau2 <- obj_prop$tau2
    
    # alpha <- min(0, llik_prop - llik_prev)
    r <- runif(1,0,1)
    # check acceptance  
    if(llik_prop < llik_prev + log(r)){
      # Update bounds
      if(gam < 0) gam_min <- gam 
      else gam_max <- gam 
      gam <- runif(1,gam_min,gam_max)
    }
    else
      break
  }
  
  return(list(llam = llam_new, llik = llik_prop, tau2 = tau2)) # log lam new
}

gibbs_sep_vdims_N_vec <- function(Yn, xm_approx, mcmc, initial, priors, v, vdims, verb = TRUE){
  
  # mean and variance process
  Xn <- as.matrix(xm_approx$x_ord[xm_approx$rev_ord_obs, , drop = FALSE])
  nm <- nrow(Xn)
  D <- ncol(Xn)
  
  xv <- Xn[ , vdims , drop = FALSE]
  nv <- nrow(xv)
  d <- length(vdims)

  N <- nrow(Yn)
  
  # create approx for original xs
  xv_approx <- create_approx(xv, xm_approx$m)
  
  llam_N <- matrix(nrow = mcmc, ncol = nm)
  theta_y <- matrix(nrow = mcmc, ncol = D)
  theta_lam <- matrix(nrow = mcmc, ncol = d)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_N[1, ] <- initial$llam

  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- g <- NULL
  llik_y_store <- llik_lam_store <- tau2 <- tau2_lam <- tau2_lam_store <- tau2_store <- NULL
  
  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y != 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    for(i in 1:d){
      obj_lam <- mcmc_theta_l_sep_vec(xv_approx, llam_N[t - 1, ], index = i, theta_cur_lam, 
                                      llik_prev = llik_lam, tau2_prev = tau2_lam, v = v,
                                      ls_check = initial$theta_check, theta_y = theta_cur_y[vdims[i]], 
                                      g0 = g[t - 1],  alpha = priors$alpha$theta_lam, 
                                      beta = priors$beta$theta_lam, l = priors$l , u = priors$u,
                                      inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                                      a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam,
                                      mean0 = initial$mean_lam, scale0 = initial$scale_lam) 

      theta_lam[t, ] <- obj_lam$theta
      theta_cur_lam <- obj_lam$theta
      llik_lam <- obj_lam$llik
      tau2_lam <- obj_lam$tau2
    }

    if(initial$noise){
      obj_g <- mcmc_g_vec_inner(llam_N[t - 1, ], xv_approx, g[t - 1], llik_prev = llik_lam,
                                tau2_prev = tau2_lam, theta = theta_lam[t, ], v = v, outer = initial$inner, 
                                calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, 
                                scale = initial$scale_lam, alpha = priors$alpha$g, beta = priors$beta$g + (t),
                                l = priors$l , u = priors$u, sep = TRUE, 
                                a = priors$a$tau2_lam, b = priors$b$tau2_lam)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
    }else g[t] <- initial$g
    
    
    for(i in 1:D){
      obj_y <- mcmc_theta_y_sep_vec(YNs2 = NULL, yn= Yn, xm_approx, A = NULL, index = i, llam_N[t - 1, ], 
                                    theta_prev = theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, 
                                    v = v, ls_check = initial$theta_check, 
                                    theta_lam = ifel(any(i == vdims), theta_cur_lam[which(i == vdims)], NULL),
                                    outer = TRUE, calc_tau2 = TRUE, mean = initial$mean_y, scale = initial$scale_y, 
                                    alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                                    a = priors$a$tau2_y, b = priors$b$tau2_y,
                                    l = priors$l , u = priors$u)
      theta_y[t, ] <- obj_y$theta
      theta_cur_y <- obj_y$theta
      llik_y <- obj_y$llik
      tau2 <- obj_y$tau2
    }
    
    llam_draw <- ess_sample_N_vec(yn = Yn, xm_approx, xv_approx, llam_N[t - 1, ],
                                  llik_prev = llik_y, theta_y = theta_y[t, ], theta_lam = theta_lam[t, ], 
                                  sep = TRUE, v = v, outer = TRUE, calc_tau2 = TRUE,
                                  mean = initial$mean_y, scale = initial$scale_y, 
                                  mean0 = initial$mean_lam, 
                                  scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam),
                                  g0 = g[t],
                                  a = priors$a$tau2_y, b = priors$b$tau2_y)
    
    llam_N[t, ] <- llam_draw$llam # lambdas that map to Xn's
    
    llik_y <- llam_draw$llik # variance process llik
    tau2 <- llam_draw$tau2

    obj_lam <-  logl_vec(llam_N[t, ], xv_approx, nugs = g[t], theta_lam[t, ], outer = initial$inner, v = v, 
                         calc_tau2 = initial$inner_tau2, sep = TRUE, mean = initial$mean_lam, 
                         scale = initial$scale_lam, a = priors$a$tau2_lam, b = priors$b$tau2_lam, 
                         latent = TRUE)
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2
    
    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)     
    tau2_store[t] <- tau2
  }
  
  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_N, tau2 = tau2_store, 
              tau2_lam = tau2_lam_store, g = g, xv = xv, xv_approx = xv_approx,
              llik_y = llik_y_store, llik_lam = llik_lam_store))
}

gibbs_iso_vdims_N_vec <- function(Yn, xm_approx, mcmc, initial, priors, v, vdims = NULL, 
                                  verb = TRUE){
  
  Xn <- as.matrix(xm_approx$x_ord[xm_approx$rev_ord_obs, , drop = FALSE])
  nm <- nrow(Xn)
  D <- ncol(Xn)
  
  xv <- Xn[ , vdims , drop = FALSE]
  nv <- nrow(xv)

  N <- length(Yn)
  
  # variance layer must be vecchia layer 
  xv_approx <- create_approx(xv, xm_approx$m)
  
  llam_N <- matrix(nrow = mcmc, ncol = nm)
  theta_y <- matrix(nrow = mcmc, ncol = 1)
  theta_lam <- matrix(nrow = mcmc, ncol = 1)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  
  llam_N[1, ] <- initial$llam

  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- NULL
  llik_y_store <- llik_lam_store <- tau2 <- tau2_lam <- tau2_lam_store <- tau2_store <- NULL

  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y != 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g

  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # theta lam w/ vecchia
    obj_lam <-  mcmc_theta_l_vec(xv_approx, llam_N[t - 1, ], theta_cur_lam, llik_prev = llik_lam,
                                 tau2_prev = tau2_lam, v = v, ls_check = initial$theta_check, 
                                 theta_y = theta_cur_y, g0 = g[t - 1], alpha = priors$alpha$theta_lam, 
                                 beta = priors$beta$theta_lam, l = priors$l , u = priors$u,
                                 inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                                 mean0 = initial$mean_lam, scale0 = initial$scale_lam,
                                 a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam)
    
    theta_lam[t] <- obj_lam$theta
    theta_cur_lam <- obj_lam$theta
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2
    
    if(initial$noise){
        obj_g <- mcmc_g_vec_inner(llam_N[t - 1, ], xv_approx, g[t - 1], llik_prev = llik_lam, tau2_prev = tau2_lam,
                                  theta = theta_lam[t, ], v = v, outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                                  mean = initial$mean_lam, scale = initial$scale_lam, alpha = priors$alpha$g, 
                                  beta = priors$beta$g + (t), l = priors$l , u = priors$u, sep = FALSE,
                                  a = priors$a$tau2_lam, b = priors$b$tau2_lam)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
    }else g[t] <- initial$g
    
    obj_y <- mcmc_theta_y_vec(YNs2 = NULL, yn = Yn, x_approx = xm_approx, A = NULL, llam_N[t - 1, ],
                              theta_prev = theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, v = v,
                              ls_check = initial$theta_check, theta_lam = theta_cur_lam, outer = initial$outer, 
                              calc_tau2 = initial$tau2, mean = initial$mean_y, scale = initial$scale_y, 
                              alpha = priors$alpha$theta_y, beta = priors$beta$theta_y,
                              a = priors$a$tau2_y, b = priors$b$tau2_y, l = priors$l , u = priors$u)
    
    theta_y[t] <- obj_y$theta
    theta_cur_y <- obj_y$theta
    llik_y <- obj_y$llik
    tau2 <- obj_y$tau2
    
    llam_draw <- ess_sample_N_vec(yn = Yn, x_approx = xm_approx, xv_approx, llam_N[t - 1, ], 
                                  llik_prev = llik_y, theta_y = theta_y[t], theta_lam = theta_lam[t], 
                                  sep = FALSE, v= v, outer = initial$outer, calc_tau2 = initial$tau2,
                                  mean = initial$mean_y, scale = initial$scale_y, a = priors$a$tau2_y, 
                                  b = priors$b$tau2_y, mean0 = initial$mean_lam, 
                                  scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam),
                                  g0 = g[t], dxv = NULL)
    
    llam_N[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2

    obj_lam <-  logl_vec(llam_N[t, ], xv_approx, nugs = g[t], theta_lam[t, ], outer = initial$inner, v = v, 
                         calc_tau2 = initial$inner_tau2, sep = FALSE, mean = initial$mean_lam, 
                         scale = initial$scale_lam, a = priors$a$tau2_lam, b = priors$b$tau2_lam, latent = TRUE)
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2
    
    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)     
    tau2_store[t] <- tau2
  }
  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_N, tau2 = tau2_store, 
              tau2_lam = tau2_lam_store, g = g, xv = xv, xv_approx = xv_approx,
              llik_y = llik_y_store, llik_lam = llik_lam_store))
}

ess_sample_N_vec <- function(yn, x_approx, xv, llam_prev_nv, llik_prev = NULL,
                             theta_y, theta_lam, sep = TRUE, v, outer, calc_tau2, 
                             mean, scale, mean0, scale0 = 1, g0, dxv = NULL, a, b) { 
  
  # previous llik
  if(is.null(llik_prev)){
    llik_prev <- logl_vec(yn, x_approx, nugs = exp(llam_prev_nv), theta=theta_y, outer = outer, v = v, 
                          calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a = a, b = b)
  }
  
  # use xv to draw priors
  llam_prior_nv <- rand_mvn_vec(xv, theta_lam, v=v, mean = mean0, g = rep(g0, nrow(xv$x_ord)), 
                                scale = scale0, sep = sep) 
  
  # angle
  gam <- runif(1,0,2*pi)
  
  # bounds
  gam_min <- gam - 2*pi
  gam_max <- gam
  
  while (1){
    
    # prop lam (length N)
    llam_new <- llam_prev_nv * cos (gam) + llam_prior_nv * sin (gam)
    
    # proposed llik
    obj_prop <- logl_vec(yn, x_approx, nugs = exp(llam_new), theta=theta_y, outer = outer, v = v, 
                         calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a = a, b= b)
    llik_prop <- obj_prop$llik
    tau2 <- obj_prop$tau2
    
    # alpha <- min(0, llik_prop - llik_prev)
    r <- runif(1,0,1)
    # check acceptance  
    if(llik_prop < llik_prev + log(r)){
      # Update bounds
      if(gam < 0) gam_min <- gam 
      else gam_max <- gam 
      gam <- runif(1, gam_min, gam_max)
    }
    else
      break
  }
  
  return(list(llam = llam_new, llik = llik_prop, tau2 = tau2)) # log lam new
}
