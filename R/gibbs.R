# ---------------  Gibbs: HetGP + HomGP ----------------------------------------
# gibbs_sep: Gibbs sampler for separable HetGP
# gibbs_iso: Gibbs sampler for isotropic HetGP
# gibbs_sep_vdims: Gibbs for Sep HetGP where variance changed in some dims
# gibbs_iso_vdims: Gibbs for Iso HetGP where variance changed in some dims
# gibbs_sep_hom: Gibbs for Sep HomGP 
# gibbs_iso_hom: Gibbs for Iso HomGP 
# ------------------------------------------------------------------------------

gibbs_sep <- function(YNs2 = NULL, Yn, Xn, A, mcmc, initial, priors, v, verb = TRUE){

  n <- nrow(Xn)
  D <- ncol(Xn)
  N <- sum(A)
  
  # initial and assign
  llam_samples <- matrix(nrow = mcmc, ncol = n)
  theta_y <- matrix(nrow = mcmc, ncol = D)
  theta_lam <- matrix(nrow = mcmc, ncol = D)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_samples[1, ] <- initial$llam
  
  llik_y_store <- llik_lam_store <- tau2_store <- tau2_lam_store <- rep(NA, mcmc)
  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- NULL
  tau2 <- tau2_lam <- g <- NULL

  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y!= 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g
# stop("test")
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # theta lam updates
    for(i in 1:D){
      obj_lam <- mcmc_theta_l_sep(Xn, llam_samples[t-1, ], index = i, theta_cur_lam, 
                                  llik_prev = llik_lam, tau2_prev = tau2_lam, v = v,
                                  ls_check = initial$theta_check, theta_y = theta_cur_y[i], 
                                  g0 = g[t - 1], alpha = priors$alpha$theta_lam, 
                                  beta = priors$beta$theta_lam, l = priors$l , u = priors$u,
                                  a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam,
                                  inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                                  mean0 = initial$mean_lam, scale0 = initial$scale_lam)
      
      theta_lam[t, ] <- obj_lam$theta
      theta_cur_lam <- obj_lam$theta
      llik_lam <- obj_lam$llik
      tau2_lam <- obj_lam$tau2
    }

    if(initial$noise){
      obj_g <- mcmc_g_sep_inner(llam_samples[t-1, ], Xn, g[t - 1], llik_prev = llik_lam,
                                tau2_prev = tau2_lam, theta = theta_lam[t, ], v = v, 
                                outer = initial$inner, calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, 
                                scale = initial$scale_lam, alpha = priors$alpha$g, 
                                beta = priors$beta$g + t, l = priors$l , u = priors$u,
                                a = priors$a$tau2_lam, b = priors$b$tau2_lam)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2

    }else g[t] <- initial$g
    
    # theta y updates
    for(i in 1:D){
      obj_y <- mcmc_theta_y_sep(YNs2 = YNs2, yn= Yn, xn=Xn, A = A, index = i, llam_samples[t-1, ], 
                                theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, v = v,
                                ls_check = initial$theta_check, theta_lam = theta_cur_lam[i],
                                outer = initial$outer, calc_tau2 = initial$tau2, mean = initial$mean_y, 
                                scale = initial$scale_y, alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                                l = priors$l, u = priors$u, a = priors$a$tau2_y, b = priors$b$tau2_y)
        
        theta_y[t, ] <- obj_y$theta
        theta_cur_y <- obj_y$theta
        llik_y <- obj_y$llik
        tau2 <- obj_y$tau2
    }
    
    # Log Lam updates
    llam_draw <- ess_sample(YNs2= YNs2, yn=Yn, xn=Xn, llam_samples[t - 1, ], A = A, llik_prev = llik_y,
                            theta_y = theta_y[t, ], theta_lam = theta_lam[t, ],  sep = TRUE, v = v, 
                            outer = initial$outer, calc_tau2 = initial$tau2, mean0 = initial$mean_lam, 
                            scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam),
                            g0 = g[t], a = priors$a$tau2_y, b = priors$b$tau2_y,
                            mean = initial$mean_y, scale = initial$scale_y)
    llam_samples[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2

    obj_lam <-  logl(llam_samples[t, ], Xn, nugs = g[t], theta_lam[t, ], v= v, outer = initial$inner, 
                     calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, scale = initial$scale_lam,
                     a = priors$a$tau2_lam, b = priors$b$tau2_lam)
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2

    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)     
    tau2_store[t] <- tau2
  }
  
  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_samples, g = g,
              tau2 = tau2_store, tau2_lam = tau2_lam_store, llik_y = llik_y_store, llik_lam = llik_lam_store))
}

gibbs_iso <- function(YNs2 = NULL, Yn, Xn, A, v, mcmc, initial, priors, verb = TRUE){
  
  n <- nrow(Xn)
  D <- ncol(Xn)
  N <- sum(A)
  
  # Assign space and initial values
  llam_samples <- matrix(nrow = mcmc, ncol = n)
  theta_y <- matrix(nrow = mcmc, ncol = 1)
  theta_lam <- matrix(nrow = mcmc, ncol = 1)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_samples[1, ] <- initial$llam

  llik_y_store <- llik_lam_store <- tau2_lam_store <- tau2_store <- rep(NA, mcmc)
  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- NULL
  tau2 <- tau2_lam <- g <- NULL

  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y!= 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g
  
  # Pre-calculate distance for Isotropic
  dx_n <- sq_dist(Xn)
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # theta lam updates
    obj_lam <- mcmc_theta_l(dx_n, llam_samples[t-1, ], theta_cur_lam, llik_prev = llik_lam, tau2_prev = tau2_lam, 
                            v = v, ls_check = initial$theta_check, theta_y = theta_cur_y,
                            g0 = g[t - 1], alpha = priors$alpha$theta_lam, beta = priors$beta$theta_lam, 
                            l = priors$l , u = priors$u, inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                            mean0 = initial$mean_lam, scale0 = initial$scale_lam,
                            a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam)

    theta_lam[t] <- obj_lam$theta
    theta_cur_lam <- obj_lam$theta
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2

    if(initial$noise){
      obj_g <- mcmc_g_inner(llam_samples[t-1, ], dx_n, g[t - 1], llik_prev = llik_lam, tau2_prev = tau2_lam,
                            theta = theta_lam[t], v = v, outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                            mean = initial$mean_lam, scale = initial$scale_lam, alpha = priors$alpha$g, 
                            beta = priors$beta$g + t, l = priors$l , u = priors$u,
                            a = priors$a$tau2_y, b = priors$b$tau2_y)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
    }else g[t] <- initial$g
    
    # theta y updates
    obj_y <- mcmc_theta_y(YNs2, Yn, dx_n, A, llam_samples[t-1, ], theta_cur_y, llik_prev = llik_y, 
                          tau2_prev = tau2, v = v, ls_check = initial$theta_check, theta_lam = theta_cur_lam,
                          outer = initial$outer, calc_tau2 = initial$tau2, mean = initial$mean_y, 
                          scale = initial$scale_y, alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                          l = priors$l, u = priors$u, a = priors$a$tau2_y, b = priors$b$tau2_y)
    theta_y[t] <- obj_y$theta
    theta_cur_y <- obj_y$theta
    llik_y <- obj_y$llik
    tau2 <- obj_y$tau2

    # llam updates 
    llam_draw <- ess_sample(YNs2, Yn, Xn, llam_samples[t - 1, ], A, llik_prev = llik_y, theta_y[t, ], 
                            theta_lam[t, ],  sep = FALSE, v = v, outer = initial$outer, calc_tau2 = initial$tau2,
                            scale0 =  tau2_lam, g0 = g[t], dx_n = dx_n, mean0 = initial$mean_lam,
                            a = priors$a$tau2_y, b = priors$b$tau2_y, mean = initial$mean_y, scale = initial$scale_y)
    
    llam_samples[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2

    # Update log llikelihood value for lambda layer
    obj_lam <-  logl_iso(llam_samples[t, ], dx_n, nugs = g[t], theta_lam[t, ], v = v, outer = initial$inner, 
                         calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, scale = initial$scale_lam, 
                         a = priors$a$tau2_lam, b = priors$b$tau2_lam)
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2
    
    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_store[t] <- tau2
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)     
  }
  
  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_samples, g = g,
              tau2 = tau2_store, tau2_lam = tau2_lam_store, llik_y = llik_y_store, llik_lam = llik_lam_store))
}

gibbs_sep_vdims <- function(YNs2 = NULL, Yn, xm, A, mcmc, initial, priors, v, vdims, reps_vdims, verb = TRUE){
  
  # mean layer
  nm <- nrow(xm)
  D <- ncol(xm)
  
  # variance layer
  xv <- reps_vdims$X0
  nv <- nrow(xv)
  d <- length(vdims)

  N <- sum(A)
  
  # Assign space and initial values
  llam_nv <- matrix(nrow = mcmc, ncol = nv) # only for unique nv X inputs
  theta_y <- matrix(nrow = mcmc, ncol = D)
  theta_lam <- matrix(nrow = mcmc, ncol = d)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_nv[1, ] <- initial$llam_nv
  
  llik_y_store <- llik_lam_store <- tau2_store <- tau2_lam_store <- rep(NA, mcmc)
  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- g <- NULL
  tau2 <- tau2_lam <- llam_n <- NULL

  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y!= 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # Map nv lambdas to n unique Xn 
    llam_n[reps_vdims$Z] <- rep(llam_nv[t - 1, ], reps_vdims$mult)
    
    # theta lam based on nv values
    for(i in 1:d){
      obj_lam <- mcmc_theta_l_sep(xv, llam_nv[t - 1, ], index = i, theta_cur_lam, llik_prev = llik_lam, 
                                  tau2_prev = tau2_lam, v = v, ls_check = initial$theta_check, 
                                  theta_y = theta_cur_y[vdims[i]], g0 = g[t - 1], alpha = priors$alpha$theta_lam, 
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
      obj_g <- mcmc_g_sep_inner(llam_nv[t - 1, ], xv, g[t - 1], llik_prev = llik_lam, tau2_prev = tau2_lam,
                                theta = theta_lam[t, ], v = v, outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                                mean = initial$mean_lam, scale = initial$scale_lam, alpha = priors$alpha$g, 
                                beta = priors$beta$g + t, l = priors$l , u = priors$u,
                                a = priors$a$tau2_lam, b = priors$b$tau2_lam)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
      
    }else g[t] <- initial$g
    
    # theta y updates
    for(i in 1:D){
      
      obj_y <- mcmc_theta_y_sep(YNs2 = YNs2, yn= Yn, xn= xm, A = A, index = i, llam_n, theta_cur_y, 
                                llik_prev = llik_y, tau2_prev = tau2, v = v,
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
    
    # Ess for log lam updates
    llam_draw <- ess_sample_vdims(YNs2= YNs2, yn=Yn, xn=xm, vdims = vdims, llam_nv[t - 1, ], A = A,
                                  llik_prev = llik_y, theta_y = theta_y[t, ], theta_lam = theta_lam[t, ], sep =T, 
                                  v = v, outer = initial$outer, calc_tau2 = initial$tau2, 
                                  scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam), 
                                  g0 = g[t], r0 = reps_vdims,
                                  a = priors$a$tau2_y, b = priors$b$tau2_y,
                                  mean0 = initial$mean_lam, mean = initial$mean_y, scale = initial$scale_y)
    
    llam_nv[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2
    
    obj_lam <- logl(llam_nv[t, ], xv, nugs = g[t], theta_lam[t, ], v= v, outer = initial$inner, 
                    calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, scale = initial$scale_lam, 
                    a = priors$a$tau2_lam, b = priors$b$tau2_lam)
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2

    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_store[t] <- tau2
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)     
  }
  
  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_nv, g = g,
              tau2 = tau2_store, tau2_lam = tau2_lam_store, llik_y = llik_y_store, llik_lam = llik_lam_store))
}

gibbs_iso_vdims <- function(YNs2 = NULL, Yn, xm, A, mcmc, initial, priors, v, vdims, reps_vdims, verb = TRUE){
  
  # Mean layer
  nm <- nrow(xm)
  D <- ncol(xm)
  
  xv <- reps_vdims$X0
  nv <- nrow(xv)
  d <- length(vdims)

  N <- sum(A)
  
  # Assign and initialise
  llam_nv <- matrix(nrow = mcmc, ncol = nv) # nv for unique inputs in X[ , vdims]
  theta_y <- matrix(nrow = mcmc, ncol = 1)
  theta_lam <- matrix(nrow = mcmc, ncol = 1)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_nv[1, ] <- initial$llam_nv
  
  llik_y_store <- llik_lam_store <- tau2_store <- tau2_lam_store <- rep(NA, mcmc)
  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- g <- NULL
  tau2 <- tau2_lam <- llam_n <- NULL
  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y!= 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g
  # Precal distance
  dxv <- sq_dist(xv)
  dxm <- sq_dist(xm)
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # Map nv unique llams to n inputs
    llam_n[reps_vdims$Z] <- rep(llam_nv[t - 1, ], reps_vdims$mult)
    
    # theta lam layer with nv X's
    obj_lam <- mcmc_theta_l(dxv, llam_nv[t-1, ], theta_cur_lam, llik_prev = llik_lam, tau2_prev = tau2_lam,
                            v = v, ls_check = initial$theta_check, theta_y = theta_cur_y, g0 = g[t - 1], 
                            alpha = priors$alpha$theta_lam, beta = priors$beta$theta_lam, l = priors$l , 
                            u = priors$u, inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                            a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam,
                            mean0 = initial$mean_lam, scale0 = initial$scale_lam)
    
    theta_lam[t] <- obj_lam$theta
    theta_cur_lam <- obj_lam$theta
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2
    
    if(initial$noise){
      obj_g <- mcmc_g_inner(llam_nv[t-1, ], dxv, g[t - 1], llik_prev = llik_lam, tau2_prev = tau2_lam,
                            theta = theta_lam[t], v = v, outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                            mean = initial$mean_lam, scale = initial$scale_lam, alpha = priors$alpha$g, 
                            beta = priors$beta$g + t, l = priors$l , u = priors$u,
                            a = priors$a$tau2_y, b = priors$b$tau2_y)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
    }else g[t] <- initial$g

    # theta y update
    obj_y <- mcmc_theta_y(YNs2, Yn, dxm, A, llam_n, theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, v = v, 
                          ls_check = initial$theta_check, theta_lam = theta_cur_lam,
                          outer = initial$outer, calc_tau2 = initial$tau2,
                          mean = initial$mean_y, scale = initial$scale_y, alpha = priors$alpha$theta_y, 
                          a = priors$a$tau2_y, b = priors$b$tau2_y,
                          beta = priors$beta$theta_y, l = priors$l, u = priors$u)
    theta_y[t] <- obj_y$theta
    theta_cur_y <- obj_y$theta
    llik_y <- obj_y$llik
    tau2 <- obj_y$tau2
    
    # Log lam ESS sampling
    llam_draw <- ess_sample_vdims(YNs2= YNs2, yn= Yn, xn= xm, vdims = vdims, llam_nv[t - 1, ], A,
                                  llik_prev = llik_y,  tau2_prev = tau2_lam, theta_y[t, ], theta_lam[t, ], 
                                  sep = FALSE, v = v, outer = initial$outer, calc_tau2 = initial$tau2,
                                  scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam), 
                                  g0 = g[t], r0 = reps_vdims, dx_n = dxm, dx_lam = dxv, 
                                  mean0 = initial$mean_lam, a = priors$a$tau2_y, b = priors$b$tau2_y,
                                  mean = initial$mean_y, scale = initial$scale_y)
    
    llam_nv[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2

    # lam layer likelihood update
    obj_lam <- logl_iso(llam_nv[t, ], dxv, nugs = g[t], theta_lam[t, ], v = v, outer = initial$inner, 
                        calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, scale = initial$scale_lam,
                        a = priors$a$tau2_lam, b = priors$b$tau2_lam)
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2

    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_store[t] <- tau2 
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)     
  }
  
  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_nv, g = g,
              tau2 = tau2_store, tau2_lam = tau2_lam_store, llik_y = llik_y_store, llik_lam = llik_lam_store))
}

gibbs_sep_hom <- function(YNs2 = NULL, Yn, Xn, A, mcmc, initial, priors, v, verb = TRUE){
  
  # Set up and initializing
  n <- nrow(Xn)
  D <- ncol(Xn)
  N <- sum(A)
  
  g_y <- matrix(nrow = mcmc, ncol = 1)
  theta_y <- matrix(nrow = mcmc, ncol = D)

  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  g_y[1, ] <- initial$g
  
  llik_y_store <- tau2_store <- rep(NA, mcmc)
  llik_y <- obj_y <- obj_g <- llam_draw <- tau2 <- NULL
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  tau2 <- tau2_store[1] <- initial$scale_y
  
  if(initial$outer) initial$scale_y <- 1
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)

    # Update theta y
    for(i in 1:D){
      
      obj_y <- mcmc_theta_sep_hom(YNs2 = YNs2, yn= Yn, xn=Xn, A = A, index = i, g_y[t - 1], 
                                theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, v = v,
                                outer = initial$outer, calc_tau2 = initial$tau2,
                                mean = initial$mean_y, scale = initial$scale_y,
                                a = priors$a$tau2_y, b = priors$b$tau2_y,
                                alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                                l = priors$l, u = priors$u)
      
      theta_y[t, ] <- obj_y$theta
      theta_cur_y <- obj_y$theta
      llik_y <- obj_y$llik
      tau2 <- obj_y$tau2
    }

    # Update g
    if(initial$noise & g_y[t-1] > 1e-8){
      obj_g <- mcmc_g_sep(YNs2 = YNs2, yn = Yn, xn = Xn, A = A, g_prev = g_y[t - 1], theta_y[t, ],
                          llik_prev = llik_y,  tau2_prev = tau2, v = v, outer = initial$outer, 
                          calc_tau2 = initial$tau2, alpha = priors$alpha$g, beta = priors$beta$g, 
                          l = priors$l, u = priors$u, a = priors$a$tau2_y, b = priors$b$tau2_y, 
                          mean = initial$mean_y, scale = initial$scale_y) 
      g_y[t] <- obj_g$g
      llik_y <- obj_g$llik
      tau2 <- obj_g$tau2
    } else if(initial$noise & g_y[t-1] < 1e-8) g_y[t] <- g_y[t-1]
    else g_y[t] <- initial$g
    
    llik_y_store[t] <- llik_y
    tau2_store[t] <- tau2
  }
  
  return(list(theta_y = theta_y, g_y = g_y, tau2 = tau2_store, llik_y = llik_y_store))
}

gibbs_iso_hom <- function(YNs2 = NULL, Yn, Xn, A, mcmc, initial, priors, v, verb = TRUE){  
  
  # Set up and intialization
  n <- nrow(Xn)
  D <- ncol(Xn)
  N <- sum(A)
  
  g_y <- matrix(nrow = mcmc, ncol = 1)
  theta_y <- matrix(nrow = mcmc, ncol = 1)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  g_y[1, ] <- initial$g
  
  llik_y_store <- tau2_store <- rep(NA, mcmc)
  llik_y <- obj_y <- obj_g <- tau2 <- NULL
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  tau2 <- tau2_store[1] <- initial$scale_y
  
  if(initial$outer) initial$scale_y <- 1
  
  # Precalculate distances
  dx_n <- sq_dist(Xn)
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # Update theta y
    obj_y <- mcmc_theta_hom(YNs2, Yn, dx_n, A, g_y[t - 1], theta_cur_y, llik_prev = llik_y,
                            tau2_prev = tau2, v = v, outer = initial$outer, calc_tau2 = initial$tau2,
                            mean = initial$mean_y, scale = initial$scale_y, alpha = priors$alpha$theta_y, 
                            beta = priors$beta$theta_y, l = priors$l, u = priors$u,
                            a = priors$a$tau2_y, b = priors$b$tau2_y)
    theta_y[t] <- obj_y$theta
    theta_cur_y <- obj_y$theta
    llik_y <- obj_y$llik
    tau2 <- obj_y$tau2
    
    # update nugget
    if(initial$noise){
      obj_g <- mcmc_g(YNs2= YNs2, yn=Yn, dx_n, g_y[t - 1], A = A, theta_y = theta_y[t, ], 
                      llik_prev = llik_y, tau2_prev = tau2, v = v, outer = initial$outer, 
                      calc_tau2 = initial$tau2, mean = initial$mean_y, scale = initial$scale_y,
                      alpha = priors$alpha$g, beta = priors$beta$g, l = priors$l, u = priors$u,
                      a = priors$a$tau2_y, b = priors$b$tau2_y)
      
      g_y[t] <- obj_g$g
      llik_y <- obj_g$llik 
      tau2 <- obj_g$tau2
    } else g_y[t] <- initial$g
  
    llik_y_store[t] <- llik_y
    tau2_store[t] <- tau2
  }
  
  return(list(theta_y = theta_y, g_y = g_y, tau2 = tau2_store, llik_y = llik_y_store))
}

