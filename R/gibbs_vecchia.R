# ---------------  Gibbs: HetGP + HomGP  Vecchia -------------------------------
# gibbs_sep_vec: Gibbs sampler for separable HetGP w/ vec
# gibbs_iso_vec: Gibbs sampler for isotropic HetGP w/ vec
# gibbs_sep_vdims_vec: Gibbs for Sep HetGP where variance changed in some dims w/ vec
# gibbs_iso_vdims_vec: Gibbs for Iso HetGP where variance changed in some dims w/ vec
# gibbs_sep_vec_hom: Gibbs for Sep HomGP w/ vecchia
# gibbs_iso_vec_hom: Gibbs for Iso HomGP w/ vecchia
# ------------------------------------------------------------------------------

gibbs_sep_vec <- function(YNs2 = NULL, Yn, x_approx, A, mcmc, initial, priors, v, verb = TRUE){
  
  # Initial setup
  n <- nrow(x_approx$x_ord)
  D <- ncol(x_approx$x_ord)
  N <- sum(A)
  
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
  if(initial$scale_y != 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g

  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # theta lam update
    for(i in 1:D){
      obj_lam <- mcmc_theta_l_sep_vec(x_approx, llam_samples[t-1, ], index = i, theta_cur_lam, 
                                      llik_prev = llik_lam, tau2_prev = tau2_lam, v = v,
                                      ls_check = initial$theta_check, theta_y = theta_cur_y[i], 
                                      g0 = g[t - 1], alpha = priors$alpha$theta_lam, 
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
      obj_g <- mcmc_g_vec_inner(llam_samples[t-1, ], x_approx, g[t - 1], llik_prev = llik_lam, 
                                tau2_prev = tau2_lam, theta = theta_lam[t, ], v = v, outer = initial$inner, 
                                calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, 
                                scale = initial$scale_lam, alpha = priors$alpha$g, beta = priors$beta$g + (t),
                                l = priors$l , u = priors$u, a = priors$a$tau2_lam, b = priors$b$tau2_lam, 
                                sep = TRUE)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
    } else g[t] <- initial$g

    # theta y updates
    for(i in 1:D){
      obj_y <- mcmc_theta_y_sep_vec(YNs2 = YNs2, yn= Yn, x_approx, A, index = i, llam_samples[t-1, ],
                                    theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, v = v, 
                                    ls_check = initial$theta_check, theta_lam = theta_cur_lam[i],
                                    outer = initial$outer, calc_tau2 = initial$tau2,
                                    mean = initial$mean_y, scale = initial$scale_y,
                                    alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                                    l = priors$l , u = priors$u, a = priors$a$tau2_y, b = priors$b$tau2_y)
      theta_y[t, ] <- obj_y$theta
      theta_cur_y <- obj_y$theta
      llik_y <- obj_y$llik
      tau2 <- obj_y$tau2
    }
    
    # ESS sample for lambda
    llam_draw <- ess_sample_vec(YNs2, yn = Yn, x_approx = x_approx, llam_samples[t - 1, ], A = A,
                                llik_prev = llik_y, theta_y = theta_y[t, ], theta_lam = theta_lam[t, ], 
                                sep = TRUE, v = v, outer = initial$outer, calc_tau2 = initial$tau2, 
                                mean = initial$mean_y, scale = initial$scale_y, mean0 = initial$mean_lam, 
                                scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam),
                                g0 = g[t], a = priors$a$tau2_y, b = priors$b$tau2_y)
    
    llam_samples[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2

    # llik lambda update
    obj_lam <-  logl_vec(llam_samples[t, ], x_approx, nugs = g[t], theta_lam[t, ],
                         outer = initial$inner, v = v, calc_tau2 = initial$inner_tau2, 
                         mean = initial$mean_lam, sep = TRUE, scale = initial$scale_lam, 
                         a = priors$a$tau2_lam, b= priors$b$tau2_lam, latent = TRUE)
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

gibbs_iso_vec <- function(YNs2 = NULL, Yn, x_approx, A, mcmc, 
                          initial, priors, v, reps_vdims, verb = TRUE){

  # Initial Setup
  n <- nrow(x_approx$x_ord)
  N <- sum(A)

  llam_samples <- matrix(nrow = mcmc, ncol = n)
  theta_y <- matrix(nrow = mcmc, ncol = 1)
  theta_lam <- matrix(nrow = mcmc, ncol = 1)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_samples[1, ] <- initial$llam
  
  llik_y_store <- llik_lam_store <- tau2_lam_store <- tau2_store <- rep(NA, mcmc)
  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- g <- NULL
  tau2 <- tau2_lam <- NULL
  
  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y != 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # Update theta lambda
    obj_lam <-  mcmc_theta_l_vec(x_approx, llam_samples[t-1, ], theta_cur_lam, 
                                 llik_prev = llik_lam, tau2_prev = tau2_lam, v = v,
                                 ls_check = initial$theta_check, theta_y = theta_cur_y,
                                 g0 = g[t-1], alpha = priors$alpha$theta_lam, 
                                 beta = priors$beta$theta_lam, l = priors$l , u = priors$u,
                                 inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                                 mean0 = initial$mean_lam, scale0 = initial$scale_lam,
                                 a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam)
    
    theta_lam[t] <- obj_lam$theta
    theta_cur_lam <- obj_lam$theta
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2

    if(initial$noise){
      obj_g <- mcmc_g_vec_inner(llam_samples[t-1, ], x_approx, g[t - 1], llik_prev = llik_lam,
                                tau2_prev = tau2_lam, theta = theta_lam[t, ], v = v, outer = initial$inner, 
                                calc_tau2 = initial$inner_tau2, mean = initial$mean_lam, 
                                scale = initial$scale_lam, alpha = priors$alpha$g, beta = priors$beta$g + t,
                                l = priors$l, u = priors$u, a = priors$a$tau2_lam, b = priors$b$tau2_lam, sep = TRUE)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
      
    }else g[t] <- initial$g
    
    # Theta y update
    obj_y <- mcmc_theta_y_vec(YNs2, yn = Yn, x_approx = x_approx, A = A, llam_samples[t-1, ],
                              theta_prev = theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, v= v,
                              ls_check = initial$theta_check, theta_lam = theta_cur_lam,
                              outer = initial$outer, calc_tau2 = initial$tau2,
                              mean = initial$mean_y, scale = initial$scale_y,
                              alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                              a = priors$a$tau2_y, b = priors$b$tau2_y,
                              l = priors$l , u = priors$u)
    
    theta_y[t] <- obj_y$theta
    theta_cur_y <- obj_y$theta
    llik_y <- obj_y$llik
    tau2 <- obj_y$tau2
    
    # ESS for log lambda
    llam_draw <- ess_sample_vec(YNs2, yn = Yn, x_approx = x_approx, llam_samples[t - 1, ], A,
                                llik_prev = llik_y, theta_y = theta_y[t], theta_lam = theta_lam[t], 
                                sep = FALSE, v= v, outer = initial$outer, calc_tau2 = initial$tau2,
                                mean = initial$mean_y, scale = initial$scale_y,
                                mean0 = initial$mean_lam, 
                                scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam),
                                g0 = g[t],
                                a = priors$a$tau2_y, b = priors$b$tau2_y)
    
    llam_samples[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2

    # Update llik for theta after new lambdas drawn
    obj_lam <-  logl_vec(llam_samples[t, ], x_approx, nugs = g[t], theta_lam[t, ],
                          outer = initial$inner, v = v, calc_tau2 = initial$inner_tau2, sep = FALSE, 
                          mean = initial$mean_lam, scale = initial$scale_lam, 
                         a = priors$a$tau2_lam, b = priors$b$tau2_lam, latent = TRUE)
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

# Change latent arg

gibbs_sep_vdims_vec <- function(YNs2 = NULL, Yn, xm_approx, A, vdims, mcmc, 
                                initial, priors, v, reps_vdims, 
                                xv_approx = NULL, vecchia_var, verb = TRUE){
  
  # initial setup - mean layer
  nm <- nrow(xm_approx$x_ord)
  D <- ncol(xm_approx$x_ord)
  
  # Variance layer
  xv <- reps_vdims$X0
  nv <- nrow(xv)
  d <- length(vdims)

  N <- sum(A)
  
  # Check if vecchia needed for variance layer
  llam_nv <- matrix(nrow = mcmc, ncol = nv)
  theta_y <- matrix(nrow = mcmc, ncol = D)
  theta_lam <- matrix(nrow = mcmc, ncol = d)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_nv[1, ] <- initial$llam_nv
  
  llik_y_store <- llik_lam_store <- tau2_lam_store <- tau2_store <- rep(NA, mcmc)
  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw <- g <- NULL
  tau2 <- tau2_lam <- llam_n <- NULL

  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y != 1) stop("broken") # reset scale
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  initial$scale_y <- initial$scale_lam <- 1

  g[1] <- initial$g
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    # if(t %% 10 == 0 ) print(c(t, tau2_lam))
    
    # Map nv unique lambdas to n unique inputs
    llam_n[reps_vdims$Z] <- rep(llam_nv[t - 1, ],  reps_vdims$mult)
    # if(tau2_lam > 500) stop("check")
    for(i in 1:d){
      
      # theta lam update 
      if(vecchia_var){ # using vecchia
        obj_lam <- mcmc_theta_l_sep_vec(xv_approx, llam_nv[t - 1, ], index = i, 
                                        theta_cur_lam, llik_prev = llik_lam, tau2_prev = tau2_lam, v = v,
                                        ls_check = initial$theta_check, theta_y = theta_cur_y[vdims[i]], 
                                        g0 = g[t - 1],  alpha = priors$alpha$theta_lam, 
                                        beta = priors$beta$theta_lam, l = priors$l , u = priors$u,
                                        inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                                        mean0 = initial$mean_lam, scale0 = initial$scale_lam,
                                        a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam) 
      } else{ # no vecchia
        obj_lam <- mcmc_theta_l_sep(xv, llam_nv[t - 1, ], index = i, theta_cur_lam, llik_prev = llik_lam, 
                                    tau2_prev = tau2_lam, v = v, ls_check = initial$theta_check, 
                                    theta_y = theta_cur_y[vdims[i]], g0 = g[t - 1],  alpha = priors$alpha$theta_lam, 
                                    beta = priors$beta$theta_lam, l = priors$l , u = priors$u,
                                    inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                                    mean0 = initial$mean_lam, scale0 = initial$scale_lam,
                                    a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam)
      }
      
      theta_lam[t, ] <- obj_lam$theta
      theta_cur_lam <- obj_lam$theta
      llik_lam <- obj_lam$llik
      tau2_lam <- obj_lam$tau2
    }
    
    if(initial$noise){
      if(vecchia_var)
        obj_g <- mcmc_g_vec_inner(llam_nv[t - 1, ], xv_approx, g[t - 1], llik_prev = llik_lam, 
                                  tau2_prev = tau2_lam, theta = theta_lam[t, ],
                                  v = v, outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                                  mean = initial$mean_lam, scale = initial$scale_lam,
                                  alpha = priors$alpha$g, beta = priors$beta$g + (t),
                                  l = priors$l , u = priors$u, sep = TRUE,
                                  a = priors$a$tau2_lam, b = priors$b$tau2_lam)
      else
        obj_g <- mcmc_g_sep_inner(llam_nv[t - 1, ], xv, g[t - 1], llik_prev = llik_lam,
                              tau2_prev = tau2_lam, theta = theta_lam[t, ], v = v, 
                              outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                              mean = initial$mean_lam, scale = initial$scale_lam,
                              alpha = priors$alpha$g, beta = priors$beta$g + (t),
                              l = priors$l , u = priors$u, 
                              a = priors$a$tau2_lam, b = priors$b$tau2_lam)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
    }else g[t] <- initial$g
    
    # theta y update
    for(i in 1:D){ 
      obj_y <- mcmc_theta_y_sep_vec(YNs2 = YNs2, yn= Yn, xm_approx, A, index = i, llam_n, theta_prev = theta_cur_y, 
                                    llik_prev = llik_y, tau2_prev = tau2, v = v, ls_check = initial$theta_check, 
                                    theta_lam = ifel(any(i == vdims), theta_cur_lam[which(i == vdims)], NULL),
                                    outer = TRUE, calc_tau2 = TRUE, mean = initial$mean_y, scale = initial$scale_y, 
                                    alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                                    l = priors$l , u = priors$u, a = priors$a$tau2_y, b = priors$b$tau2_y)
      theta_y[t, ] <- obj_y$theta
      theta_cur_y <- obj_y$theta
      llik_y <- obj_y$llik
      tau2 <- obj_y$tau2
    }
    
    # ess sampling for new lambdas
    llam_draw <- ess_sample_vec_vdims(YNs2, yn = Yn, xm_approx, xv_approx, vdims, llam_nv[t - 1, ], A,
                                      llik_prev = llik_y, theta_y = theta_y[t, ], theta_lam = theta_lam[t, ], 
                                      sep = TRUE, v = v, outer = TRUE, calc_tau2 = TRUE, mean = initial$mean_y, 
                                      scale = initial$scale_y, mean0 = initial$mean_lam, 
                                      scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam), 
                                      g0 = g[t], r0 = reps_vdims, vec_var = vecchia_var,
                                      a = priors$a$tau2_y, b = priors$b$tau2_y)
    
    llam_nv[t, ] <- llam_draw$llam 
    llik_y <- llam_draw$llik # variance process llik
    tau2 <- llam_draw$tau2

    # calculate log llik for lambda's 
    if(vecchia_var) 
      obj_lam <-  logl_vec(llam_nv[t, ], xv_approx, nugs = g[t], theta_lam[t, ],
                            outer = initial$inner, v = v, calc_tau2 = initial$inner_tau2,
                            mean = initial$mean_lam, scale = initial$scale_lam, 
                           a = priors$a$tau2_lam, b = priors$b$tau2_lam, sep = TRUE, latent = TRUE)
    else
      obj_lam <-  logl(llam_nv[t, ], xv, nugs = g[t], theta_lam[t, ],
                       outer = initial$inner, v = v, calc_tau2 = initial$inner_tau2, 
                       mean = initial$mean_lam, scale = initial$scale_lam, 
                       a = priors$a$tau2_lam, b = priors$b$tau2_lam)
    
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2

    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_store[t] <- tau2
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)
  }

  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_nv, tau2 = tau2_store,
              tau2_lam = tau2_lam_store, g = g, xv = xv, xv_approx = xv_approx, llik_y = llik_y_store, 
              llik_lam = llik_lam_store))
}

gibbs_iso_vdims_vec <- function(YNs2 = NULL, Yn, xm_approx, A, vdims, mcmc, 
                                initial, priors, v, reps_vdims, 
                                xv_approx = NULL, vecchia_var, verb = TRUE){
  
  # Setup - Mean layer
  nm <- nrow(xm_approx$x_ord)
  D <- ncol(xm_approx$x_ord)
  
  xv <- reps_vdims$X0
  nv <- nrow(xv)
  d <- length(vdims)

  N <- sum(A)
  
  # check if vec needed for var layer
  if(!vecchia_var) dxv <- sq_dist(xv)
  
  llam_nv <- matrix(nrow = mcmc, ncol = nv)
  theta_y <- matrix(nrow = mcmc, ncol = 1)
  theta_lam <- matrix(nrow = mcmc, ncol = 1)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  theta_cur_lam <- theta_lam[1, ] <- initial$theta_lam
  llam_nv[1, ] <- initial$llam_nv
  
  llik_y_store <- llik_lam_store <- tau2_store <- tau2_lam_store <- rep(NA, mcmc)
  llik_y <- llik_lam  <- obj_y <- obj_lam <- llam_draw  <- g <-  NULL
  tau2 <- llam_n <- tau2_lam <- NULL

  tau2_lam <- tau2_lam_store[1] <- initial$tau2_lam
  tau2 <- tau2_store[1] <- initial$tau2_y
  
  if(initial$inner_tau2 & initial$scale_lam != 1) stop("broken") # reset scale for prof ll
  if(initial$scale_y != 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  llik_lam <- llik_lam_store[1] <- initial$llik_lam
  
  g[1] <- initial$g
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # Map nv lambdas from X[ , vdims] to n unique inputs Xn
    llam_n[reps_vdims$Z] <- rep(llam_nv[t - 1, ], reps_vdims$mult)
    
    # update theta lam
    if(vecchia_var){ # if vecchia
      obj_lam <-  mcmc_theta_l_vec(xv_approx, llam_nv[t - 1, ], theta_cur_lam, llik_prev = llik_lam, 
                                   tau2_prev = tau2_lam, v = v, ls_check = initial$theta_check, 
                                   theta_y = theta_cur_y, g0 = g[t - 1], alpha = priors$alpha$theta_lam, 
                                   beta = priors$beta$theta_lam, l = priors$l , u = priors$u,
                                   inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                                   mean0 = initial$mean_lam, scale0 = initial$scale_lam,
                                   a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam)
    } else{ # if non vec
      obj_lam <-  mcmc_theta_l(dxv, llam_nv[t - 1, ], theta_cur_lam, llik_prev = llik_lam,
                               tau2_prev = tau2_lam, v = v, ls_check = initial$theta_check, 
                               theta_y = theta_cur_y, g0 = g[t - 1], alpha = priors$alpha$theta_lam, 
                               beta = priors$beta$theta_lam, l = priors$l , u = priors$u,
                               inner = initial$inner, calc_inner_tau2 = initial$inner_tau2,
                               mean0 = initial$mean_lam, scale0 = initial$scale_lam,
                               a0 = priors$a$tau2_lam, b0 = priors$b$tau2_lam)
    }
    
    theta_lam[t] <- obj_lam$theta
    theta_cur_lam <- obj_lam$theta
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2
    
    if(initial$noise){
      if(vecchia_var)
        obj_g <- mcmc_g_vec_inner(llam_nv[t - 1, ], xv_approx, g[t - 1], llik_prev = llik_lam,
                                  tau2_prev = tau2_lam, theta = theta_lam[t, ], v = v, 
                                  outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                                  mean = initial$mean_lam, scale = initial$scale_lam,
                                  alpha = priors$alpha$g, beta = priors$beta$g + (t),
                                  l = priors$l , u = priors$u, sep = FALSE,
                                  a = priors$a$tau2_lam, b = priors$b$tau2_lam)
      else 
        obj_g <- mcmc_g_inner(llam_nv[t - 1, ], dxv, g[t - 1], llik_prev = llik_lam,
                              tau2_prev = tau2_lam, theta = theta_lam[t, ], v = v, 
                              outer = initial$inner, calc_tau2 = initial$inner_tau2, 
                              mean = initial$mean_lam, scale = initial$scale_lam,
                              alpha = priors$alpha$g, beta = priors$beta$g + (t),
                              l = priors$l , u = priors$u, 
                              a = priors$a$tau2_lam, b = priors$b$tau2_lam)
      
      g[t] <- obj_g$g
      llik_lam <- obj_g$llik
      tau2_lam <- obj_g$tau2
    }else g[t] <- initial$g
    

    # update theta y
    obj_y <- mcmc_theta_y_vec(YNs2, yn = Yn, x_approx = xm_approx, A = A, llam_n,
                              theta_prev = theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, v = v,
                              ls_check = initial$theta_check, theta_lam = theta_cur_lam,
                              outer = initial$outer, calc_tau2 = initial$tau2, 
                              mean = initial$mean_y, scale = initial$scale_y, 
                              alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                              a = priors$a$tau2_y, b = priors$b$tau2_y,
                              l = priors$l , u = priors$u)
    
    theta_y[t] <- obj_y$theta
    theta_cur_y <- obj_y$theta
    llik_y <- obj_y$llik
    tau2 <- obj_y$tau2
    
    # ESS sample for log lambdas
    llam_draw <- ess_sample_vec_vdims(YNs2, yn = Yn, x_approx = xm_approx, xv_approx, vdims, llam_nv[t - 1, ], A,
                                      llik_prev = llik_y, theta_y = theta_y[t], theta_lam = theta_lam[t], 
                                      sep = FALSE, v= v, outer = initial$outer, calc_tau2 = initial$tau2, 
                                      mean = initial$mean_y, scale = initial$scale_y, mean0 = initial$mean_lam,
                                      scale0 = ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam), g0 = g[t], a = priors$a$tau2_y, 
                                      b = priors$b$tau2_y, r0 = reps_vdims, vec_var = vecchia_var)
    
    llam_nv[t, ] <- llam_draw$llam
    llik_y <- llam_draw$llik
    tau2 <- llam_draw$tau2
    
    # recalc likelihood for lambda layer
    if(vecchia_var)
      obj_lam <-  logl_vec(llam_nv[t, ], xv_approx, nugs = g[t], theta_lam[t, ], outer = initial$inner, 
                           v = v, calc_tau2 = initial$inner_tau2, sep = FALSE, 
                           mean = initial$mean_lam, scale = initial$scale_lam,
                           a = priors$a$tau2_lam, b = priors$b$tau2_lam, latent = TRUE)
    else 
      obj_lam <-  logl_iso(llam_nv[t, ], dxv, nugs = g[t], theta_lam[t, ], outer = initial$inner, 
                           v = v, calc_tau2 = initial$inner_tau2,
                           mean = initial$mean_lam,  scale = initial$scale_lam,
                           a = priors$a$tau2_lam, b = priors$b$tau2_lam)
    
    llik_lam <- obj_lam$llik
    tau2_lam <- obj_lam$tau2

    llik_y_store[t] <- llik_y
    llik_lam_store[t] <- llik_lam
    
    tau2_store[t] <- tau2
    tau2_lam_store[t] <- ifel(initial$prof_ll_lam, tau2_lam, initial$scale_lam)
  }
  
  return(list(theta_lam = theta_lam, theta_y = theta_y, llam_samples = llam_nv, tau2 = tau2_store, 
              tau2_lam = tau2_lam_store, g = g, xv = xv, xv_approx = xv_approx, llik_y = llik_y_store, 
              llik_lam = llik_lam_store))
}

gibbs_sep_vec_hom <- function(YNs2= NULL, Yn, x_approx, A, mcmc, initial, priors, v, verb = TRUE){
  
  # Initial setup
  n <- nrow(x_approx$x_ord)
  D <- ncol(x_approx$x_ord)
  N <- sum(A)
  
  g_y <- matrix(nrow = mcmc, ncol = 1)
  theta_y <- matrix(nrow = mcmc, ncol = D)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  g_y[1, ] <- initial$g
  
  llik_y_store <- tau2_store <- rep(NA, mcmc)
  llik_y  <- obj_y <- obj_g <- tau2 <- NULL
  
  tau2 <- tau2_store[1] <- initial$tau2_y
  if(initial$tau2 == TRUE & initial$scale != 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # theta y updates
    for(i in 1:D){
      obj_y <- mcmc_theta_sep_vec_hom(YNs2 = YNs2, yn= Yn, x_approx, A, index = i, g_y[t-1, ],
                                      theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, v = v, 
                                      outer = initial$outer, calc_tau2 = initial$tau2, 
                                      mean = initial$mean_y, scale = initial$scale_y,
                                      alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                                      l = priors$l , u = priors$u, a = priors$a$tau2_y, b = priors$b$tau2_y)
      theta_y[t, ] <- obj_y$theta
      theta_cur_y <- obj_y$theta
      llik_y <- obj_y$llik
      tau2 <- obj_y$tau2
    }
    
    # nugget update
    if(initial$noise & g_y[t-1] > 1e-8){
      obj_g <- mcmc_g_vec(YNs2, yn = Yn, x_approx, A, g_y[t - 1, ], llik_prev = llik_y, tau2_prev = tau2, 
                          theta_y = theta_y[t, ], v = v, outer = initial$outer, calc_tau2 = initial$tau2,
                          mean = initial$mean_y, scale = initial$scale_y,
                          alpha = priors$alpha$g, beta = priors$beta$g, 
                          l = priors$l , u = priors$u, a = priors$a$tau2_y, b = priors$b$tau2_y, sep = TRUE)
      
      g_y[t, ] <- obj_g$g
      llik_y <- obj_g$llik
      tau2 <- obj_g$tau2
    } else if(initial$noise & g_y[t-1] <= 1e-8) g_y[t] <- g_y[t-1]
    else g_y[t] <- initial$g
    
    llik_y_store[t] <- llik_y
    tau2_store[t] <- tau2
  }
  
  return(list(theta_y = theta_y, g_y = g_y, tau2 = tau2_store, llik_y = llik_y_store))
}

gibbs_iso_vec_hom <- function(YNs2 = NULL, Yn, x_approx, A, mcmc, initial, priors, v, verb = TRUE){
  
  # Initial setup
  n <- nrow(x_approx$x_ord)
  D <- ncol(x_approx$x_ord)
  N <- sum(A)
  
  g_y <- matrix(nrow = mcmc, ncol = 1)
  theta_y <- matrix(nrow = mcmc, ncol = D)
  
  theta_cur_y <- theta_y[1, ] <- initial$theta_y
  g_y[1, ] <- initial$g
  
  llik_y_store <- tau2_store <- rep(NA, mcmc)
  llik_y  <- obj_y <- obj_g <- tau2 <- NULL
  
  tau2 <- tau2_store[1] <- initial$tau2_y
  if(initial$tau2 == TRUE & initial$scale != 1) stop("broken") # reset scale 
  
  llik_y <- llik_y_store[1] <- initial$llik_y
  
  for (t in 2:mcmc) {
    
    if(t %% 500 == 0 & verb == TRUE) print(t)
    
    # update theta y
    obj_y <- mcmc_theta_vec_hom(YNs2, yn = Yn, x_approx = x_approx, A = A, g_y[t - 1, ],
                              theta_prev = theta_cur_y, llik_prev = llik_y, tau2_prev = tau2, 
                              v= v, outer = initial$outer, calc_tau2 = initial$tau2,
                              mean = initial$mean_y, scale = initial$scale_y,
                              alpha = priors$alpha$theta_y, beta = priors$beta$theta_y, 
                              l = priors$l , u = priors$u, a = priors$a$tau2_y, b = priors$b$tau2_y)
    
    theta_y[t] <- obj_y$theta
    theta_cur_y <- obj_y$theta
    llik_y <- obj_y$llik
    tau2 <- obj_y$tau2
    
    # update nugget
    if(initial$noise){
      obj_g <- mcmc_g_vec(YNs2, yn = Yn, x_approx, A, g_y[t - 1, ], llik_prev = llik_y, 
                          tau2_prev = tau2, theta_y = theta_y[t, ], v = v, outer = initial$outer, 
                          calc_tau2 = initial$tau2, mean = initial$mean_y, scale = initial$scale_y,
                          alpha = priors$alpha$g, beta = priors$beta$g, 
                          l = priors$l , u = priors$u, a = priors$a$tau2_y, b = priors$b$tau2_y, sep = FALSE)
      
      g_y[t, ] <- obj_g$g
      llik_y <- obj_g$llik
      tau2 <- obj_g$tau2
    } else g_y[t] <- initial$g
    
    llik_y_store[t] <- llik_y 
    tau2_store[t] <- tau2
  }
  
  return(list(theta_y = theta_y, g_y = g_y, tau2 = tau2_store, llik_y = llik_y_store))
}

# ------------------------------------------------------------------------------