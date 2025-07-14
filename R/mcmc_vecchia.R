# ----------------------  MCMC HomGP + HetGP Vecchia ---------------------------

# ------- HetGP Vecchia --------------------------------------------------------
# ess_sample_vec: ESS sampling for Lambda layer w/ vecchia
# mcmc_theta_y_sep_vec: MCMC sampling for theta y (sep case)
# mcmc_theta_l_sep_vec: MCMC sampling for theta lam (sep case)
# mcmc_theta_y_vec: MCMC sampling for theta y (isotropic)
# mcmc_theta_l_vec: MCMC sampling for theta l (isotropic)
# ess_sample_vec_vdims: ESS for Lambda layer when lambdas change in some dims

# ------- HomGP Vecchia --------------------------------------------------------
# mcmc_g_vec: MCMC sampling for g isotropic w/ vecchia 
# mcmc_theta_sep_vec_hom: MCMC sampling for theta y sep w/ vecchia 
# mcmc_theta_vec_hom: MCMC sampling for theta y isotropic w/vecchia
# ------------------------------------------------------------------------------

# -------------------- HetGP ---------------------------------------------------
ess_sample_vec <- function(YNs2 = NULL, yn, x_approx, llam_prev, A, llik_prev = NULL,
                           theta_y, theta_lam, sep = TRUE, v, outer, calc_tau2, 
                           mean, scale, mean0, scale0, g0, a, b) {
  # -----------------
  
  # Prev Likelihood
  if(is.null(llik_prev)){ # Woodbury
    if(!is.null(YNs2)){
      obj_prev <- loglw_vec(YNs2, yn, x_approx, A= A, exp(llam_prev), 
                            theta=theta_y, outer = outer, v = v, calc_tau2 = calc_tau2, 
                            sep = sep,  mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
    } else{ # Full N
      obj_prev <- logl_vec(yn, x_approx, nugs = exp(llam_prev), theta=theta_y, 
                           outer = outer, v = v, calc_tau2 = calc_tau2, sep = sep, 
                           mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
    }
  }  

  # Propose a prior
  llam_prior <- rand_mvn_vec(x_approx, theta_lam, v = v, mean = mean0,
                             g = rep(g0, nrow(x_approx$x_ord)), scale = scale0, sep = sep)

  # angle
  gam <- runif(1,0,2*pi)
  
  # bounds
  gam_min <- gam - 2*pi
  gam_max <- gam
  
  count=0 # If ESS takes too long
  while (1){
    
    count = count + 1
    if(count > 50) message("too many iterations")
    
    # New log-lambda 
    llam_new <- llam_prev * cos (gam) + llam_prior * sin (gam)
    
    # l_new <- exp(llam_new)
    # l_new[l_new == 0] <- 1e-8
    
    # Woodbury new llik
    if(!is.null(YNs2)) {
      obj_prop <- loglw_vec(YNs2, yn, x_approx, A = A, exp(llam_new), 
                            theta=theta_y, outer = outer, v = v, calc_tau2 = calc_tau2, 
                            sep = sep, mean = mean, scale = scale, a= a, b = b)
      llik_prop <- obj_prop$llik
      tau2 <- obj_prop$tau2
    } else{ # Full N
      obj_prop <- logl_vec(yn, x_approx, nugs = exp(llam_new), theta=theta_y, outer = outer, v = v, 
                           calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
      llik_prop <- obj_prop$llik
      tau2 <- obj_prop$tau2
    }
    
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

mcmc_theta_y_sep_vec <- function(YNs2= NULL, yn= NULL, x_approx, A= NULL, index, llam_prev, theta_prev, 
                                 llik_prev = NULL, tau2_prev = NULL, v, ls_check = FALSE, theta_lam = NULL, 
                                 outer, calc_tau2, mean , scale, alpha, beta, l , u, a, b){
  
  theta_star <- theta_prev
  theta_star[index] <- runif(1, min = l * theta_prev[index]/u , max = u * theta_prev[index]/l)
  
  # check if variance layer theta > mean theta
  if(ls_check == TRUE & !is.null(theta_lam)){
      if(theta_star[index] > theta_lam) 
        return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
    } 
  
  # Woodbury llik
  if(!is.null(YNs2)){
    if(is.null(llik_prev)) {
      obj_prev <- loglw_vec(YNs2, yn, x_approx, A, exp(llam_prev), theta = theta_prev, outer = outer, v = v, 
                            calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    # else{
    #   obj_prev <- loglw_vec(YNs2, yn, x_approx, A, exp(llam_prev), theta = theta_prev, outer = outer, v = v,
    #                         calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
    #   if(abs(llik_prev - obj_prev$llik) > 1e-8) stop('llik mismatch')
    #   if(abs(tau2_prev - obj_prev$tau2) > 1e-8) stop("tau2 mismatch")
    # }
    obj_prop <- loglw_vec(YNs2, yn, x_approx, A, exp(llam_prev), theta = theta_star, outer = outer, v = v, 
                          calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  } else{ # full N
    if(is.null(llik_prev)) {
      obj_prev <- logl_vec(out_vec= yn, approx = x_approx, nugs = exp(llam_prev), theta= theta_prev, outer = outer, 
                           v = v, calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    obj_prop <- logl_vec(out_vec= yn, approx = x_approx, nugs = exp(llam_prev), theta= theta_star, outer = outer, 
                         v = v, calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  
  # Priors
  prior_prev <- dgamma(theta_prev[index], shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(theta_star[index], shape = alpha, rate = beta, log = TRUE)
  
  # Accept/ Reject
  l_ratio <- llik_prop + prior_prop - llik_prev - prior_prev 
  rn <- runif(1, 0, 1)

  threshold <- l_ratio + log(theta_prev[index]) - log(theta_star[index]) 
  
  if(log(rn) < threshold) 
    return(list(theta = theta_star, llik = llik_prop, tau2 = tau2_prop))
  else
    return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
}

mcmc_theta_l_sep_vec <- function(x_approx, llam_prev, index, theta_prev, llik_prev = NULL, tau2_prev = NULL,
                                 v, ls_check = FALSE, theta_y, g0, alpha, beta, l , u, inner, calc_inner_tau2, 
                                 mean0, scale0 = 1, a0, b0){
  
  theta_star <- theta_prev
  theta_star[index] <- runif(1, min = l * theta_prev[index]/u , max = u * theta_prev[index]/l)
  
  # check if variance layer theta > mean theta
  if(ls_check == TRUE)
    if(theta_star[index] < theta_y) 
      return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
  
  # Likelihood
  if(is.null(llik_prev)) {
    obj_prev <- logl_vec(llam_prev, x_approx, nugs= g0, theta = theta_prev, outer = inner, v= v, 
                         calc_tau2 = calc_inner_tau2, sep = TRUE, mean = mean0, 
                         scale = scale0, a = a0, b = b0, latent = TRUE)
    llik_prev <- obj_prev$llik
    tau2_prev <- obj_prev$tau2
  }
  # else{
  #   obj_prev <- logl_vec(llam_prev, x_approx, nugs= g0, theta = theta_prev, outer = inner, v= v, 
  #                        calc_tau2 = calc_inner_tau2, sep = TRUE, mean = mean0, 
  #                        scale = scale0, a = a0, b = b0, latent = T)
  #   if(abs(llik_prev - obj_prev$llik) > 1e-8) stop('llik mismatch')
  #   if(abs(tau2_prev - obj_prev$tau2) > 1e-8) stop("tau2 mismatch")
  # }
  obj_prop <- logl_vec(llam_prev, x_approx, nugs= g0, theta = theta_star, outer = inner,
                       v= v, calc_tau2 = calc_inner_tau2, sep = TRUE, mean = mean0, scale = scale0,
                       a = a0, b = b0, latent = TRUE)
  llik_prop <- obj_prop$llik
  tau2_prop <- obj_prop$tau2
  
  # Priors
  prior_prev <- dgamma(theta_prev[index], shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(theta_star[index], shape = alpha, rate = beta, log = TRUE)
  
  # Accept/Reject
  l_ratio <- llik_prop + prior_prop - llik_prev - prior_prev 
  rn <- runif(1, 0, 1)
  # r <- exp(l_ratio) * (theta_prev[index]/theta_star[index])
  
  threshold <- l_ratio + log(theta_prev[index]) - log(theta_star[index]) 
  
  # if(runif(1, 0, 1) < exp(l_ratio) * (theta_prev[index]/theta_star[index]))
  if(log(rn) < threshold)
    return(list(theta = theta_star, llik = llik_prop, tau2 = tau2_prop))
  else
    return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
}

mcmc_theta_y_vec <- function(YNs2= NULL, yn= NULL, x_approx, A=NULL, llam_prev, theta_prev, llik_prev = NULL, 
                             tau2_prev = NULL, v, ls_check = FALSE, theta_lam, outer, calc_tau2, 
                             mean, scale, a, b, alpha, beta, l , u){
  
  theta_star <- theta_prev
  theta_star <- runif(1, min = l * theta_prev/u , max = u * theta_prev/l)
  
  # check if variance layer theta > mean theta
  if(ls_check == TRUE)
    if(theta_star > theta_lam) 
      return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
  
  if(!is.null(YNs2)){ # woodbury log likelihood
    if(is.null(llik_prev)) {
      obj_prev <- loglw_vec(YNs2, yn, x_approx, A, exp(llam_prev), theta = theta_prev, 
                            outer = outer, v = v, calc_tau2 = calc_tau2, sep = FALSE, 
                            mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    # else{
    #   obj_prev <- loglw_vec(YNs2, yn, x_approx, A, exp(llam_prev), theta = theta_prev, 
    #                         outer = outer, v = v, calc_tau2 = calc_tau2, sep = FALSE, 
    #                         mean = mean, scale = scale, a= a, b = b)
    #   if(abs(llik_prev - obj_prev$llik) > 1e-8) stop('llik mismatch')
    #   if(abs(tau2_prev - obj_prev$tau2) > 1e-8) stop("tau2 mismatch")
    # }
    obj_prop <- loglw_vec(YNs2, yn, x_approx, A, exp(llam_prev), theta = theta_star, 
                          outer = outer, v = v, calc_tau2 = calc_tau2, sep = FALSE, 
                          mean = mean, scale = scale, a= a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  else{ # full N
    if(is.null(llik_prev)) {
      obj_prev <- logl_vec(out_vec= yn, approx = x_approx, nugs = exp(llam_prev), theta= theta_prev, 
                           outer = outer, v = v, calc_tau2 = calc_tau2, sep = FALSE, 
                           mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    obj_prop <- logl_vec(out_vec= yn, approx = x_approx, nugs = exp(llam_prev), theta= theta_star, 
                         outer = outer, v = v, calc_tau2 = calc_tau2, sep = FALSE, 
                         mean = mean, scale = scale, a= a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  
  # Priors
  prior_prev <- dgamma(theta_prev, shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(theta_star, shape = alpha, rate = beta, log = TRUE)
  
  # Accept/Reject
  l_ratio <- llik_prop + prior_prop - llik_prev - prior_prev 
  rn <- runif(1, 0, 1)
  # r <- exp(l_ratio) * (theta_prev[index]/theta_star[index])
  
  threshold <- l_ratio + log(theta_prev) - log(theta_star) 
  
  #if(runif(1, 0, 1) < exp(l_ratio) * (theta_prev[index]/theta_star[index]))
  if(log(rn) < threshold) 
    return(list(theta = theta_star, llik = llik_prop, tau2 = tau2_prop))
  else
    return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
}

mcmc_theta_l_vec <- function(x_approx, llam_prev, theta_prev, llik_prev = NULL, tau2_prev = NULL,
                             v, ls_check = FALSE, theta_y, g0, alpha, beta, l , u, inner, 
                             calc_inner_tau2, mean0, scale0 = 1, a0, b0){
  
  theta_star <- NULL
  theta_star <- runif(1, min = l * theta_prev/u , max = u * theta_prev/l)
  
  # check if variance layer theta > mean theta
  if(ls_check == TRUE)
    if(theta_star < theta_y) 
      return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
  
  # Log likelihood calculations
  if(is.null(llik_prev)) {
    obj_prev <- logl_vec(llam_prev, x_approx, g0, theta_prev, outer = inner, v = v,
                         calc_tau2 = calc_inner_tau2, sep = FALSE, mean = mean0, scale = scale0, 
                         a = a0, b = b0, latent = TRUE)
    llik_prev <- obj_prev$llik
    tau2_prev <- obj_prev$tau2
  }
  # else{
  #   obj_prev <- logl_vec(llam_prev, x_approx, g0, theta_prev, outer = inner, v = v,
  #                        calc_tau2 = calc_inner_tau2, sep = FALSE, mean = mean0, scale = scale0, 
  #                        a = a0, b = b0, latent = T)
  #   if(abs(llik_prev - obj_prev$llik) > 1e-8) stop('llik mismatch')
  #   if(abs(tau2_prev - obj_prev$tau2) > 1e-8) stop("tau2 mismatch")
  # }
  
  obj_prop <- logl_vec(llam_prev, x_approx, g0, theta_star, outer = inner, v = v,
                       calc_tau2 = calc_inner_tau2, sep = FALSE, mean = mean0, scale = scale0,
                       a = a0, b = b0, latent = TRUE)
  llik_prop <- obj_prop$llik
  tau2_prop <- obj_prop$tau2
  
  # Priors
  prior_prev <- dgamma(theta_prev, shape= alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(theta_star, shape = alpha, rate = beta, log = TRUE)
  
  # Accept/Reject
  l_ratio <- llik_prop + prior_prop - llik_prev - prior_prev 
  rn <- runif(1, 0, 1)
  # r <- exp(l_ratio) * (theta_prev[index]/theta_star[index])
  
  threshold <- l_ratio + log(theta_prev) - log(theta_star) 
  
  #if(runif(1, 0, 1) < exp(l_ratio) * (theta_prev[index]/theta_star[index]))
  if(log(rn) < threshold) 
    return(list(theta = theta_star, llik = llik_prop, tau2 = tau2_prop))
  else
    return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
}

# -------------------- Het Vdims -----------------------------------------------

ess_sample_vec_vdims <- function(YNs2 = NULL, yn, x_approx, xv = NULL, vdims, llam_prev_nv, A, llik_prev = NULL, 
                                 theta_y, theta_lam, sep = TRUE, v, outer, calc_tau2, mean, scale = 1, 
                                 mean0, scale0 = 1, g0, r0 = NULL, vec_var = TRUE, a, b) { 
  
  # Order according to Xn unique inputs
  llam_prev <- llam_new <- rep(NA, length(yn))
  llam_prev[r0$Z] <- rep(llam_prev_nv, r0$mult)
  
  # Previous llik
  if(is.null(llik_prev)){
    if(!is.null(YNs2)){ # woodbury
      obj_prev <- loglw_vec(YNs2, yn, x_approx, A= A, exp(llam_prev), theta=theta_y, outer = outer, v = v, 
                            calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
    } else{ # full N
      obj_prev <- logl_vec(YNs2, x_approx, nugs = exp(llam_prev), theta=theta_y,outer = outer, v = v, 
                           calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
    }
  }
  
  if(!vec_var) { # if no vecchia for vdims 
    xv <- r0$X0
    llam_prior_nv <- mvn_prior(xv, v, theta_lam, g0, dx = NULL, mean0, scale0 = scale0, sep = TRUE)
  }
  else{ # if vecchia for vdims
    llam_prior_nv <- rand_mvn_vec(xv, theta_lam, v=v, mean = mean0,
                               g = rep(g0, nrow(xv$x_ord)), scale = scale0, sep = sep) 
  }
  
  # Map prior draw to Xn inputs
  # llam_prior[r0$Z] <- rep(llam_prior_nv, r0$mult)
  
  # angle
  gam <- runif(1,0,2*pi)
  
  # bounds
  gam_min <- gam - 2*pi
  gam_max <- gam
  
  count=0
  while (1){
    
    count = count + 1
    if(count > 50) stop("too many iterations")
    
    llam_new_nv <- llam_prev_nv * cos (gam) + llam_prior_nv * sin (gam)
    llam_new[r0$Z] <- rep(llam_new_nv, r0$mult)
    # llam_newN <- rep(llam_new, A)
    
    # Woodbury new llik
    if(!is.null(YNs2)) {
      obj_prop <- loglw_vec(YNs2, yn, x_approx, A = A, exp(llam_new), theta=theta_y, outer = outer, 
                            v = v, calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
      llik_prop <- obj_prop$llik
      tau2 <- obj_prop$tau2
    } else{ # Full N llik
      obj_prop <- logl_vec(YNs2, x_approx, nugs = exp(llam_new), theta=theta_y, outer = outer, v = v, 
                           calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
      llik_prop <- obj_prop$llik
      tau2 <- obj_prop$tau2
    }
    
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
  
  return(list(llam = llam_new_nv, llik = llik_prop, tau2 = tau2)) # log lam new
}

# ------------------- HomGP ----------------------------------------------------

mcmc_theta_vec_hom <- function(YNs2= NULL, yn= NULL, x_approx, A=NULL, g_y, theta_prev, 
                               llik_prev = NULL, tau2_prev = NULL, v, outer, calc_tau2, mean, scale = 1, 
                               alpha, beta, l , u, a, b){
  
  theta_star <- theta_prev
  theta_star <- runif(1, min = l * theta_prev/u , max = u * theta_prev/l)
  
  # Likelihoods
  if(!is.null(YNs2)){ # Woodbury
    if(is.null(llik_prev)) {
      obj_prev <- loglw_vec(YNs2, yn, x_approx, A, g_y, theta = theta_prev, outer = outer, v = v, 
                            calc_tau2 = calc_tau2, sep = FALSE, mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    # else{
    #   obj_prev <- loglw_vec(YNs2, yn, x_approx, A, g_y, theta = theta_prev, outer = outer, v = v, 
    #                         calc_tau2 = calc_tau2, sep = FALSE, mean = mean, scale = scale, a= a, b = b)
    #   if(abs(llik_prev - obj_prev$llik) > 1e-8) stop('llik mismatch')
    #   if(abs(tau2_prev - obj_prev$tau2) > 1e-8) stop("tau2 mismatch")
    # }
      
    obj_prop <- loglw_vec(YNs2, yn, x_approx, A, g_y, theta = theta_star, outer = outer, v = v, 
                          calc_tau2 = calc_tau2, sep = FALSE, mean = mean, scale = scale, a= a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  else{ # full N
    if(is.null(llik_prev)) {
      obj_prev <- logl_vec(out_vec= yn, approx = x_approx, nugs = g_y, theta= theta_prev, outer = outer, v = v, 
                           calc_tau2 = calc_tau2, sep = FALSE,  mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
      
    obj_prop <- logl_vec(out_vec= yn, approx = x_approx, nugs = g_y, theta= theta_star, outer = outer, v = v, 
                         calc_tau2 = calc_tau2, sep = FALSE, mean = mean, scale = scale, a= a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  
  # Priors
  prior_prev <- dgamma(theta_prev, shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(theta_star, shape = alpha, rate = beta, log = TRUE)
  
  # Accept/Reject
  l_ratio <- llik_prop + prior_prop - llik_prev - prior_prev 
  rn <- runif(1, 0, 1)
  # r <- exp(l_ratio) * (theta_prev[index]/theta_star[index])
  
  threshold <- l_ratio + log(theta_prev) - log(theta_star) 
  
  #if(runif(1, 0, 1) < exp(l_ratio) * (theta_prev[index]/theta_star[index]))
  if(log(rn) < threshold) 
    return(list(theta = theta_star, llik = llik_prop, tau2 = tau2_prop))
  else
    return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
}

mcmc_theta_sep_vec_hom <- function(YNs2= NULL, yn= NULL, x_approx, A= NULL, index, g_y, theta_prev, 
                                   llik_prev = NULL, tau2_prev = NULL, v, outer, calc_tau2, mean, scale = 1, 
                                   alpha, beta, l , u, a, b){
  
  theta_star <- theta_prev
  theta_star[index] <- runif(1, min = l * theta_prev[index]/u , max = u * theta_prev[index]/l)
  
  if(!is.null(YNs2)){ # woodbury llik
    if(is.null(llik_prev)) {
      obj_prev <- loglw_vec(YNs2, yn, x_approx, A, g_y, theta = theta_prev, outer = outer, v = v, 
                             calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    # else{
    #   obj_prev <- loglw_vec(YNs2, yn, x_approx, A, g_y, theta = theta_prev, outer = outer, v = v, 
    #                         calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
    #   if(abs(llik_prev - obj_prev$llik) > 1e-8) stop('llik mismatch')
    #   if(abs(tau2_prev - obj_prev$tau2) > 1e-8) stop("tau2 mismatch")
    # }
    obj_prop <- loglw_vec(YNs2, yn, x_approx, A, g_y, theta = theta_star, outer = outer, v = v, 
                          calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  } else{ # full N
    if(is.null(llik_prev)){
      obj_prev <- logl_vec(out_vec= yn, approx = x_approx, nugs = g_y, theta= theta_prev, outer = outer, 
                            v = v, calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    } 
    obj_prop <- logl_vec(out_vec= yn, approx = x_approx, nugs = g_y, theta= theta_star, outer = outer, 
                          v = v, calc_tau2 = calc_tau2, sep = TRUE, mean = mean, scale = scale, a= a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  
  # Priors
  prior_prev <- dgamma(theta_prev[index], shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(theta_star[index], shape = alpha, rate = beta, log = TRUE)
  
  # Accept/Reject
  l_ratio <- llik_prop + prior_prop - llik_prev - prior_prev 
  rn <- runif(1, 0, 1)
  # r <- exp(l_ratio) * (theta_prev[index]/theta_star[index])
  
  threshold <- l_ratio + log(theta_prev[index]) - log(theta_star[index]) 
  
  #if(runif(1, 0, 1) < exp(l_ratio) * (theta_prev[index]/theta_star[index]))
  if(log(rn) < threshold) 
    return(list(theta = theta_star, llik = llik_prop, tau2 = tau2_prop))
  else
    return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
}

mcmc_g_vec <- function(YNs2= NULL, yn= NULL, x_approx, A=NULL, g_prev, theta_y, llik_prev = NULL, 
                       tau2_prev = NULL, v, outer, calc_tau2, mean, scale = 1, alpha, beta, l , u, sep, a, b){
  
  g_star <- runif(1, min = l * g_prev/u , max = u * g_prev/l)
  eps <- sqrt(.Machine$double.eps)
  
  # Woodbury llik
  if(!is.null(YNs2)){
    if(is.null(llik_prev)) {
      obj_prev <- loglw_vec(YNs2, yn, x_approx, A, g_prev, theta = theta_y, outer = outer, v = v, 
                            calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    # else{
    #       obj_prev <- loglw_vec(YNs2, yn, x_approx, A, g_prev, theta = theta_y, outer = outer, v = v, 
    #                             calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
    #       if(abs(llik_prev - obj_prev$llik) > 1e-8) stop('llik mismatch')
    #       if(abs(tau2_prev - obj_prev$tau2) > 1e-8) stop("tau2 mismatch")
    #     }
    obj_prop <- loglw_vec(YNs2, yn, x_approx, A, g_star, theta = theta_y, outer = outer, v = v, 
                          calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
  }
  else{ # full N
    if(is.null(llik_prev)) {
      obj_prev <- logl_vec(out_vec= yn, approx = x_approx, nugs = g_prev, theta= theta_y, outer = outer, 
                           v = v, calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    obj_prop <- logl_vec(out_vec= yn, approx = x_approx, nugs = g_star, theta= theta_y, outer = outer, 
                         v = v, calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b)
  }
  
  llik_prop <- obj_prop$llik
  
  # Priors
  prior_prev <- dgamma(g_prev - eps, shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(g_star - eps, shape = alpha, rate = beta, log = TRUE)
  
  # Accept/Reject
  l_ratio <- llik_prop + prior_prop - llik_prev - prior_prev 
  rn <- runif(1, 0, 1)
  # r <- exp(l_ratio) * (theta_prev[index]/theta_star[index])
  
  threshold <- l_ratio + log(g_prev) - log(g_star) 
  
  #if(runif(1, 0, 1) < exp(l_ratio) * (theta_prev[index]/theta_star[index]))
  if(log(rn) < threshold) {
    # print("accepted")
    return(list(g = g_star, llik = llik_prop, tau2 = obj_prop$tau2))
  }
  else
    return(list(g = g_prev, llik = llik_prev, tau2 = tau2_prev))
}

mcmc_g_vec_inner <- function(llam, x_approx, g_prev, theta, llik_prev = NULL, tau2_prev = NULL, 
                       v, outer, calc_tau2, mean, scale = 1, alpha, beta, l , u, sep, a, b){
  
  eps <- sqrt(.Machine$double.eps)
  g_star <- runif(1, min = l * g_prev/u , max = u * g_prev/l)
  
  if(is.null(llik_prev)) {
      obj_prev <- logl_vec(out_vec= llam, approx = x_approx, nugs = g_prev, theta= theta, outer = outer, v = v, 
                           calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b, latent = TRUE)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
  }
  # else{
  #   obj_prev <- logl_vec(out_vec= llam, approx = x_approx, nugs = g_prev, theta= theta, outer = outer, v = v, 
  #                        calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b, latent = T)
  #   if(abs(llik_prev - obj_prev$llik) > 1e-8) stop('llik mismatch')
  #   if(abs(tau2_prev - obj_prev$tau2) > 1e-8) stop("tau2 mismatch")
  # }
  obj_prop <- logl_vec(out_vec= llam, approx = x_approx, nugs = g_star, theta= theta, outer = outer, 
                       v = v, calc_tau2 = calc_tau2, sep = sep, mean = mean, scale = scale, a= a, b = b, latent = TRUE)
  
  llik_prop <- obj_prop$llik
  
  # Priors
  prior_prev <- dgamma(g_prev - eps, shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(g_star - eps, shape = alpha, rate = beta, log = TRUE)
  
  # Accept/Reject
  l_ratio <- llik_prop + prior_prop - llik_prev - prior_prev 
  rn <- runif(1, 0, 1)
  # r <- exp(l_ratio) * (theta_prev[index]/theta_star[index])
  
  threshold <- l_ratio + log(g_prev) - log(g_star) 
  
  #if(runif(1, 0, 1) < exp(l_ratio) * (theta_prev[index]/theta_star[index]))
  if(log(rn) < threshold) 
    return(list(g = g_star, llik = llik_prop, tau2 = obj_prop$tau2))
  else
    return(list(g = g_prev, llik = llik_prev, tau2 = tau2_prev))
}
