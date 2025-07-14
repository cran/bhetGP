# -------------------------  MCMC HomGP + HetGP --------------------------------

# ------- HetGP ----------------------------------------------------------------
# ess_sample: ESS sampling for Lambda layer
# mcmc_theta_y_sep: MCMC sampling for theta y (sep case)
# mcmc_theta_l_sep: MCMC sampling for theta lam (sep case)
# mcmc_theta_y: MCMC sampling for theta y (isotropic)
# mcmc_theta_l: MCMC sampling for theta l (isotropic)
# ess_sample_vdims: ESS for Lambda layer when lambdas change in some dims

# ------- HomGP ----------------------------------------------------------------
# mcmc_g: MCMC sampling for g isotropic
# mcmc_g_sep: MCMC sampling for g seperable (use logl)
# mcmc_theta_sep_hom: MCMC sampling for theta y sep 
# mcmc_theta_hom: MCMC sampling for theta y isotropic
# ------------------------------------------------------------------------------

# ------------ HetGP functions -------------------------------------------------

ess_sample <- function(YNs2=NULL, yn, xn, llam_prev, A = NULL, llik_prev = NULL, 
                       theta_y, theta_lam, sep = TRUE, v, outer, calc_tau2, mean, scale = 1, a, b, 
                       mean0, scale0 = 1, g0, dx_n = NULL) { 
  
  # Calculate previous llik
  if(is.null(llik_prev)){
    if(!is.null(YNs2)){ # woodbury version
      if(sep) llik_prev <- loglw(ys2=YNs2, yn=yn, xn=xn, exp(llam_prev), 
                                 A=A, v = v, theta = theta_y, outer = outer, 
                                 calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)$llik
      else llik_prev <- loglw_iso(ys2=YNs2, yn=yn, dx_n, exp(llam_prev), 
                                  A=A, v = v, theta = theta_y, outer = outer, 
                                  calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)$llik
    } else { # full N
      if(sep) llik_prev <- logl(yn, xn, exp(llam_prev), theta_y, v = v, outer = outer, 
                                calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)$llik
      else llik_prev <- logl_iso(yn, dx_n, exp(llam_prev), theta_y, v = v, outer = outer, 
                                 calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)$llik
    } 
  }

  # draw a prior
  llam_prior <- mvn_prior(x = xn, v, theta_lam, g0, dx = dx_n, mean0, scale0, sep = sep)

  # angle
  gam <- runif(1, 0 , 2*pi)
  
  # bounds
  gam_min <- gam - 2*pi
  gam_max <- gam
  
  while (1){
    
    llam_new <- llam_prev * cos (gam) + llam_prior * sin (gam)

    if(!is.null(YNs2)) {
      if(sep) obj_prop <- loglw(ys2=YNs2, yn=yn, xn=xn, exp(llam_new), A= A, v = v, theta = theta_y, 
                                outer = outer, calc_tau2 = calc_tau2,  mean = mean, scale = scale, 
                                a = a, b = b)
      else obj_prop <- loglw_iso(ys2=YNs2, yn=yn, dx_n, exp(llam_new), A= A, v = v, theta = theta_y, 
                                 outer = outer, calc_tau2 = calc_tau2,  mean = mean, scale = scale, 
                                 a = a, b = b)
      llik_prop <- obj_prop$llik
      tau2 <- obj_prop$tau2
    }
    else {
      
      if(sep)
        obj_prop <- logl(yn, xn, nugs = exp(llam_new), theta_y, v = v, outer = outer, 
                         calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      else 
        obj_prop <- logl_iso(yn, dx_n, nugs = exp(llam_new), theta_y, v = v, outer = outer, 
                             calc_tau2 = calc_tau2,  mean = mean, scale = scale, a = a, b = b)
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
  
  return(list(llam = llam_new, llik = llik_prop, tau2 = tau2))
}

mcmc_theta_y_sep <- function(YNs2= NULL, yn= NULL, xn, A=NULL, index, llam_prev, theta_prev, llik_prev = NULL, 
                             tau2_prev = NULL, v, ls_check = FALSE, theta_lam = NULL, outer, 
                             calc_tau2, mean = mean, scale = scale, alpha, beta, l , u, a, b){
  
  # return(list(theta = theta_prev, llik = NULL, tau2 = NULL))
  theta_star <- theta_prev
  theta_star[index] <- runif(1, min = l * theta_prev[index]/u , max = u * theta_prev[index]/l)
  
  # If variance theta > mean theta, then check
  if(ls_check == TRUE & !is.null(theta_lam)){
      if(theta_star[index] > theta_lam) 
        return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
  } 
  
  # Likelihood previous and current
  if(!is.null(YNs2)){ # woodbury
    if(is.null(llik_prev)) {
      obj_prev <- loglw(YNs2, yn, xn, exp(llam_prev), A = A, v = v, theta= theta_prev, outer = outer, 
                        calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }

    obj_prop <- loglw(YNs2, yn, xn, exp(llam_prev), A = A, v = v, theta= theta_star, outer = outer, 
                      calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  else{ # full N
    if(is.null(llik_prev)) {
      obj_prev <- logl(yn, x= xn, nugs = exp(llam_prev), theta= theta_prev, v = v, outer = outer, 
                        calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
      
    obj_prop <- logl(yn, x= xn, nugs = exp(llam_prev), theta= theta_star, v = v, outer = outer, 
                     calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  
  # Priors
  prior_prev <- dgamma(theta_prev[index], shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(theta_star[index], shape = alpha, rate = beta, log = TRUE)
  
  # Accept reject
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

mcmc_theta_l_sep <- function(xlam, llam_prev, index, theta_prev, llik_prev = NULL, tau2_prev = NULL,
                             v, ls_check = FALSE, theta_y, g0, alpha, beta, l , u, inner, 
                             calc_inner_tau2, mean0, scale0, a0, b0){
  
  theta_star <- theta_prev
  theta_star[index] <- runif(1, min = l * theta_prev[index]/u , max = u * theta_prev[index]/l)
  
  # If variance theta > mean theta, tehn check
  if(ls_check == TRUE)
    if(theta_star[index] < theta_y) 
      return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
  
  if(is.null(llik_prev)) {
    obj_prev <- logl(llam_prev, xlam, g0, theta_prev, v = v, outer = inner, calc_tau2 = calc_inner_tau2, 
                      mean = mean0, scale = scale0, a = a0, b = b0)
    llik_prev <- obj_prev$llik
    tau2_prev <- obj_prev$tau2
  }
  # if(llik_prev != logl(llam_prev, xlam, g0, theta_prev, v = v, outer = inner, tau2 = inner_tau2, 
  #                      mean = mean0, scale = scale0, a = a0, b = b0)$llik) stop ("mismatch")
  
  obj_prop <- logl(llam_prev, xlam, g0, theta_star, v = v, outer = inner, calc_tau2 = calc_inner_tau2, 
                    mean = mean0, scale = scale0, a = a0, b = b0)
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

mcmc_theta_y <- function(YNs2= NULL, yn= NULL, dx_n = NULL, A=NULL, llam_prev, theta_prev, llik_prev = NULL, 
                         tau2_prev = NULL, v, ls_check = FALSE, theta_lam, outer, calc_tau2, 
                         mean, scale = 1, alpha, beta, l , u, a, b){
  
  theta_star <- theta_prev
  theta_star <- runif(1, min = l * theta_prev/u , max = u * theta_prev/l)
  
  # If variance theta > mean theta, tehn check
  if(ls_check == TRUE)
    if(theta_star > theta_lam) 
      return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
  
  # Likelihood calc
  if(!is.null(YNs2)){
    if(is.null(llik_prev)) {
      obj_prev <- loglw_iso(YNs2, yn, dx_n, exp(llam_prev), A= A, v = v, theta= theta_prev, outer = outer, 
                            calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    obj_prop <- loglw_iso(YNs2, yn, dx_n, exp(llam_prev), A= A, v = v, theta= theta_star, outer = outer,
                          calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  else{
    if(is.null(llik_prev)) {
      obj_prev <- logl_iso(yn, dx_n, nugs = exp(llam_prev), theta= theta_prev, v = v, outer = outer, 
                           calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    obj_prop <- logl_iso(yn, dx_n, nugs = exp(llam_prev), theta= theta_star, v = v, outer = outer, 
                         calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
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

mcmc_theta_l <- function(dx_lam, llam_prev, theta_prev, llik_prev = NULL, tau2_prev = NULL, v, 
                         ls_check = FALSE, theta_y, g0, alpha, beta, l , u, inner, calc_inner_tau2, 
                         mean0, scale0, a0, b0){
  
  theta_star <- theta_prev
  theta_star <- runif(1, min = l * theta_prev/u , max = u * theta_prev/l)
  
  # If variance theta > mean theta, then check
  if(ls_check == TRUE)
    if(theta_star < theta_y) 
      return(list(theta = theta_prev, llik = llik_prev, tau2 = tau2_prev))
  
  # Likelihood
  if(is.null(llik_prev)) {
    obj_prev <- logl_iso(llam_prev, dx_lam, g0, theta_prev, v = v, outer = inner, 
                         calc_tau2 = calc_inner_tau2, mean = mean0, scale = scale0, a = a0, b = b0)
    llik_prev <- obj_prev$llik
    tau2_prev <- obj_prev$tau2
  }
  
  obj_prop <- logl_iso(llam_prev, dx_lam, g0, theta_star, v = v, outer = inner, 
                       calc_tau2 = calc_inner_tau2, mean = mean0, scale = scale0, a = a0, b = b0)
  llik_prop <- obj_prop$llik
  tau2_prop <- obj_prop$tau2
  
  # Prios
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

# ------------ Vdims function --------------------------------------------------

ess_sample_vdims <- function(YNs2=NULL, yn, xn, vdims, llam_prev_nv, A, llik_prev, tau2_prev = NULL,
                             theta_y, theta_lam, sep = TRUE, v, outer, calc_tau2, mean, scale, mean0, scale0, g0, a, b,
                             r0 = NULL, dx_n = NULL, dx_lam =NULL) { 
  
  # map nv lambdas to Xn locations
  llam_prev <- llam_new <- rep(NA, length(yn))
  llam_prev[r0$Z] <- rep(llam_prev_nv, r0$mult)
  
  # Calculate previous llik
  if(is.null(llik_prev)){
    if(!is.null(YNs2)){
      if(sep) llik_prev <- loglw(ys2=YNs2, yn=yn, xn=xn, exp(llam_prev), A=A, v = v, theta = theta_y, outer = outer, 
                                 calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)$llik
      else llik_prev <- loglw_iso(ys2=YNs2, yn=yn, dx_n, exp(llam_prev), A=A, v = v, theta = theta_y, outer = outer, 
                                  calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)$llik
    } 
    else {
      if(sep) llik_prev <- logl(yn, xn, nugs = exp(llam_prev), theta_y, v = v, outer = outer, 
                                calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)$llik
      else llik_prev <- logl_iso(yn, dx_n, nugs = exp(llam_prev), theta_y, v = v, outer = outer, 
                                 calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)$llik
    }
  }
  warning("check if llik match")
  
  # Extract unique nv X's
  xv <- r0$X0
  
  # Draw nv unique lambdas
  llam_prior_nv <- mvn_prior(xv, v, theta_lam, g0, dx = NULL, mean0, scale0, sep = sep)
  
  # angle
  gam <- runif(1,0,2*pi)
  
  # bounds
  gam_min <- gam - 2*pi
  gam_max <- gam
  
  while (1){
    
    # Calculate new lambdas and map to Xn
    llam_new_nv <- llam_prev_nv * cos(gam) + llam_prior_nv * sin(gam)
    llam_new[r0$Z] <- rep(llam_new_nv, r0$mult)

    if(!is.null(YNs2)) {
      if(sep) obj_prop <- loglw(ys2=YNs2, yn=yn, xn=xn, exp(llam_new), A= A, v = v, theta = theta_y, outer = outer, 
                                calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b) 
      else obj_prop <- loglw_iso(ys2=YNs2, yn=yn, dx_n, exp(llam_new), A= A, v = v, theta = theta_y, outer = outer, 
                                 calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prop <- obj_prop$llik
      tau2 <- obj_prop$tau2
    }
    
    else {
      if(sep)
        obj_prop <- logl(yn, xn, nugs = exp(llam_new), theta_y, v = v, outer = outer, calc_tau2 = calc_tau2, 
                         mean = mean, scale = scale, a = a, b = b)
      else 
        obj_prop <- logl_iso(yn, dx_n, nugs = exp(llam_new), theta_y, v = v, outer = outer, calc_tau2 = calc_tau2, 
                             mean = mean, scale = scale, a = a, b = b)
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
      gam <- runif(1,gam_min,gam_max)
    }
    else
      break
  }
  
  return(list(llam = llam_new_nv, llik = llik_prop, tau2 = tau2)) # log lam new
}

# ------------ HomGP functions -------------------------------------------------

mcmc_g <- function(YNs2 = NULL, yn, dx_n, A = NULL, g_prev, theta_y, llik_prev = NULL,
                   tau2_prev = NULL, v, outer, calc_tau2,  mean, scale, alpha, beta, l , u, a, b){
  
  eps <- sqrt(.Machine$double.eps)
  g_star <- runif(1, min = l * g_prev/u , max = u * g_prev/l)
  
  # If woodbury
  if(!is.null(YNs2)){
    if(is.null(llik_prev)) {
      obj_prev <- loglw_iso(YNs2, yn, dx_n, g_prev, A= A, v = v, theta= theta_y, outer = outer, 
                            calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    obj_prop <- loglw_iso(YNs2, yn, dx_n, g_star, A= A, v = v, theta= theta_y, outer = outer,
                           calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
  }
  else{
    if(is.null(llik_prev)) {
      obj_prev <- logl_iso(yn, dx_n, nugs = g_prev, theta= theta_y, v = v, outer = outer, 
                           calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    obj_prop <- logl_iso(yn, dx_n, nugs = g_star, theta= theta_y, v = v, outer = outer, 
                         calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
  }
  
  llik_prop <- obj_prop$llik
  
  # priors
  prior_prev <- dgamma(g_prev - eps, shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(g_star - eps, shape = alpha, rate = beta, log = TRUE)
  
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


mcmc_g_sep <- function(YNs2 = NULL, yn, xn, A = NULL, g_prev, theta_y, llik_prev = NULL, tau2_prev = NULL,
                       v, outer, calc_tau2, mean, scale, alpha, beta, l , u, a, b){
  
  eps <- sqrt(.Machine$double.eps)
  g_star <- runif(1, min = l * g_prev/u , max = u * g_prev/l)
  
  # Woodbury version
  if(!is.null(YNs2)){
    if(is.null(llik_prev)) {
      obj_prev <- loglw(YNs2, yn, xn, g_prev, A = A, v = v, theta= theta_y, outer = outer, 
                        calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }

    obj_prop <- loglw(YNs2, yn, xn, g_star, A = A, v = v, theta= theta_y, outer = outer, 
                      calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
  }
  else{
    if(is.null(llik_prev)) {
      obj_prev <- logl(yn, x= xn, nugs = g_prev, theta= theta_y, v = v, outer = outer, 
                       calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    obj_prop <- logl(yn, x= xn, nugs = g_star, theta= theta_y, v = v, outer = outer, 
                     calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
  }
  
  llik_prop <- obj_prop$llik
  
  # Prios
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

mcmc_g_inner <- function(yn, dx_n, g_prev, theta, llik_prev = NULL, tau2_prev = NULL, v, 
                         outer, calc_tau2, mean, scale, alpha, beta, l , u, a, b){
  
  eps <- sqrt(.Machine$double.eps)
  g_star <- runif(1, min = l * g_prev/u , max = u * g_prev/l)
  
  if(is.null(llik_prev)){
    obj_prev <- logl_iso(yn, dx_n, nugs = g_prev, theta= theta, v = v, outer = outer, 
                         calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
    llik_prev <- obj_prev$llik
    tau2_prev <- obj_prev$tau2
  }
  
  obj_prop <- logl_iso(yn, dx_n, nugs = g_star, theta= theta, v = v, outer = outer, 
                       calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
  
  llik_prop <- obj_prop$llik
  
  # Prios
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

mcmc_g_sep_inner <- function(yn, xn, g_prev, theta, llik_prev = NULL, tau2_prev = NULL, v, 
                       outer, calc_tau2, mean, scale, alpha, beta, l , u, a, b){
  
  g_star <- runif(1, min = l * g_prev/u , max = u * g_prev/l)
  tau2_prev <- NULL
  
  if(is.null(llik_prev)){
      obj_prev <- logl(yn, x= xn, nugs = g_prev, theta= theta, v = v, outer = outer, 
                       calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
  }
  
  obj_prop <- logl(yn, x= xn, nugs = g_star, theta= theta, v = v, outer = outer, 
                   calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
  
  llik_prop <- obj_prop$llik
  
  # Prios
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

mcmc_theta_sep_hom <- function(YNs2= NULL, yn= NULL, xn, A=NULL, index, g_y, theta_prev, llik_prev = NULL,
                               tau2_prev = NULL, v, outer, calc_tau2, mean, scale, alpha, beta, l , u, a, b){
  
  theta_star <- theta_prev
  theta_star[index] <- runif(1, min = l * theta_prev[index]/u , max = u * theta_prev[index]/l)
  
  # If woodbury
  if(!is.null(YNs2)){
    if(is.null(llik_prev)) {
      obj_prev <- loglw(YNs2, yn, xn, g_y, A = A, v = v, theta= theta_prev, outer = outer, 
                        calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    # if(abs(llik_prev - loglw(YNs2, yn, xn, g_y, A = A, v = v, theta= theta_prev, outer = outer, 
    #                     calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)$llik) > 1e-8) 
    #   stop('mismatch')
    obj_prop <- loglw(YNs2, yn, xn, g_y, A = A, v = v, theta= theta_star, outer = outer, 
                      calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  else{
    if(is.null(llik_prev)) {
      obj_prev <- logl(yn, x= xn, nugs = g_y, theta= theta_prev, v = v, outer = outer, 
                       calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    }
    obj_prop <- logl(yn, x= xn, nugs = g_y, theta= theta_star, v = v, outer = outer, 
                     calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  
  # Priors
  prior_prev <- dgamma(theta_prev[index], shape = alpha, rate = beta, log = TRUE)
  prior_prop <- dgamma(theta_star[index], shape = alpha, rate = beta, log = TRUE)
  
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

mcmc_theta_hom <- function(YNs2= NULL, yn= NULL, dx_n = NULL, A=NULL, g_y, theta_prev, llik_prev = NULL,
                           tau2_prev = NULL, v, ls_check = FALSE, theta_lam, outer, calc_tau2, mean, scale, 
                           alpha, beta, l , u, a, b){
  
  theta_star <- theta_prev
  theta_star <- runif(1, min = l * theta_prev/u , max = u * theta_prev/l)

  if(!is.null(YNs2)){
    if(is.null(llik_prev)){
      obj_prev <- loglw_iso(YNs2, yn, dx_n, g_y, A= A, v = v, theta= theta_prev, outer = outer, 
                            calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    } 
    obj_prop <- loglw_iso(YNs2, yn, dx_n, g_y, A= A, v = v, theta= theta_star, outer = outer, 
                          calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
    llik_prop <- obj_prop$llik
    tau2_prop <- obj_prop$tau2
  }
  
  else{
    if(is.null(llik_prev)){
      obj_prev <- logl_iso(yn, dx_n, nugs = g_y, theta= theta_prev, v = v, outer = outer, 
                           calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
      llik_prev <- obj_prev$llik
      tau2_prev <- obj_prev$tau2
    } 
    obj_prop <- logl_iso(yn, dx_n, nugs = g_y, theta= theta_star, v = v, outer = outer, 
                         calc_tau2 = calc_tau2, mean = mean, scale = scale, a = a, b = b)
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
