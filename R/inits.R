# maybe use lam directly for small problems and llam for bigger? check.
# make 3 functions

# inits_pred: predictions from hetGP
inits_het <- function(reps, v){

  out <- list()
  x <- reps$X0
  y <- reps$Z0
  ylist <- unlist(reps$Zlist)
  
  X = list(X0 = x, Z0 = y, mult = reps$mult)
  
  if(v == 999) covtype = "Gaussian"
  else if(v == 1000 + 1.5 || v == 1.5) covtype = "Matern3_2"
  else if(v == 1000 + 2.5 || v == 2.5) covtype = "Matern5_2"
  
  fit <- hetGP::mleHetGP(X, ylist, covtype = covtype, settings = list(linkThetas = "none",
                                                               checkHom = FALSE))
  
  # Update how we pick the values
  pred <- predict(fit, reps$X0)
  if(covtype != "Gaussian") {
    tg <- fit$theta_g^2
    ty <- fit$theta^2
  }else{
    tg <- fit$theta_g
    ty <- fit$theta
  }
  
  out <- list(reps = reps, ty = ty, tg = tg, mean0 = fit$nmean, scale0 = fit$nu_hat_var, 
              llam_pred = log(pred$nugs/fit$nu_hat), mu_y = pred$mean, scale = fit$nu_hat,
              g = fit$g)
  return(out)
}

inits_sV <- function(reps, v){

  out <- list()
  x <- reps$X0
  y <- reps$Z0
  
  if(v == 999) v1 <- 4.5
  if(v == 1001.5 || v == 1002.5) v1 <- v - 1000
  
  sVfit<- fit_scaled(y, x, nu = v1, ms = min(30, length(y) - 1))
  sVpred <- predictions_scaled(sVfit, x, m = min(100, length(y) - 1))
  
  out <- list(reps = reps, ty = sVfit$covparms[2:(ncol(x) +1)]^2,
              tg = sVfit$covparms[2:(ncol(x) +1)]^2, g = sVfit$covparms[ncol(x) +2],
              mu_y = sVpred, scale = sVfit$covparms[1])
  
  return(out)
}

inits_bhetgp <- function(reps, v, vec){
  
  n <- length(reps$Z0)
  
  if(!vec)# || ncol(reps$X0) == 1)
    out <- inits_het(reps, v = v)
  else 
    out <- inits_sV(reps, v = v) 

  lam <- unlist(sapply(1:length(reps$Zlist), function(i) mean((reps$Zlist[[i]] - out$mu_y[i])^2/out$scale)))
  
  # lamN <- (reps$Z - rep(out$mu_y, times = reps$mult))^2 # This gives N est.
  # lamN <-  lamN / out$scale  # scale correctly
  # lam_init <- drop(hetGP:::fast_tUY2(reps$mult, lamN))/(reps$mult) # average to n est. 
  out$llam <- log(lam)
  
  sm <- smoother(out)
  if(vec){
    out$tg <- sm$tg
    out$scale_lam <- sm$scale_lam
  }
  out$llam_w <- sm$llam_s
  out$mean0 <- mean(out$llam_w)
  
  return(out)
}  

inits_emp <- function(reps){
  Sy2 <- unlist(lapply(reps$Zlist, var))
  return(Sy2)
}

inits_flat <- function(reps){
  n <- nrow(reps$X0)
  d <- ncol(reps$X0)
  
  out <- list(reps = reps, ty = rep(1, d), tg = rep(2, d), mean0 = 0, scale0 = 1, 
              llam_w = rep(log(var(reps$Z) * 0.1), n), mu_y = 0, scale = 1, g = 1e-5)
}

# find a GP
smoother <- function(init){
  
  if(ncol(init$reps$X0) == 1)
    ls <- stats::loess(init$llam~init$reps$X0)$fitted
  else{
    sV <- fit_scaled(init$llam, init$reps$X0, ms = min(30, length(init$reps$Z0) - 1))
    ls <- predictions_scaled(sV, init$reps$X0, m = min(100, length(init$reps$Z0) - 1))
    init$tg <- sV$covparms[2:(ncol(init$reps$X0) +1)]^2
    init$scale_lam <- sV$covparms[1]
  }
  init$llam_s <- ls
  
  return(init)
}

inits_hom <- function(reps, v){
  
  out <- list()
  x <- reps$X0
  y <- reps$Z0
  ylist <- unlist(reps$Zlist)
  
  X = list(X0 = x, Z0 = y, mult = reps$mult)
  
  if(v == 999) covtype = "Gaussian"
  else if(v == 1000 + 1.5 || v == 1.5) covtype = "Matern3_2"
  else if(v == 1000 + 2.5 || v == 2.5) covtype = "Matern5_2"
  
  fit <- hetGP::mleHomGP(X, ylist, covtype = covtype)
  pred <- predict(fit, reps$X0)
  
  if(covtype != "Gaussian") {
    ty <- fit$theta^2
  }else{
    ty <- fit$theta
  }
  
  out <- list(reps = reps, ty = ty, mu_y= pred$mean, scale = fit$nu_hat, g = fit$g)
  
  return(out)
}

inits_bhomgp <- function(reps, v, vec){
  
  if(!vec)# || ncol(reps$X0) == 1)
    out <- inits_hom(reps, v = v)
  else
    out <- inits_sV(reps, v = v)
  return(out)
}  

# Still working on the init
inits_vdims <- function(reps, reps_vdims, v, vec){
  
  out <- list()
  
  if(!vec){
    out <- inits_het(reps, v)
    fitv <- inits_het(reps_vdims, v)
  
    out$tg <- fitv$tg
    out$mean0 <- fitv$mean0
    out$scale0 <- fitv$scale0
    out$llam_pred <- fitv$llam_pred
  }else{
    out <- inits_sV(reps, v)
    # if(ncol(reps_vdims$X0) == 1){
    #   fitv <- inits_het(reps_vdims, v)
    # }
    # else{
      fitv <- inits_sV(reps_vdims, v)
    # }
  }
  
  lam <- unlist(sapply(1:length(reps_vdims$Zlist), 
                       function(i) mean((reps_vdims$Zlist[[i]] - fitv$mu_y[i])^2/out$scale)))
  
  # lam <- (reps_vdims$Z - rep(fitv$mu_y, times = reps_vdims$mult))^2 # This gives N est.
  # lam <-  lam / out$scale  # scale correctly
  # lam_init <- drop(hetGP:::fast_tUY2(reps_vdims$mult, lam))/(reps_vdims$mult) # average to n est. 
  out$llam <- log(lam)
    
  # if(ncol(reps_vdims$X0) == 1)
  #   ls <- loess(out$llam~reps_vdims$X0)$fitted
  # else{
    sV <- fit_scaled(out$llam, reps_vdims$X0, ms = min(30, length(reps_vdims$Z0) - 1))
    ls <- predictions_scaled(sV, reps_vdims$X0, m = min(100, length(reps_vdims$Z0) - 1))
    out$tg <- sV$covparms[2:(ncol(reps_vdims$X0) +1)]^2
    out$scale0 <- sV$covparms[1]
    out$mean0 <- sV$betahat
  # }
  
  out$llam_w <- ls
  
  return(out)
}

