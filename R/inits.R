# maybe use lam directly for small problems and llam for bigger? check.
# make 3 functions

#### Theta note###

# theta(mk) = sqrt(theta(p)/2 * nu)
# theta(p) = 2 * nu * (theta(mk)^2)

# inits_pred: predictions from hetGP
inits_het <- function(reps, v){

  out <- list()
  x <- reps$X0
  y <- reps$Z0
  ylist <- unlist(reps$Zlist)
  
  X = list(X0 = x, Z0 = y, mult = reps$mult)
  
  if(v == 999) {covtype = "Gaussian"; v = 4.5}
  else if(v == 1000 + 1.5 || v == 1.5) {covtype = "Matern3_2"; v =1.5}
  else if(v == 1000 + 2.5 || v == 2.5) {covtype = "Matern5_2"; v =2.5}
  
  fit <- hetGP::mleHetGP(X, ylist, covtype = covtype)#, 
                         # settings = list(linkThetas = "none", checkHom = FALSE))
  
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
              g = fit$g, v = v)
  return(out)
}

inits_sV <- function(reps, v){
  
  # dim
  out <- list()
  x <- reps$X0
  y <- reps$Z0

  if(v == 999) v1 <- 4.5
  if(v == 1001.5 || v == 1002.5) v1 <- v - 1000
  if(v == 1.5 || v == 2.5) v1 <- v
  
  sVfit <- fit_scaled(y, x, nu = v1, ms = min(30, length(y) - 1))
  sVpred <- predictions_scaled(sVfit, x, m = min(100, length(y) - 1))
  
  out <- list(reps = reps, ty = 2 * v1 * (sVfit$covparms[2:(ncol(x) +1)])^2,
              tg = 2 * v1 * (sVfit$covparms[2:(ncol(x) +1)])^2, 
              g = sVfit$covparms[ncol(x) +2],
              mu_y = sVpred, scale = sVfit$covparms[1], v = v1)
  
  return(out)
}

inits_bhetgp <- function(reps, v, vec){
  
  n <- length(reps$Z0)
  if(!vec)
    out <- inits_het(reps, v = v)
  else 
    out <- inits_sV(reps, v = v) 

  lam <- unlist(sapply(1:length(reps$Zlist), function(i) mean((reps$Zlist[[i]] - out$mu_y[i])^2/out$scale)))
  out$llam <- log(lam)
  
  sm <- smoother(out)
  if(vec){
    out$tg <- sm$tg
    out$scale_lam <- sm$scale_lam
  }
  out$llam_w <- sm$llam_s
  out$mean0 <- mean(out$llam_w)
  
  out$tg[out$tg > 5] <- 5
  out$tg[out$tg < 0.001] <- 0.001
  
  out$ty[out$ty > 5] <- 5
  out$ty[out$ty < 1e-3] <- 1e-3
  
  out$llam_w[out$llam_w < -10] <- -10
  
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
  
  if(length(init$llam) > ncol(init$reps$X0) * 100){
    ih <- sample(1:length(init$llam), ncol(init$reps$X0) * 100) 
  }else{
    ih <- 1:length(init$llam)
  }
  
  hgp <- mleHomGP(init$reps$X0[ih, ], init$llam[ih])
  p <- predict(hgp, init$reps$X0)
  init$llam_s <- p$mean
  init$tg <- hgp$theta
  init$scale_lam <- hgp$nu_hat
  
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
  if(!vec)
    out <- inits_hom(reps, v = v)
  else
    out <- inits_sV(reps, v = v)
  return(out)
}  

# Still working on the init
inits_vdims <- function(reps, reps_vdims, v, vec, vdims){
  out <- list()
  if(!vec){
    out <- inits_het(reps, v)
  }else{
    out <- inits_sV(reps, v)
  }
  
  muy <- out$mu_y
  lam <- unlist(sapply(1:length(reps$Zlist), 
                        function(i) mean((reps$Zlist[[i]] - muy[i])^2/out$scale)))
  
  # We want log of means of lam ---> gives better scores.
  r_vdims <- find_reps(reps$X0[, vdims], 1:nrow(reps$X0))
  lam <- unlist(sapply(1:length(r_vdims$Zlist), function(i) mean(lam[r_vdims$Zlist[[i]]])))

  out$llam <- log(lam)

  # r_vdims <- find_reps(reps$X0[, vdims], 1:nrow(reps$X0))
  # out$llam <- unlist(sapply(1:length(r_vdims$Zlist), function(i) mean(out$llam[r_vdims$Zlist[[i]]])))

  # add smoother function
  if(length(out$llam) > ncol(r_vdims$X0) * 100){
    ih <- sample(1:length(out$llam), ncol(r_vdims$X0) * 100) 
  }else{
    ih <- 1:length(out$llam)
  }
  
  hgp <- mleHomGP(r_vdims$X0[ih, ], out$llam[ih])
  p <- predict(hgp, r_vdims$X0)
  out$llam_w <- p$mean
  out$tg <- hgp$theta
  out$scale_lam <- hgp$nu_hat
  
  out$mean0 <- mean(out$llam_w)
  
  out$tg[out$tg > 5] <- 5
  out$tg[out$tg < 0.001] <- 0.001
  
  out$ty[out$ty > 5] <- 5
  out$ty[out$ty < 1e-3] <- 1e-3
  
  out$llam_w[out$llam_w < -10] <- -10
  
  return(out)
}

# 
# inits_vdims <- function(reps, reps_vdims, v, vec, vdims){
#   out <- list()
#   if(!vec){
#     out <- inits_het(reps, v)
#   }else{
#     out <- inits_sV(reps, v)
#   }
#   
#   muy <- out$mu_y
#   ymap <- find_reps(reps$X0[, vdims, drop = FALSE], reps$Z0)
#   
#   lam <- unlist(sapply(1:length(ymap$Zlist), 
#                        function(i) mean((ymap$Zlist[[i]] - muy[i])^2/out$scale)))
#   out$llam <- log(lam)
#   
#   # add smoother function
#   if(length(out$llam) > ncol(reps_vdims$X0) * 100){
#     ih <- sample(1:length(out$llam), ncol(reps_vdims$X0) * 100) 
#   }else{
#     ih <- 1:length(out$llam)
#   }
#   
#   hgp <- mleHomGP(reps_vdims$X0[ih, ], out$llam[ih])
#   p <- predict(hgp, reps_vdims$X0)
#   out$llam_w <- p$mean
#   out$tg <- hgp$theta
#   out$scale_lam <- hgp$nu_hat
#   
#   out$mean0 <- mean(out$llam_w)
#   
#   out$tg[out$tg > 5] <- 5
#   out$tg[out$tg < 0.001] <- 0.001
#   
#   out$ty[out$ty > 5] <- 5
#   out$ty[out$ty < 1e-3] <- 1e-3
#   
#   out$llam_w[out$llam_w < -10] <- -10
#   
#   # sV <- fit_scaled(out$llam, reps_vdims$X0, ms = min(30, length(reps_vdims$Z0) - 1), nu = out$v)
#   # ls <- predictions_scaled(sV, reps_vdims$X0, m = min(100, length(reps_vdims$Z0) - 1))
#   # 
#   # out$tg <- 2 * out$v * sV$covparms[2:(ncol(reps_vdims$X0) +1)]^2
#   # out$scale0 <- sV$covparms[1]
#   # out$mean0 <- sV$betahat
#   # out$llam_w <- ls
#   
#   return(out)
# }

