# ------------ Functions ------------------------------------------------------
# predict_vec: non vecchia predictions for hetgp w/ vecchia
# predict_vec_vdims: non vecchia preds for hetgp w/ vdims w/ vecchia
# predict_vec_hom: homGP preds w/ vecchia
# gp_vec_lam: calculates pred mean and s2 for lambda layer w/ vecchia
# gp_vec_y: calculates pred mean and s2 for y layer w/ vecchia
# gp_vec_hom: calclulates pred mean and s2 for homgp w/ vecchia
# -----------------------------------------------------------------------------

predict_vec <- function(object, x_new, m = object$m, settings) {

  tic <- proc.time()[[3]]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  object$m_pred <- m

  if (settings$return_all & !settings$lite)
    stop("return_all only offered when lite = TRUE")

  if (!is.null(settings$ordering_new)) {
    if (settings$lite) message("ordering_new is only relevant when lite = FALSE")
    test <- check_ordering(settings$ordering_new, nrow(x_new))
  }

  # Pre-calculate nearest neighbors for x
  if (settings$lite) {
    NN_x_new <- FNN::get.knnx(object$x, x_new, m)$nn.index
    x_approx <-  object$x_approx
  } else {
    NN_x_new <- NULL
    x_approx <- add_pred_to_approx(object$x_approx, x_new, m, settings$ordering_new)
  }

  if (settings$cores == 1) { # run serial for loop

    mean_y <- mean_lam <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (settings$lite) {
      if(settings$interval == "pi" ||settings$interval == "both")
        s2_sum_y <- s2_sum_lam <- rep(0, times = n_new)
      if(settings$interval == "ci" ||settings$interval == "both")
        s2_sum_y_ci <- s2_sum_lam_ci <- rep(0, times = n_new)
      if (settings$return_all){
        if(settings$interval == "pi" ||settings$interval == "both")
          s2_y <- s2_lam <- matrix(nrow = n_new, ncol = object$nmcmc)
        if(settings$interval == "ci" ||settings$interval == "both")
          s2_y_ci <- s2_lam_ci <- matrix(nrow = n_new, ncol = object$nmcmc) 
      } 
    } else{
      if(settings$interval == "pi" ||settings$interval == "both")
        sigma_sum_y <- sigma_sum_lam <- matrix(0, nrow = n_new, ncol = n_new)
      if(settings$interval == "ci" ||settings$interval == "both")
        sigma_sum_y_ci <- sigma_sum_lam_ci <- matrix(0, nrow = n_new, ncol = n_new)
    } 
    for (t in 1:object$nmcmc) {

      plam <- gp_vec_lam(object$llam_samples[t, ], x_approx, x_new, g = object$g[t],
                     theta = object$theta_lam[t, ], tau2_0 = object$tau2_lam[t], mean0 = settings$mean_lam,
                     lite = settings$lite, NNarray_pred = NN_x_new,
                     sep = settings$sep, v = settings$v, m = object$m_pred)

      mean_lam[, t] <- plam$mean
        
      py <- gp_vec_y(object$y, x_approx, x_new, lam = exp(object$llam_samples[t, ]), 
                     lam_new = ifel(settings$ub, exp(plam$mean + qnorm(0.95) * plam$s2),exp(plam$mean)),
                     theta = object$theta_y[t, ], tau2 =object$tau2[t], A = object$A,
                     mu_y = settings$mean_y, lite = settings$lite, NNarray_pred = NN_x_new,
                     sep = settings$sep, v = settings$v, m = object$m_pred)
      
      mean_y[, t] <- py$mean

      if (settings$lite) {
        if(settings$interval == "pi" ||settings$interval == "both"){
          s2_sum_lam <- s2_sum_lam + plam$s2
          s2_sum_y <- s2_sum_y + py$s2
        }
        if(settings$interval == "ci" ||settings$interval == "both"){
          s2_sum_lam_ci <- s2_sum_lam_ci + plam$s2_ci
          s2_sum_y_ci <- s2_sum_y_ci + py$s2_ci
        }
        if (settings$return_all){
          if(settings$interval == "pi" ||settings$interval == "both"){
            s2_lam[, t] <- plam$s2
            s2_y[, t] <- py$s2 # storing individual variances for each posterior sample
          }
          if(settings$interval == "ci" ||settings$interval == "both"){
            s2_lam_ci[, t] <- plam$s2_ci
            s2_y_ci[, t] <- py$s2_ci
          }
        } 
      } else{
        if(settings$interval == "pi" ||settings$interval == "both"){
          sigma_sum_y <- sigma_sum_y + py$sigma # storing main covariance matrix
          sigma_sum_lam <- sigma_sum_lam + plam$sigma 
        }
        if(settings$interval == "ci" ||settings$interval == "both"){
          sigma_sum_y_ci <- sigma_sum_y_ci + py$sigma_ci # storing main covariance matrix
          sigma_sum_lam_ci <- sigma_sum_lam_ci + plam$sigma_ci
        }
      }
    } # end of t for loop

  } else { # run in parallel using foreach

    if(settings$verb) print("in parallel")

    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, settings$cores, labels = FALSE)))
    if (settings$cores > detectCores())
      warning("cores is greater than available nodes")

    cl <- makeCluster(settings$cores)
    registerDoParallel(cl)

    thread <- c()
    result <- foreach(thread = 1:settings$cores) %dopar% {
      out <- list()
      out$mean_lam <- out$mean_y <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      if (settings$lite) {
        if(settings$interval == "pi" ||settings$interval == "both")
          out$s2_sum_lam <- out$s2_sum_y <- rep(0, times = n_new)
        if(settings$interval == "ci" ||settings$interval == "both")
          out$s2_sum_lam_ci <-  out$s2_sum_y_ci <- rep(0, times = n_new)
        if (settings$return_all){
          if(settings$interval == "pi" ||settings$interval == "both")
            out$s2_lam <- out$s2_y <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
          if(settings$interval == "ci" ||settings$interval == "both")
            out$s2_lam_ci <- out$s2_y_ci <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
        } 
      } else{
        if(settings$interval == "pi" ||settings$interval == "both")
          out$sigma_sum_lam <- out$sigma_sum_y <- matrix(0, nrow = n_new, ncol = n_new)
        if(settings$interval == "ci" ||settings$interval == "both")
          out$sigma_sum_lam_ci <- out$sigma_sum_y_ci <- matrix(0, nrow = n_new, ncol = n_new)
      } 

      # llam_draw <- yp_draw <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      j <- 1
      for (t in chunks[[thread]]) {

        plam <- gp_vec_lam(object$llam_samples[t, ], x_approx, x_new, g = object$g[t],
                           theta = object$theta_lam[t, ], tau2_0 = object$tau2_lam[t],
                           mean0 = settings$mean_lam, lite = settings$lite, NNarray_pred = NN_x_new,
                           sep = settings$sep, v = settings$v, m = object$m_pred)

        out$mean_lam[, j] <- plam$mean

        # can use llam_draw[, j] but will take time.
        py <- gp_vec_y(object$y, x_approx, x_new, lam = exp(object$llam_samples[t, ]), 
                       lam_new = ifel(settings$ub, exp(plam$mean + qnorm(0.95) * plam$s2),exp(plam$mean)), # exp(out$mean_lam[, j]),
                       theta = object$theta_y[t, ], tau2 =object$tau2[t], A = object$A,
                       mu_y = settings$mean_y, lite = settings$lite, NNarray_pred = NN_x_new,
                       sep = settings$sep, v = settings$v, m = object$m_pred)

        out$mean_y[, j] <- py$mean

        if (settings$lite) {
          if(settings$interval == "pi" ||settings$interval == "both"){
            out$s2_sum_lam <- out$s2_sum_lam + plam$s2
            out$s2_sum_y <- out$s2_sum_y + py$s2
          }
          if(settings$interval == "ci" ||settings$interval == "both"){
            out$s2_sum_lam_ci <- out$s2_sum_lam_ci + plam$s2_ci
            out$s2_sum_y_ci <- out$s2_sum_y_ci + py$s2_ci
          }
          if (settings$return_all){
            if(settings$interval == "pi" ||settings$interval == "both"){
              out$s2_lam[, j] <- plam$s2
              out$s2_y[, j] <- py$s2 # storing individual variances for each posterior sample
            }
            if(settings$interval == "ci" ||settings$interval == "both"){
              out$s2_lam_ci[, j] <- plam$s2_ci
              out$s2_y_ci[, j] <- py$s2_ci # storing individual variances for each posterior sample
            }
          } 
        } else{
          if(settings$interval == "pi" ||settings$interval == "both"){
            out$sigma_sum_lam <- out$sigma_sum_lam + plam$sigma
            out$sigma_sum_y <- out$sigma_sum_y + py$sigma # storing main covariance matrix
          }
          if(settings$interval == "ci" ||settings$interval == "both"){
            out$sigma_sum_lam_ci <- out$sigma_sum_lam_ci + plam$sigma_ci
            out$sigma_sum_y_ci <- out$sigma_sum_y_ci + py$sigma_ci # storing main covariance matrix
          }
        } 
        j <- j + 1
      } # end of t for loop
      return(out)
    } # end of foreach loop

    stopCluster(cl)

    mean_lam <- do.call(cbind, lapply(result, with, eval(parse(text = "mean_lam"))))
    mean_y <- do.call(cbind, lapply(result, with, eval(parse(text = "mean_y"))))

    if (settings$lite) {
      if(settings$interval == "pi" ||settings$interval == "both"){
        s2_sum_lam <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum_lam"))))
        s2_sum_y <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum_y"))))
      }
      if(settings$interval == "ci" ||settings$interval == "both"){
        s2_sum_lam_ci <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum_lam_ci"))))
        s2_sum_y_ci <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum_y_ci"))))
      }
      if (settings$return_all){
        if(settings$interval == "pi" ||settings$interval == "both"){
          s2_lam <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_lam"))))
          s2_y <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_y"))))
        }
        if(settings$interval == "ci" ||settings$interval == "both"){
          s2_lam_ci <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_lam_ci"))))
          s2_y_ci <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_y_ci"))))
        }
      } 
    } else{
      if(settings$interval == "pi" ||settings$interval == "both"){
        sigma_sum_lam <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum_lam"))))
        sigma_sum_y <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum_y"))))
      }
      if(settings$interval == "ci" ||settings$interval == "both"){
        sigma_sum_lam_ci <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum_lam_ci"))))
        sigma_sum_y_ci <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum_y_ci"))))
      }
    } 
  } # end of else statement

  # Add variables to the output list
  object$mean_lnugs <- rowMeans(mean_lam)
  object$mean <- rowMeans(mean_y)
  
  if (settings$lite) {
    if(settings$interval == "pi" ||settings$interval == "both"){
      object$s2_lnugs <- s2_sum_lam / object$nmcmc + apply(mean_lam, 1, var)
      object$s2_y <- s2_sum_y / object$nmcmc + apply(mean_y, 1, var)
    }
    if(settings$interval == "ci" ||settings$interval == "both"){
      object$s2_lnugs_ci <- s2_sum_lam_ci / object$nmcmc + apply(mean_lam, 1, var)
      object$s2_y_ci <- s2_sum_y_ci / object$nmcmc + apply(mean_y, 1, var)
    }
    if (settings$return_all) {
      object$mean_lnugs_all <- mean_lam
      object$mean_all <- mean_y
      if(settings$interval == "pi" ||settings$interval == "both"){
        object$s2_lnugs_all <- s2_lam
        object$s2_all <- s2_y
      }
      if(settings$interval == "ci" ||settings$interval == "both"){
        object$s2_lnugs_all_ci <- s2_lam_ci
        object$s2_all_ci <- s2_y_ci
      }
    }
  } else{
    if(settings$interval == "pi" ||settings$interval == "both"){
      object$Sigma_nugs <- sigma_sum_lam / object$nmcmc + cov(t(mean_lam))
      object$Sigma <- sigma_sum_y / object$nmcmc + cov(t(mean_y))
    }
    if(settings$interval == "ci" ||settings$interval == "both"){
      object$Sigma_nugs_ci <- sigma_sum_lam_ci / object$nmcmc + cov(t(mean_lam))
      object$Sigma_ci <- sigma_sum_y_ci / object$nmcmc + cov(t(mean_y))
    }
  } 

  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  
  if(!settings$sample){
    del <- c( "x", "y", "ys2", "nmcmc", "reps", "settings", "mappings", "m", "x_approx",
              "theta_lam", "theta_y", "llam_samples", "g", "tau2", "tau2_lam", "llik_y", "llik_lam")
    object <- object[!names(object) %in% del]
    return(object)
  } 
  
  return(object)
}

predict_vec_vdims <- function(object, x_new, m = object$m, settings, mapping = object$mappings) {

  tic <- proc.time()[[3]]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  object$m_pred <- m

  if(!is.null(mapping$reps_vdims)){
    xv <- mapping$reps_vdims$X0
    Av <- mapping$reps_vdims$mult
    map <- mapping$reps_vdims$Z

    reps_vdims_new <- find_reps(x_new[, settings$vdims], 1:n_new)
    xv_new <- reps_vdims_new$X0
    Av_new <- reps_vdims_new$mult
    map_new <- reps_vdims_new$Z
  }else{
    xv <- object$xv
    xv_new <- x_new[, settings$vdims, drop = FALSE]
  }

  nv_new <- nrow(xv_new)

  if(nrow(xv) > 300) {
    vecchia_var <- TRUE
    if (settings$lite){
      NN_xv_new <- FNN::get.knnx(object$xv, xv_new, m)$nn.index
      xv_n <- object$xv_approx
    }
    else{
      NN_xv_new <- NULL
      xv_n <- add_pred_to_approx(object$xv_approx, xv_new, m, settings$ordering_new)
    }
  }
  else vecchia_var <- FALSE

  if (settings$return_all & !settings$lite)
    stop("return_all only offered when lite = TRUE")

  if (!is.null(settings$ordering_new)) {
    if (settings$lite) message("ordering_new is only relevant when lite = FALSE")
    test <- check_ordering(settings$ordering_new, nrow(x_new))
  }

  # Pre-calculate nearest neighbors for x
  if (settings$lite) {
    NN_x_new <- FNN::get.knnx(object$x, x_new, m)$nn.index
    x_approx <-  object$x_approx
  } else {
    NN_x_new <- NULL
    x_approx <- add_pred_to_approx(object$x_approx, x_new, m, settings$ordering_new)
  }

  if (settings$cores == 1) { # run serial for loop

    mean_y <- matrix(nrow = n_new, ncol = object$nmcmc)
    # yp_draw <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (settings$lite) {
      s2_sum_y <- rep(0, times = n_new)
      if (settings$return_all) s2_y <- matrix(nrow = n_new, ncol = object$nmcmc)
    } else sigma_sum_y <- matrix(0, nrow = n_new, ncol = n_new)
    
    mean_lam <- matrix(nrow = n_new, ncol = object$nmcmc)
    llam_draw <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (settings$lite) {
      s2_sum_lam <- rep(0, times = n_new)
      if (settings$return_all) s2_lam <- matrix(nrow = n_new, ncol = object$nmcmc)
    } else sigma_sum_lam <- matrix(0, nrow = n_new, ncol = n_new)

    for (t in 1:object$nmcmc) {

      if(vecchia_var){
        plam <- gp_vec_lam(object$llam_samples[t, ], xv_n, xv_new, g = object$g[t], theta = object$theta_lam[t, ],
                           mean0 = settings$mean_lam, tau2_0 = object$tau2_lam[t],
                           lite = settings$lite, NNarray_pred = NN_xv_new,
                           sep = settings$sep, v = settings$v, m = object$m_pred)
      }else{
        plam <- gp_lam(object$llam_samples[t, ], xv, xv_new, g = object$g[t],
                       theta = object$theta_lam[t, ], mu_0 = settings$mean_lam, 
                       tau2_0 =  object$tau2_lam[t], sep = settings$sep, 
                       v = settings$v, lite = settings$lite)
      }

      if(!is.null(mapping$reps_vdims)){
        mean_np_lam[map_new] <- rep(plam$mean, Av_new)
        s2_np_lam[map_new] <- rep(plam$s2, Av_new)
        llam_iter <- ifel(settings$ub, mean_np_lam + qnorm(0.95) * s2_np_lam, mean_np_lam)
        llam_n[map] <- rep(object$llam_samples[t, ], Av)

      }else{
        mean_np_lam <- plam$mean
        s2_np_lam <- plam$s2
        llam_iter <- ifel(settings$ub, plam$mean + qnorm(0.95) * plam$s2, plam$mean)
        llam_n <- object$llam_samples[t, ]
      }

      mean_lam[, t] <- plam$mean

      if (settings$lite) {
        s2_sum_lam <- s2_sum_lam + plam$s2
        if (settings$return_all) s2_lam[, t] <- plam$s2
      } else sigma_sum_lam <- sigma_sum_lam + plam$sigma

      py <- gp_vec_y(object$y, x_approx, x_new, lam = exp(llam_n), 
                     lam_new = exp(llam_iter), # exp(llam_draw[, t]),
                     theta = object$theta_y[t, ], tau2 = object$tau2[t], lite = settings$lite,
                     NNarray_pred = NN_x_new, mu_y = settings$mean_y,
                     sep = settings$sep, v = settings$v, m = object$m_pred)
      
      mean_y[, t] <- py$mean
      if (settings$lite) {
        s2_sum_y <- s2_sum_y + py$s2
        if (settings$return_all) s2_y[, t] <- py$s2 # storing individual variances for each posterior sample
      } else sigma_sum_y <- sigma_sum_y + py$sigma # storing main covariance matrix
      
      # yp_draw[, t] <- mvtnorm::rmvnorm(1, mean = mean_y[, t], sigma = diag (py$s2))
    } # end of t for loop

  } else { # run in parallel using foreach

    if(settings$verb) print("in parallel")

    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, settings$cores, labels = FALSE)))
    if (settings$cores > detectCores())
      warning("cores is greater than available nodes")

    cl <- makeCluster(settings$cores)
    registerDoParallel(cl)

    thread <- NULL
    result <- foreach(thread = 1:settings$cores) %dopar% {
      out <- list()
      out$mean_y <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      out$mean_lam <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      if (settings$lite) {
        out$s2_sum_y <- rep(0, times = n_new)
        out$s2_sum_lam <- rep(0, times = n_new)
        if (settings$return_all) {
          out$s2_y <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
          out$s2_lam <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
        }
      } else {
        out$sigma_sum_y <- matrix(0, nrow = n_new, ncol = n_new)
        out$sigma_sum_lam <- matrix(0, nrow = n_new, ncol = n_new)
      }

      # llam_draw <- yp_draw <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      j <- 1
      for (t in chunks[[thread]]) {
        
        if(vecchia_var){
          plam <- gp_vec_lam(object$llam_samples[t, ], xv_n, xv_new, g = object$g[t], theta = object$theta_lam[t, ],
                             mean0 = settings$mean_lam, tau2_0 = object$tau2_lam[t],
                             lite = settings$lite, NNarray_pred = NN_xv_new,
                             sep = settings$sep, v = settings$v, m = object$m_pred)
        }else{
          plam <- gp_lam(object$llam_samples[t, ], xv, xv_new, g = object$g[t],
                         theta = object$theta_lam[t, ], mu_0 = settings$mean_lam, 
                         tau2_0 =  object$tau2_lam[t], sep = settings$sep, 
                         v = settings$v, lite = settings$lite)
        }
        
        if(!is.null(mapping$reps_vdims)){
          mean_np_lam[map_new] <- rep(plam$mean, Av_new)
          s2_np_lam[map_new] <- rep(plam$s2, Av_new)
          llam_iter <- ifel(settings$ub, mean_np_lam + qnorm(0.95) * s2_np_lam, mean_np_lam)
          llam_n[map] <- rep(object$llam_samples[t, ], Av)
          
        }else{
          mean_np_lam <- plam$mean
          s2_np_lam <- plam$s2
          llam_iter <- ifel(settings$ub, plam$mean + qnorm(0.95) * plam$s2, plam$mean)
          llam_n <- object$llam_samples[t, ]
        }
        
        mean_lam[, j] <- plam$mean
        
        if (settings$lite) {
          s2_sum_lam <- s2_sum_lam + plam$s2
          if (settings$return_all) s2_lam[, j] <- plam$s2
        } else sigma_sum_lam <- sigma_sum_lam + plam$sigma
        
        py <- gp_vec_y(object$y, x_approx, x_new, lam = exp(llam_n), 
                       lam_new = exp(llam_iter), # exp(llam_draw[, t]),
                       theta = object$theta_y[t, ], tau2 = object$tau2[t], lite = settings$lite,
                       NNarray_pred = NN_x_new, mu_y = settings$mean_y,
                       sep = settings$sep, v = settings$v, m = object$m_pred)
        
        mean_y[, j] <- py$mean
        if (settings$lite) {
          s2_sum_y <- s2_sum_y + py$s2
          if (settings$return_all) s2_y[, j] <- py$s2 # storing individual variances for each posterior sample
        } else sigma_sum_y <- sigma_sum_y + py$sigma # storing main covariance matrix
        
        j <- j + 1
        # yp_draw[, t] <- mvtnorm::rmvnorm(1, mean = mean_y[, t], sigma = diag (py$s2))
      } # end of t for loop
    } # end of foreach loop

    stopCluster(cl)
    
    # Group elements out of the list
    mean_y <- do.call(cbind, lapply(result, with, eval(parse(text = "mean_y"))))
    mean_lam <- do.call(cbind, lapply(result, with, eval(parse(text = "mean_lam"))))

    if (settings$lite) {
      s2_sum_y <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum_y"))))
      s2_sum_lam <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum_lam"))))
      if (settings$return_all) {
        s2_y <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_y"))))
        s2_lam <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_lam"))))
      }
    } else {
      sigma_sum_y <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum_y"))))
      sigma_sum_lam <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum_lam"))))
    }
  } # end of else statement

  # Add variables to the output list
  object$mean <- rowMeans(mean_y)
  if (settings$lite) {
    object$s2_y <- s2_sum_y / object$nmcmc + apply(mean_y, 1, var)
    if (settings$return_all) {
      object$mean_all <- mean_y
      object$s2_all <- s2_y
    }
  } else object$Sigma <- sigma_sum_y / object$nmcmc + cov(t(mean_y))
  
  object$mean_lnugs <- rowMeans(mean_lam)
  if (settings$lite) {
    object$s2_lnugs <- s2_sum_lam / object$nmcmc + apply(mean_lam, 1, var)
    if (settings$return_all) {
      object$mean_lnugs_all <- mean_lam
      object$s2_lnugs_all <- s2_lam
    }
  } else object$Sigma_nugs <- sigma_sum_lam / object$nmcmc + cov(t(mean_lam))

  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)

  return(object)
}

predict_vec_hom <- function(object, x_new, m = object$m, settings) {

  tic <- proc.time()[[3]]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  object$m_pred <- m

  if (settings$return_all & !settings$lite)
    stop("return_all only offered when lite = TRUE")

  if (!is.null(settings$ordering_new)) {
    if (settings$lite) message("ordering_new is only relevant when lite = FALSE")
    test <- check_ordering(settings$ordering_new, nrow(x_new))
  }

  # Pre-calculate nearest neighbors for x
  if (settings$lite) {
    NN_x_new <- FNN::get.knnx(object$x, x_new, m)$nn.index
    x_approx <-  object$x_approx
  } else {
    NN_x_new <- NULL
    x_approx <- add_pred_to_approx(object$x_approx, x_new, m, settings$ordering_new)
  }

  if (settings$cores == 1) { # run serial for loop

    mean_y <- matrix(nrow = n_new, ncol = object$nmcmc)

    if (settings$lite) {
      if(settings$interval == "pi" ||settings$interval == "both")
        s2_sum_y <- rep(0, times = n_new)
      if(settings$interval == "ci" ||settings$interval == "both")
        s2_sum_y_ci <- rep(0, times = n_new)
      if (settings$return_all){
        s2_y <- matrix(nrow = n_new, ncol = object$nmcmc)
        if(settings$interval == "ci" ||settings$interval == "both")
          s2_y_ci <- matrix(nrow = n_new, ncol = object$nmcmc)
      } 
    } else {
        if(settings$interval == "pi" ||settings$interval == "both")
          sigma_sum_y <- matrix(0, nrow = n_new, ncol = n_new)
        if(settings$interval == "ci" ||settings$interval == "both")
          sigma_sum_y_ci <- matrix(0, nrow = n_new, ncol = n_new)
      }
    # yp_draw <- matrix(nrow = n_new, ncol = object$nmcmc)

    for (t in 1:object$nmcmc) {

      py <- gp_vec_hom(object$y, x_approx, x_new, A = object$A, nugs = object$g_y[t],
                       theta = object$theta_y[t, ], tau2 = object$tau2[t], lite = settings$lite,
                       NNarray_pred = NN_x_new, sep = settings$sep, v = settings$v, m = object$m_pred,
                       mu_y = settings$mean_y)

      mean_y[, t] <- py$mean
      if (settings$lite) {
        if(settings$interval == "pi" ||settings$interval == "both")
          s2_sum_y <- s2_sum_y + py$s2
        if(settings$interval == "ci" ||settings$interval == "both")
          s2_sum_y_ci <- s2_sum_y_ci + py$s2_ci
        if (settings$return_all){
          if(settings$interval == "pi" ||settings$interval == "both")
            s2_y[, t] <- py$s2 # storing individual variances for each posterior sample
          if(settings$interval == "ci" ||settings$interval == "both")
            s2_y_ci[, t] <- py$s2_ci
        } 
      } else{
        if(settings$interval == "pi" ||settings$interval == "both")
          sigma_sum_y <- sigma_sum_y + py$sigma # storing main covariance matrix
        if(settings$interval == "ci" ||settings$interval == "both")
          sigma_sum_y_ci <- sigma_sum_y_ci + py$sigma_ci # storing main covariance matrix
      } 

    } # end of t for loop

  } else { # run in parallel using foreach
    
    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, settings$cores, labels = FALSE)))
    if (settings$cores > detectCores())
      warning("cores is greater than available nodes")

    cl <- makeCluster(settings$cores)
    registerDoParallel(cl)

    thread <- NULL
    result <- foreach(thread = 1:settings$cores) %dopar% {
      out <- list()
      out$mean_y <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      if (settings$lite) {
        if(settings$interval == "pi" ||settings$interval == "both")
          out$s2_sum_y <- rep(0, times = n_new)
        if(settings$interval == "ci" ||settings$interval == "both")
          out$s2_sum_y_ci <- rep(0, times = n_new)
        if (settings$return_all) {
          if(settings$interval == "pi" ||settings$interval == "both")
            out$s2_y <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
          if(settings$interval == "ci" ||settings$interval == "both")
            out$s2_y_ci <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
        }
      } else {
        if(settings$interval == "pi" ||settings$interval == "both")
          out$sigma_sum_y <- matrix(0, nrow = n_new, ncol = n_new)
        if(settings$interval == "ci" ||settings$interval == "both")
          out$sigma_sum_y_ci <- matrix(0, nrow = n_new, ncol = n_new)
      }

      j <- 1
      for (t in chunks[[thread]]) {

        py <- gp_vec_hom(object$y, x_approx, x_new, A = object$A, nugs = object$g_y[t],
                         theta = object$theta_y[t, ], tau2 = object$tau2[t], lite = settings$lite,
                         NNarray_pred = NN_x_new, sep = settings$sep, v = settings$v, m = object$m_pred,
                         mu_y = settings$mean_y)
        
        out$mean_y[, j] <- py$mean
        if (settings$lite) {
          if(settings$interval == "pi" ||settings$interval == "both")
            out$s2_sum_y <- out$s2_sum_y + py$s2
          if(settings$interval == "ci" ||settings$interval == "both")
            out$s2_sum_y_ci <- out$s2_sum_y_ci + py$s2_ci
          if (settings$return_all){
            if(settings$interval == "pi" ||settings$interval == "both")
              out$s2_y[, j] <- py$s2
            if(settings$interval == "ci" ||settings$interval == "both")
              out$s2_y_ci[, j] <- py$s2_ci
          }  # storing individual variances for each posterior sample
        } else{
          if(settings$interval == "pi" ||settings$interval == "both")
            out$sigma_sum_y <- out$sigma_sum_y + py$sigma # storing main covariance matrix
          if(settings$interval == "ci" ||settings$interval == "both")
            out$sigma_sum_y_ci <- out$sigma_sum_y_ci + py$sigma_ci
        }  # storing main covariance matrix

        # yp_draw[, j] <- rmvnorm(1, mean = py$mean, sigma = diag (py$s2))
        j <- j + 1
      } # end of t for loop
      return(out)
    } # end of foreach loop

    stopCluster(cl)
    mean_y <- do.call(cbind, lapply(result, with, eval(parse(text = "mean_y"))))

    if (settings$lite) {
      if(settings$interval == "pi" ||settings$interval == "both")
        s2_sum_y <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum_y"))))
      if(settings$interval == "ci" ||settings$interval == "both")
        s2_sum_y_ci <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum_y_ci"))))
      if (settings$return_all) {
        if(settings$interval == "pi" ||settings$interval == "both")
          s2_y <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_y"))))
        if(settings$interval == "ci" ||settings$interval == "both")
          s2_y_ci <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_y_ci"))))
      }
    } else {
      if(settings$interval == "pi" ||settings$interval == "both")
        sigma_sum_y <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum_y"))))
      if(settings$interval == "ci" ||settings$interval == "both")
        sigma_sum_y_ci <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum_y_ci"))))
    }
  } # end of else statement

  # Add variables to the output list
  object$mean <- rowMeans(mean_y)

  if (settings$lite) {
    if(settings$interval == "pi" ||settings$interval == "both")
      object$s2_y <- s2_sum_y / object$nmcmc + apply(mean_y, 1, var)
    if(settings$interval == "ci" ||settings$interval == "both")
      object$s2_y_ci <- s2_sum_y_ci / object$nmcmc + apply(mean_y, 1, var)
    if (settings$return_all) {
      object$mean_all <- mean_y
      if(settings$interval == "pi" ||settings$interval == "both")
        object$s2_all <- s2_y
      if(settings$interval == "ci" ||settings$interval == "both")
        object$s2_all_ci <- s2_y_ci
    }
  } else{
    if(settings$interval == "pi" ||settings$interval == "both")
      object$Sigma <- sigma_sum_y / object$nmcmc + cov(t(mean_y))
    if(settings$interval == "ci" ||settings$interval == "both")
      object$Sigma_ci <- sigma_sum_y_ci / object$nmcmc + cov(t(mean_y))
  } 

  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)

  return(object)
}

gp_vec_lam <- function(out_vec, x_approx = NULL, x_new, g = 1e-5,
                       theta, tau2_0 = 1, mean0, a0, b0, lite = TRUE,
                       NNarray_pred = NULL, sep = FALSE, v, m = NULL){

  # g_vec <- rep(g, length(out_vec) + nrow(x_new))

  if(!lite){

    g_vec <- rep(g, length(out_vec) + nrow(x_new))
    out_vec_o <- out_vec[x_approx$ord[x_approx$observed]] - mean0
    U_mat <- create_U(x_approx, g_vec, theta = theta, v= v, sep = sep)# obs U
    Upp <- U_mat[!x_approx$observed, !x_approx$observed] # pred U
    Uppinv <- Matrix::solve(Upp, sparse = TRUE) # pred U inv

    Winv <- Matrix::crossprod(Uppinv) #  Upp inv T Upp in

    Uop <- U_mat[x_approx$observed, !x_approx$observed] # U obs X U pred
    UopT_out_vec <- Matrix::crossprod(Uop, out_vec_o) # Uop transpoe Y
    mu_ordered <- -Matrix::crossprod(Uppinv, UopT_out_vec) # - Up T (inv) Uop  T Y

    mean <- mean0 + mu_ordered[x_approx$rev_ord_pred] # OG order back
    
    sigma_ci <- sigma <- as.matrix( (tau2_0) * Winv[x_approx$rev_ord_pred, x_approx$rev_ord_pred])
    diag(sigma_ci) <- diag(sigma) - (tau2_0 * g_vec[-c(1:length(out_vec))])

    return(list(mean = mean, sigma = sigma, s2 = diag(sigma), 
                sigma_ci = sigma_ci, s2_ci = diag(sigma_ci)))
  }

  else{
    n_new <- nrow(x_new)
    # # don't want it ordered in any way. original X's
    Xn <- x_approx$x_ord[x_approx$rev_ord_obs, , drop = FALSE]
    mean <- vector(length = n_new)
    s2 <- s2_ci <- vector(length = n_new)
    g_vec <- rep(g, ncol(NNarray_pred) + 1)
    
    # prior_mean <- rep(0, m)
    for (i in 1:n_new) {
      
      NN <- NNarray_pred[i, ]
      x_combined <- rbind(Xn[NN, , drop = FALSE], x_new[i, , drop = FALSE])
      
      if(v == 999){
        if(sep) K <- Exp2SepVec(x_combined, x_combined, 1, theta, g_vec)
        else K <- Exp2vec(sq_dist(x_combined), 1, theta, g_vec)
      }else if(v > 1000){
        if(sep) K <- MaternProdSepVec(x_combined, x_combined, 1, theta, g_vec, (v- 1000))
        else K <- MaternVec(sq_dist(x_combined), 1, theta, g_vec, (v - 1000))
      }
      else{
        if(sep) K <- MaternSepVec(x_combined, x_combined, 1, theta, g_vec, v)
        else K <- MaternVec(sq_dist(x_combined), 1, theta, g_vec, v)
      }

      L <- t(chol(K))

      mean[i] <- mean0 + L[m + 1, 1:m] %*% forwardsolve(L[1:m, 1:m], out_vec[NN] - mean0)
      s2[i] <- tau2_0 * (L[m + 1, m + 1] ^ 2)
      s2_ci[i] <- s2[i] - (tau2_0 * g)
    }
    
    return(list(mean = mean, s2 = s2, sigma = NULL, s2_ci = s2_ci))
  }
}

gp_vec_y <- function(out_vec, x_approx = NULL, x_new, A = NULL, lam = 1e-5,
                     lam_new = 1e-5, theta, tau2 = 1, mu_y, a, b, lite = TRUE,
                     NNarray_pred = NULL, sep = FALSE, v, m = NULL){

  if(!is.null(A)) lam <- lam/A
  lam_vec <- c(lam, lam_new)

  if(!lite){
    # stop("check here")
    lam_vec <- lam_vec[x_approx$ord]

    out_vec_o <- out_vec[x_approx$ord[x_approx$observed]] - mu_y# tells you which are obs
    U_mat <- create_U(x_approx, lam_vec, theta = theta, v= v, sep = sep)# obs U
    Upp <- U_mat[!x_approx$observed, !x_approx$observed] # pred U
    Uppinv <- Matrix::solve(Upp, sparse = TRUE) # pred U inv

    Winv <- Matrix::crossprod(Uppinv) #  Upp inv T Upp in

    Uop <- U_mat[x_approx$observed, !x_approx$observed] # U obs X U pred
    UopT_out_vec <- Matrix::crossprod(Uop, out_vec_o) # Uop transpoe Y
    mu_ordered <- -Matrix::crossprod(Uppinv, UopT_out_vec) # - Up T (inv) Uop  T Y

    mean <- mu_y + mu_ordered[x_approx$rev_ord_pred] # OG order back
    sigma_ci <- sigma <- as.matrix(tau2 * Winv[x_approx$rev_ord_pred, x_approx$rev_ord_pred])
    diag(sigma_ci) <- diag(sigma) - (tau2 * lam_new)
    
    return(list(mean = mean, sigma = sigma, s2 = diag(sigma), 
                sigma_ci = sigma_ci, s2_ci = diag(sigma_ci)))
  }

  else{
    n_new <- nrow(x_new)
    # # don't want it ordered in any way. original X's
    Xn <- x_approx$x_ord[x_approx$rev_ord_obs, , drop = FALSE]
    mean <- vector(length = n_new)
    s2 <- s2_ci <- vector(length = n_new)

    for (i in 1:n_new) {
      NN <- NNarray_pred[i, ]
      x_combined <- rbind(Xn[NN, , drop = FALSE], x_new[i, , drop = FALSE])
      lam_vec <- c(lam[NN], lam_new[i])
      # eps <- rep(1e-8, length(lam_vec))
      
      if(v == 999){
        if(sep) K <- Exp2SepVec(x_combined, x_combined, 1, theta, lam_vec)
        else K <- Exp2vec(sq_dist(x_combined), 1, theta, lam_vec)
      } else if(v > 1000){
        if(sep) K <- MaternProdSepVec(x_combined, x_combined, 1, theta, lam_vec, (v - 1000))
          #K <- MaternProdSepVec(x_combined, x_combined, 1, theta, lam_vec, (v - 1000))
        else K <- MaternVec(sq_dist(x_combined), 1, theta, lam_vec, (v - 1000))
      }else{
        if(sep) K <- MaternSepVec(x_combined, x_combined, 1, theta, lam_vec, v)
        else K <- MaternVec(sq_dist(x_combined), 1, theta, lam_vec, v)
      }

      L <- t(chol(K))

      mean[i] <- mu_y + L[m + 1, 1:m] %*% forwardsolve(L[1:m, 1:m], out_vec[NN] - mu_y)
      s2[i] <- tau2 * (L[m + 1, m + 1] ^ 2)
      s2_ci[i] <- s2[i] -  (tau2 *  lam_new[i])
    }
    # stop("why is this inflated")
    return(list(mean = mean, s2 = s2, sigma = NULL, s2_ci = s2_ci))
  }
}


gp_vec_hom <- function(out_vec, x_approx = NULL, x_new, A = NULL, nugs = 1e-5,
                       theta, tau2 = 1, a, b, lite = TRUE,
                       NNarray_pred = NULL, sep = FALSE, v, m = NULL, mu_y = NULL){
  
  if(!is.null(A)) g <- nugs/A
  else g <- rep(nugs, nrow(x_approx$x_ord))

  g_new <- rep(nugs, nrow(x_new))
  g_vec <- c(g, g_new)

  if(!lite){
    g_vec <- g_vec[x_approx$ord]

    out_vec_o <- out_vec[x_approx$ord[x_approx$observed]] - mu_y

    U_mat <- create_U(x_approx, g_vec, theta = theta, v=v, sep = sep)# obs U
    Upp <- U_mat[!x_approx$observed, !x_approx$observed] # pred U
    Uppinv <- Matrix::solve(Upp, sparse = TRUE) # pred U inv

    Winv <- Matrix::crossprod(Uppinv) #  Upp inv T Upp in

    Uop <- U_mat[x_approx$observed, !x_approx$observed] # U obs X U pred
    UopT_out_vec <- Matrix::crossprod(Uop, out_vec_o) # Uop transpoe Y
    mu_ordered <- -Matrix::crossprod(Uppinv, UopT_out_vec) # - Up T (inv) Uop  T Y

    mean <- mu_y + mu_ordered[x_approx$rev_ord_pred] # OG order back
    sigma_ci <- sigma <- as.matrix(tau2 * Winv[x_approx$rev_ord_pred, x_approx$rev_ord_pred])
    diag(sigma_ci) <- diag(sigma) - (tau2 * g_new)
    
    return(list(mean = mean, sigma = sigma, s2 = diag(sigma), 
                sigma_ci = sigma_ci, s2_ci = diag(sigma_ci)))
  }

  else{
    n_new <- nrow(x_new)
    # # don't want it ordered in any way. original X's
    Xn <- x_approx$x_ord[x_approx$rev_ord_obs, , drop = FALSE]
    mean <- vector(length = n_new)
    s2_ci <- s2 <- vector(length = n_new)

    prior_mean <- rep(0, m)
    for (i in 1:n_new) {
      NN <- NNarray_pred[i, ]
      x_combined <- rbind(Xn[NN, , drop = FALSE], x_new[i, , drop = FALSE])
      g_vec <- c(g[NN], g_new[i])

      if(v == 999){
        if(sep) K <- Exp2SepVec(x_combined, x_combined, 1, theta, g_vec)
        else K <- Exp2vec(sq_dist(x_combined), 1, theta, g_vec)
      } else if(v > 1000){
        if(sep) K <- MaternSepVec(x_combined, x_combined, 1, theta, g_vec, (v - 1000))
        else K <- MaternVec(sq_dist(x_combined), 1, theta, g_vec, (v - 1000))
      }else{
        if(sep) K <- MaternSepVec(x_combined, x_combined, 1, theta, g_vec, v)
        else K <- MaternVec(sq_dist(x_combined), 1, theta, g_vec, v)
      }

      L <- t(chol(K))

      mean[i] <- mu_y + L[m + 1, 1:m] %*% forwardsolve(L[1:m, 1:m], out_vec[NN] - mu_y)
      s2[i] <- tau2 * (L[m + 1, m + 1] ^ 2)
      s2_ci[i] <- s2[i] - (tau2 * g_new[i])
    }
    return(list(mean = mean, s2 = s2, sigma = NULL, s2_ci = s2_ci))
  }
}
