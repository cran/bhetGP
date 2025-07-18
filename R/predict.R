# -------- Functions ----------------------------------------------------------
# pred_bhetGP: for bhetgp objects (vec or non vec) -> calls function
# pred_bhomGP: for bhomgp objects vec or non vec -> calls function
# predict_nonvec: non vecchia predictions for hetgp
# predict_nonvec_vdims: non vecchia preds for hetgp w/ vdims
# predict_hom: homGP preds
# gp_lam: calculates pred mean and s2 for lambda layer
# gp_yw: calculates pred mean and s2 for y layer
# gp_hom: calclulates pred mean and s2 for homgp
# -----------------------------------------------------------------------------

# Define Predict for S3 Objects -----------------------------------------------
#' @name predict
#' @title Predict posterior mean and variance/covariance
#' @description Acts on a \code{bhetgp}, \code{bhetgp_vec}, \code{bhomgp} or, 
#'     \code{bhomgp_vec} object. Calculates posterior mean and variance/covariance 
#'     over specified input locations. Optionally utilizes SNOW parallelization.
#' 
#' @details All iterations in the object are used for prediction, so samples 
#'     should be burned-in. Thinning the samples using \code{trim} will speed 
#'     up computation. Posterior moments are calculated using conditional 
#'     expectation and variance. As a default, only point-wise variance is 
#'     calculated. Full covariance may be calculated using \code{lite = FALSE}. 
#'     
#'     The posterior predictive variances are returned by default. The variance 
#'     for the mean process may be obtained by specifying \code{interval = "ci"}.
#'     \code{interval = "both"} will return both variances.
#'     
#'     SNOW parallelization reduces computation time but requires 
#'     more memory storage.
#' 
#' @param object object from \code{bhetGP} or \code{bhomGP}, 
#'        with burn-in already removed
#' @param x_new matrix of predictive input locations
#' @param lite logical indicating whether to calculate only point-wise 
#'        variances (\code{lite = TRUE}) or full covariance 
#'        (\code{lite = FALSE})
#' @param return_all logical indicating whether to return mean and point-wise
#'        variance prediction for ALL samples (only available for \code{lite = TRUE})
#' @param interval returns predictive variances by default \code{interval = "pi"}.
#'        \code{interval = "ci"} returns variances for only mean process and 
#'        \code{interval = "both"} returns both variances.
#' @param lam_ub logical uses upper 95 quantile for latent noise to 
#'        obtain predictive variances for the response. If \code{lam_ub = FALSE},
#'        the mean latent noise is used for inference.
#' @param vecchia logical uses vecchia approximation for prediction if \code{vecchia = TRUE}.
#' @param m size of Vecchia conditioning sets (only for fits with 
#'        \code{vecchia = TRUE}), defaults to the \code{m} used for MCMC
#' @param ordering_new optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x_new}, defaults to random, is applied to all layers
#'        in deeper models.
#' @param cores number of cores to utilize for SNOW parallelization
#' @param omp_cores sets cores used for OpenMP if \code{vechhia = TRUE} and 
#'        \code{lite = FALSE}. 
#' @param samples logical indicating if you want all posterior samples returned 
#'        including latent layer.
#' @param ... N/A
#'        
#' @return object of the same class with the following additional elements:
#' \itemize{
#'   \item \code{x_new}: copy of predictive input locations
#'   \item \code{m_pred}: size of predictive conditioning set if \code{vecchia = TRUE}
#'   \item \code{mean_y}: predicted posterior mean, indices correspond to 
#'         \code{x_new} locations
#'   \item \code{s2_y}: predicted point-wise variances, indices correspond to 
#'         \code{x_new} locations (only returned when \code{lite = TRUE} & 
#'         \code{interval = c("pi", "both")})
#'   \item \code{s2_y_ci}: predicted point-wise variances for the mean process, 
#'         indices correspond to \code{x_new} locations 
#'         (only returned when \code{lite = TRUE} & \code{interval = c("ci", "both")})
#'   \item \code{mean_all}: predicted posterior mean for each sample (column
#'         indices), only returned when \code{return_all = TRUE} 
#'   \item \code{s2_all} predicted point-wise variances each sample (column
#'         indices), only returned when \code{return-all = TRUE} & \code{interval = c("pi", "both")}
#'   \item \code{s2_all_ci} predicted point-wise variances for each sample (column
#'         indices), only returned when \code{return-all = TRUE} & \code{interval = c("ci", "both")}
#'   \item \code{Sigma}: predicted posterior covariance, indices correspond to 
#'         \code{x_new} locations (only returned when \code{lite = FALSE} 
#'         & \code{interval = c("pi", "both")})
#'   \item \code{Sigma_ci}: predicted posterior covariance for mean process, 
#'         indices correspond to \code{x_new} locations 
#'         (only returned when \code{lite = FALSE} & \code{interval = c("ci", "both")})
#' }
#' Additionally, if object belongs to class \code{bhetGP} or \code{bhetGP_vec}, the
#' log-noise process is also predicted for new locations \code{x_new}. The following are returned:
#' \itemize{
#'   \item \code{mean_lnugs}: predicted posterior mean for log noise process, 
#'         indices correspond to \code{x_new} locations
#'   \item \code{s2_lnugs}: predicted point-wise variances for log noise process, 
#'         indices correspond to \code{x_new} locations (only returned when 
#'         \code{lite = TRUE} & \code{interval = c("pi", "both")})
#'   \item \code{s2_lnugs_ci}: predicted point-wise variances for the log noise process, 
#'         indices correspond to \code{x_new} locations 
#'         (only returned when \code{lite = TRUE} & \code{interval = c("ci", "both")})
#'   \item \code{mean_lnugs_all}: predicted posterior mean for each sample for log 
#'         noise process (column indices), only returned when \code{return_all = TRUE} 
#'   \item \code{s2_lnugs_all} predicted point-wise variances each sample (column
#'         indices) for log noise process, only returned when \code{return-all = TRUE} 
#'         & \code{interval = c("pi", "both")}
#'   \item \code{s2_lnugs_all_ci} predicted point-wise variances for each sample (column
#'         indices) for log noise process, only returned when \code{return-all = TRUE} 
#'         & \code{interval = c("ci", "both")}
#'   \item \code{Sigma_lnugs}: predicted posterior covariance for log noise process, 
#'         indices correspond to \code{x_new} locations (only returned when \code{lite = FALSE} 
#'         & \code{interval = c("pi", "both")})
#'   \item \code{Sigma_lnugs_ci}: predicted posterior covariance for log noise process, 
#'         indices correspond to \code{x_new} locations 
#'         (only returned when \code{lite = FALSE} & \code{interval = c("ci", "both")})
#'  }
#' Computation time is added to the computation time of the existing object.
#' 
#' @rdname predict
NULL

#' @rdname predict
#' @export
predict.bhetgp <- function(object, x_new, lite = TRUE, return_all = FALSE, 
                           interval = c("pi", "ci", "both"), lam_ub = TRUE,
                           cores = 1, samples = TRUE, ...){

  # ub_lam --> 95% uber bound of lambda for PI of y.
  interval <- match.arg(interval)
  settings <- c(object$settings, lite = lite, return_all = return_all,
                cores = cores, interval = interval, samples = samples, ub = lam_ub)

  out <- list()

  if(is.null(object$settings$vdims)) out <- predict_nonvec(object, x_new, settings)
  else out <- predict_nonvec_vdims(object, x_new, settings, mapping = object$mappings)

  return(out)
}

#' @rdname predict
#' @export
predict.bhetgp_vec <- function(object, x_new, lite = TRUE, return_all = FALSE, 
                               interval = c("pi", "ci", "both"), lam_ub = TRUE,
                               vecchia = FALSE, m = object$m, ordering_new = NULL,
                               cores = 1, omp_cores = 2, samples = TRUE, ...){
  
  # ub_lam --> 95% uber bound of lambda for PI of y.
  interval <- match.arg(interval)
  settings <- c(object$settings, lite = lite, return_all = return_all, 
                ordering_new = ordering_new, cores = cores, interval = interval,
                samples = samples, ub = lam_ub)
  
  # Checks OpenMP cores for Umat parallelizing
  omp_cores <- check_cores(omp_cores)
  object$x_approx$n_cores <- omp_cores
  out <- list()
  
  if(!vecchia){ # You can do non-vec prediction for vec created object
    if(is.null(object$settings$vdims)) out <- predict_nonvec(object, x_new, settings)
    else out <- predict_nonvec_vdims(object, x_new, settings, mapping = object$mappings)
  }else{
    if(is.null(object$settings$vdims)) out <- predict_vec(object, x_new, m, settings)
    else out <- predict_vec_vdims(object, x_new, m, settings)
  }
  
  return(out)
}

#' @rdname predict
#' @export
predict.bhomgp <- function(object, x_new, lite = TRUE, return_all = FALSE,
                        interval = c("pi", "ci", "both"), cores = 1, samples = TRUE, ...){

  interval <- match.arg(interval)  
  settings <- c(object$settings, lite = lite, return_all = return_all,
                cores = cores,  interval = interval)

  out <- list()
  out <- predict_hom(object, x_new, settings)

  return(out)
}

#' @rdname predict
#' @export
predict.bhomgp_vec <- function(object, x_new, lite = TRUE, return_all = FALSE,
                           interval = c("pi", "ci", "both"),
                           vecchia = FALSE, m = object$m,  ordering_new = NULL, 
                           cores = 1, omp_cores = 2, samples = TRUE, ...){
  
  interval <- match.arg(interval)  
  settings <- c(object$settings, lite = lite, return_all = return_all,
                ordering_new = ordering_new, cores = cores,  interval = interval)
  
  # Checks OpenMP cores for Umat parallelizing
  omp_cores <- check_cores(omp_cores)
  object$x_approx$n_cores <- omp_cores
  out <- list()
  
  # can do non vec predictions for vec object
  if(!vecchia) out <- predict_hom(object, x_new, settings)
  else out <- predict_vec_hom(object, x_new, m, settings)
  
  return(out)
}

predict_nonvec <- function(object, x_new, settings){

  tic <- proc.time()[[3]]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)

  object$x_new <- x_new
  n_new <- nrow(object$x_new)

  if (settings$cores == 1) { # run serial for loop
    
    mean_y <- mean_lam <- matrix(nrow = n_new, ncol = object$nmcmc)
    yp_draw <- llam_draw <- matrix(nrow = n_new, ncol = object$nmcmc)
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

      plam <- gp_lam(object$llam_samples[t, ], object$x , x_new, g = object$g[t],
                     theta = object$theta_lam[t, ], mu_0 = settings$mean_lam, 
                     tau2_0 = object$tau2_lam[t], sep = settings$sep, v = settings$v)

      mean_lam[, t] <- plam$mean

      # draw from the full covariance and then sample llam
      llam_draw[, t] <- drop(rmvnorm(1, plam$mean, diag(plam$s2, nrow = length(plam$mean))))

      py <- gp_yw(object$y, object$x, x_new, A = object$A, lam = exp(object$llam_samples[t, ]),
                  lam_new = exp(llam_draw[, t]), theta = object$theta_y[t, ], mu = settings$mean_y, 
                  tau2 = object$tau2[t], sep = settings$sep, v = settings$v, lite = settings$lite)
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
      # yp_draw[, t] <- rmvnorm(1, mean = mean_y[, t], sigma = diag (py$s2))
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
      llam_draw <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      
      j <- 1
      for (t in chunks[[thread]]) {
        
        plam <- gp_lam(object$llam_samples[t, ], object$x , x_new, g = object$g[t],
                       theta = object$theta_lam[t, ], mu_0 = settings$mean_lam, 
                       tau2_0 = object$tau2_lam[t], sep = settings$sep, v = settings$v)

        out$mean_lam[, j] <- plam$mean
        llam_draw[, j] <- mvtnorm::rmvnorm(1, mean =  plam$mean, 
                                           sigma = diag(plam$s2, nrow = length(plam$mean)))

        py <- gp_yw(object$y, object$x, x_new, A = object$A, lam = exp(object$llam_samples[t, ]),
                    lam_new = exp(llam_draw[, j]), theta = object$theta_y[t, ], mu = settings$mean_y, 
                    tau2 = object$tau2[t], sep = settings$sep, v = settings$v, lite = settings$lite)
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
    # Group elements out of the list
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

  return(object)
}

# no interval options available yet.
predict_nonvec_vdims <- function(object, x_new, settings, mapping = object$mappings){

  tic <- proc.time()[[3]]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)

  object$x_new <- x_new
  n_new <- nrow(object$x_new)

  if(!is.null(mapping$reps_vdims)){
    xv <- mapping$reps_vdims$X0
    Av <- mapping$reps_vdims$mult
    map <- mapping$reps_vdims$Z

    reps_vdims_new <- find_reps(x_new[, settings$vdims], 1:n_new)
    xv_new <- reps_vdims_new$X0
    Av_new <- reps_vdims_new$mult
    map_new <- reps_vdims_new$Z

    mean_np_lam <- s2_np_lam <- rep(NA, length(map_new))
    llam_n <- rep(NA, length(map))
  } else{
    xv <- object$x[ , settings$vdims, drop = FALSE]
    xv_new <- x_new[ , settings$vdims, drop = FALSE]
  }

  if (settings$cores == 1) { # run serial for loop

    mean_y <- matrix(nrow = n_new, ncol = object$nmcmc)
    yp_draw <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (settings$lite) {
      s2_sum_y <- rep(0, times = n_new)
      if (settings$return_all) s2_y <- matrix(nrow = n_new, ncol = object$nmcmc)
    } else sigma_sum_y <- matrix(0, nrow = n_new, ncol = n_new)
    
    mean_lam <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (settings$lite) {
      s2_sum_lam <- rep(0, times = n_new)
      if (settings$return_all) s2_lam <- matrix(nrow = n_new, ncol = object$nmcmc)
    } else sigma_sum_lam <- matrix(0, nrow = n_new, ncol = n_new)
    llam_draw <- matrix(nrow = n_new, ncol = object$nmcmc)

    for (t in 1:object$nmcmc) {

      # use tau2_lam draws
      plam <- gp_lam(object$llam_samples[t, ], xv , xv_new, g = object$g[t], theta = object$theta_lam[t, ], 
                     mu_0 = settings$mean_lam, tau2_0 = object$tau2_lam[t], 
                     sep = settings$sep, v = settings$v)

      if(!is.null(mapping$reps_vdims)){
        mean_np_lam[map_new] <- rep(plam$mean, Av_new)
        s2_np_lam[map_new] <- rep(plam$s2, Av_new)
        llam_draw_npv <- rmvnorm(1, mean =  plam$mean, sigma = plam$sigma, checkSymmetry = FALSE)
        llam_draw[map_new, t] <- rep(llam_draw_npv, Av_new)
        llam_n[map] <- rep(object$llam_samples[t, ], Av)
      }else{
        mean_np_lam <- plam$mean
        s2_np_lam <- plam$s2
        llam_draw[, t] <-  rmvnorm(1, mean =  plam$mean, sigma = plam$sigma, checkSymmetry = FALSE)
        llam_n <- object$llam_samples[t, ]
      }

      mean_lam[, t] <- mean_np_lam

      if (settings$lite) { # This should not be plam$s2 -> that's unordered
        s2_sum_lam <- s2_sum_lam + plam$s2
        if (settings$return_all) s2_lam[, t] <- plam$s2
      } else sigma_sum_lam <- sigma_sum_lam + plam$sigma
      
      py <- gp_yw(object$y, object$x, x_new, A = object$A, lam = exp(llam_n), 
                  lam_new = exp(llam_draw[, t]), theta = object$theta_y[t, ], mu = settings$mean_y, 
                  tau2 =object$tau2[t], sep = settings$sep, v = settings$v, lite = settings$lite)
      
      mean_y[, t] <- py$mean
      if (settings$lite) {
        s2_sum_y <- s2_sum_y + py$s2
        if (settings$return_all) s2_y[, t] <- py$s2 # storing individual variances for each posterior sample
      } else sigma_sum_y <- sigma_sum_y + py$sigma # storing main covariance matrix
      
      # yp_draw[, t] <- rmvnorm(1, mean = mean_y[, t], sigma = diag (py$s2))
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

      llam_draw <- yp_draw <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      j <- 1
      for (t in chunks[[thread]]) {
        
        # use tau2_lam draws
        plam <- gp_lam(object$llam_samples[t, ], xv , xv_new, g = object$g[t], theta = object$theta_lam[t, ], 
                       mu_0 = settings$mean_lam, tau2_0 = object$tau2_lam[t], 
                       sep = settings$sep, v = settings$v)
        
        if(!is.null(mapping$reps_vdims)){
          mean_np_lam[map_new] <- rep(plam$mean, Av_new)
          s2_np_lam[map_new] <- rep(plam$s2, Av_new)
          llam_draw_npv <- rmvnorm(1, mean =  plam$mean, sigma = plam$sigma, checkSymmetry = FALSE)
          llam_draw[map_new, j] <- rep(llam_draw_npv, Av_new)
          llam_n[map] <- rep(object$llam_samples[t, ], Av)
        }else{
          mean_np_lam <- plam$mean
          s2_np_lam <- plam$s2
          llam_draw[, j] <-  rmvnorm(1, mean =  plam$mean, sigma = plam$sigma, checkSymmetry = FALSE)
          llam_n <- object$llam_samples[t, ]
        }
        
        mean_lam[, j] <- mean_np_lam
        
        if (settings$lite) { # This should not be plam$s2 -> that's unordered
          s2_sum_lam <- s2_sum_lam + plam$s2
          if (settings$return_all) s2_lam[, t] <- plam$s2
        } else sigma_sum_lam <- sigma_sum_lam + plam$sigma
        
        py <- gp_yw(object$y, object$x, x_new, A = object$A, lam = exp(llam_n), 
                    lam_new = exp(llam_draw[, j]), theta = object$theta_y[t, ], mu = settings$mean_y, 
                    tau2 =object$tau2[t], sep = settings$sep, v = settings$v, lite = settings$lite)
        
        mean_y[, j] <- py$mean
        if (settings$lite) {
          s2_sum_y <- s2_sum_y + py$s2
          if (settings$return_all) s2_y[, j] <- py$s2 # storing individual variances for each posterior sample
        } else sigma_sum_y <- sigma_sum_y + py$sigma # storing main covariance matrix
        
        # yp_draw[, j] <- rmvnorm(1, mean = py$mean, sigma = diag (py$s2))
        j <- j + 1
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

predict_hom <- function(object, x_new, settings){

  tic <- proc.time()[[3]]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)

  object$x_new <- x_new
  n_new <- nrow(object$x_new)

  if (settings$cores == 1) { # run serial for loop

    mean_y <- matrix(nrow = n_new, ncol = object$nmcmc)

    if (settings$lite) {
      if(settings$interval == "pi" ||settings$interval == "both")
        s2_sum_y <- rep(0, times = n_new)
      if(settings$interval == "ci" ||settings$interval == "both")
        s2_sum_y_ci <- rep(0, times = n_new)
      if (settings$return_all){
        if(settings$interval == "pi" ||settings$interval == "both")
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

      py <- gp_yw_hom(object$y, object$x, x_new, A = object$A, nugs = object$g_y[t],
                      theta = object$theta_y[t, ], mu = settings$mean_y, tau2 = object$tau2[t],
                      sep = settings$sep, v = settings$v, lite = settings$lite)

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

      yp_draw <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      j <- 1
      for (t in chunks[[thread]]) {
        
        py <- gp_yw_hom(object$y, object$x, x_new, A = object$A, nugs = object$g_y[t],
                        theta = object$theta_y[t, ], mu = settings$mean_y, tau2 = object$tau2[t],
                        sep = settings$sep, v = settings$v, lite = settings$lite)
        
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

    # Group elements out of the list
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
  object$mean_y <- rowMeans(mean_y)

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

# Infer Yp at new locations xp

gp_lam <- function(llam, x, x_new, g, theta, mu_0, tau2_0, a0, b0, sep, v, lite = FALSE){
  
  if(v == 999) {
    if(sep) {
      K <- Exp2Sep(x, x, 1, theta, g)
      Kcross <- Exp2Sep(x_new, x, 1, theta, 0)
      if(!lite) K_new <- Exp2Sep(x_new, x_new, 1, theta, g)
    }
    else{
      K <- Exp2(sq_dist(x), 1, theta, g)
      Kcross <- Exp2(sq_dist(x_new, x), 1, theta, 0)
      if(!lite) K_new <- Exp2(sq_dist(x_new), 1, theta, g)
    }
  }else if(v > 1000){
    if(sep) {
      K <- MaternProdSep(x, x, 1, theta, g, (v - 1000))
      Kcross <- MaternProdSep(x_new, x, 1, theta, 0, (v - 1000))
      if(!lite) K_new <- MaternProdSep(x_new, x_new, 1, theta, g, (v - 1000))
    }
    else{
      K <- Matern(sq_dist(x), 1, theta, g, (v - 1000))
      Kcross <- Matern(sq_dist(x_new, x), 1, theta, 0, (v - 1000))
      if(!lite) K_new <- Matern(sq_dist(x_new), 1, theta, g, (v - 1000))
    }
  }else{
    if(sep) {
      K <- MaternSep(x, x, 1, theta, g, v)
      Kcross <- MaternSep(x_new, x, 1, theta, 0, v)
      if(!lite) K_new <- MaternSep(x_new, x_new, 1, theta, g, v)
    }
    else{
      K <- Matern(sq_dist(x), 1, theta, g, v)
      Kcross <- Matern(sq_dist(x_new, x), 1, theta, 0, v)
      if(!lite) K_new <- Matern(sq_dist(x_new), 1, theta, g, v)
    }
  }

  Kinv <- invdet(K)$Mi
  mean <- mu_0 + Kcross %*% Kinv %*% (llam - mu_0)

  K_cross_term <- Kcross %*% Kinv %*% t(Kcross)
  
  s2 <- tau2_0 * (1 + g - diag(K_cross_term))
  s2_ci <- tau2_0 * (1 - diag(K_cross_term))
  
  if(!lite){
    sigma <- sigma_ci <- tau2_0 * (K_new - K_cross_term)
    diag(sigma_ci) <- diag(sigma_ci) - tau2_0 * g
  } 
  else {
    sigma <- sigma_ci <- NULL
  }
  return(list(mean = mean, s2 = s2, s2_ci = s2_ci, sigma = sigma, sigma_ci = sigma_ci))

}

# If A is not specified - it is the no reps case.
gp_yw <- function(y, x, x_new, A = NULL, lam, lam_new, theta, mu, tau2, a, b, sep, v, lite = FALSE){
  
  if(!is.null(A)) lam <- lam/A

  if(v == 999) {
    if(sep) {
      K <- Exp2SepVec(x, x, 1, theta, lam)
      Kcross <- Exp2Sep(x_new, x, 1, theta, 0)
      if(!lite) K_new <- Exp2SepVec(x_new, x_new, 1, theta, lam_new)
    }
    else{
      K <- Exp2vec(sq_dist(x), 1, theta, lam)
      Kcross <- Exp2(sq_dist(x_new, x), 1, theta, 0)
      if(!lite) K_new <- Exp2vec(sq_dist(x_new), 1, theta, lam_new)
    }
  }else if(v > 1000){
    if(sep){
      K <- MaternProdSepVec(x, x, 1, theta, lam, (v- 1000))
      Kcross <- MaternProdSep(x_new, x, 1, theta, 0, (v- 1000))
      if(!lite) K_new <- MaternProdSepVec(x_new, x_new, 1, theta, lam_new, (v- 1000))
    }
    else{
      K <- MaternVec(sq_dist(x), 1, theta, lam, (v- 1000))
      Kcross <- Matern(sq_dist(x_new, x), 1, theta, 0, (v- 1000))
      if(!lite) K_new <- MaternVec(sq_dist(x_new), 1, theta, lam_new, (v- 1000))
    }
  }else{
    if(sep) {
      K <- MaternSepVec(x, x, 1, theta, lam, v)
      Kcross <- MaternSep(x_new, x, 1, theta, 0, v)
      if(!lite) K_new <- MaternSepVec(x_new, x_new, 1, theta, lam_new, v)
    }
    else{
      K <- MaternVec(sq_dist(x), 1, theta, lam, v)
      Kcross <- Matern(sq_dist(x_new, x), 1, theta, 0, v)
      if(!lite) K_new <- MaternVec(sq_dist(x_new), 1, theta, lam_new, v)
    }
  }

  Kinv <- invdet(K)$Mi
  mean <- mu + Kcross %*% Kinv %*% (y - mu)
  
  K_cross_term <- Kcross %*% Kinv %*% t(Kcross)
  s2 <- tau2 * (1 + lam_new - diag(K_cross_term))
  s2_ci <- tau2 * (1 - diag(K_cross_term))

  if(!lite){
    sigma <- sigma_ci <- tau2 * (K_new - K_cross_term)
    diag(sigma_ci) <- diag(sigma_ci) - tau2 * lam_new
  } 
  else sigma <- sigma_ci <- NULL
  

  # yp_draw <- drop(rmvnorm(1, mean = mean, sigma = diag(s2, nrow = length(mean))))

  return(list(mean = mean, s2 = s2, s2_ci = s2_ci, sigma = sigma, sigma_ci = sigma_ci))

}

gp_yw_hom <- function(y, x, x_new, A = NULL, nugs, theta, mu, tau2, a, b, sep, v, lite = FALSE){
  
  if(!is.null(A)) g_vec <- nugs/A
  else g_vec <- rep(nugs, nrow(x))
  g_new <- rep(nugs, nrow(x_new))

  if(v == 999) {
    if(sep) {
      K <- Exp2SepVec(x, x, 1, theta, g_vec)
      Kcross <- Exp2Sep(x_new, x, 1, theta, 0)
      if(!lite) K_new <- Exp2SepVec(x_new, x_new, 1, theta, g_new)
    }
    else{
      K <- Exp2vec(sq_dist(x), 1, theta, g_vec)
      Kcross <- Exp2(sq_dist(x_new, x), 1, theta, 0)
      if(!lite) K_new <- Exp2vec(sq_dist(x_new), 1, theta, g_new)
    }
  }else if(v > 1000){
    if(sep) {
      K <- MaternProdSepVec(x, x, 1, theta, g_vec, v - 1000)
      Kcross <- MaternProdSep(x_new, x, 1, theta, 0, v - 1000)
      if(!lite) K_new <- MaternProdSepVec(x_new, x_new, 1, theta, g_new, v - 1000)
    }
    else{
      K <- MaternVec(sq_dist(x), 1, theta, g_vec, v - 1000)
      Kcross <- Matern(sq_dist(x_new, x), 1, theta, 0, v - 1000)
      if(!lite) K_new <- MaternVec(sq_dist(x_new), 1, theta, g_new, v - 1000)
    }
  }else{
    if(sep) {
      K <- MaternSepVec(x, x, 1, theta, g_vec, v)
      Kcross <- MaternSep(x_new, x, 1, theta, 0, v)
      if(!lite) K_new <- MaternSepVec(x_new, x_new, 1, theta, g_new, v)
    }
    else{
      K <- MaternVec(sq_dist(x), 1, theta, g_vec, v)
      Kcross <- Matern(sq_dist(x_new, x), 1, theta, 0, v)
      if(!lite) K_new <- MaternVec(sq_dist(x_new), 1, theta, g_new, v)
    }
  }

  Kinv <- invdet(K)$Mi
  mean <- mu + Kcross %*% Kinv %*% (y - mu)

  K_cross_term <- Kcross %*% Kinv %*% t(Kcross)
  s2 <- tau2 * (1 + g_new - diag(K_cross_term))
  s2_ci <- tau2 * (1 - diag(K_cross_term))

  if(!lite){
    sigma <- sigma_ci <- tau2 * (K_new - K_cross_term)
    diag(sigma_ci) <- diag(sigma_ci) - tau2 *g_new 
  } 
  else sigma <- sigma_ci <- NULL

  # yp_draw <- drop(rmvnorm(1, mean = mean, sigma = diag(s2, nrow = length(mean))))

  return(list(mean = mean, s2 = s2, s2_ci = s2_ci, sigma = sigma, sigma_ci = sigma_ci))

}
