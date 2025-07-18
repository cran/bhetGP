#' @name plot
#' @title Plots object from \code{bhetGP} package
#' 
#' @description Acts on a \code{bhetgp}, \code{bhetgp_vec}, \code{bhomgp} or,
#'     \code{bhomgp_vec} object.  Generates trace plots for log likelihood of 
#'     mean and noise process, length scales of corresponding processes,
#'     scale parameters and the nuggets.
#'     Generates plots of hidden layers for one-dimensional inputs. Generates
#'     plots of the posterior mean and estimated 90\% prediction intervals for 
#'     one-dimensional inputs; generates heat maps of the posterior mean and 
#'     point-wise variance for two-dimensional inputs.
#'     
#' @details Trace plots are useful in assessing burn-in.  If there are too
#'     many hyperparameters to plot them all, then it is most useful to 
#'     visualize the log likelihood (e.g., \code{plot(fit$ll, type = "l")}).
#'
#' @param x object of class \code{bhetgp}, \code{bhetgp_vec}, \code{bhomgp}, 
#'        or \code{bhomgp_vec}
#' @param trace logical indicating whether to generate trace plots (default is
#'        TRUE if the object has not been through \code{predict})
#' @param predict logical indicating whether to generate posterior predictive 
#'        plot (default is TRUE if the object has been through \code{predict})
#' @param verb logical indicating whether to print plot.
#' @param ... N/A
#' 
#' @return ...N/A
#' 
#' @rdname plot
NULL

#' @rdname plot
#' @export
plot.bhetgp <- function(x, trace = NULL, predict = NULL, verb = TRUE,...) {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if(verb){
    if (is.null(x$mean)) {
      if (is.null(trace)) trace <- TRUE
      if (is.null(predict)) predict <- FALSE
    } else {
      if (is.null(trace)) trace <- FALSE
      if (is.null(predict)) predict <- TRUE
    }
    
    if(trace){
      if(!x$settings$sep){
        par(mfrow = c(2, 3), mar = c(5, 4, 2, 2))
        plot(x$theta_lam, type = 'l', lwd = 1, main = "Trace Plot (theta_lam)")
        plot(x$theta_y, type = 'l', lwd = 1, main = "Trace Plot (theta_y)")
        plot(x$tau2_lam, type = 'l', lwd = 1, main = "Trace Plot (tau2 lam)")
        plot(x$tau2, type = 'l', lwd = 1, main = "Trace Plot (tau2)")
        plot(x$llik_lam, type = 'l', lwd = 1, main = "Trace Plot (llik lam)")
        plot(x$llik_y, type = 'l', lwd = 1, main = "Trace Plot (llik y)")
        
        if(ncol(x$x) == 1) {
          par(mfrow = c(1, 1), mar = c(5, 4, 2, 2))
          matplot(t(x$llam_samples), lty = 1, type = 'l')
        }
      }else{
        par(mfrow = c(2, ncol(x$theta_y)), mar = c(5, 4, 2, 2))
        for (i in 1:ncol(x$theta_y))
          plot(x$theta_y[, i], type = 'l', ylab = 'theta_y', xlab = 'Iteration',
               main = paste0('Trace Plot of theta_y[', i, ']'))
        for (i in 1:ncol(x$theta_lam))
          plot(x$theta_lam[, i], type = 'l', ylab = 'theta_lam', xlab = 'Iteration',
               main = paste0('Trace Plot of theta_lam[', i, ']'))
        
        par(mfrow = c(1, 4), mar = c(5, 4, 2, 2))
        plot(x$tau2_lam, type = 'l', ylab = 'tau2 lam', xlab = 'Iteration',
             main = paste0('Trace Plot of tau2 lam'))
        plot(x$tau2, type = 'l', ylab = 'tau2', xlab = 'Iteration',
             main = paste0('Trace Plot of tau2'))
        plot(x$llik_lam, type = 'l', ylab = 'llik lam', xlab = 'Iteration',
             main = paste0('Trace Plot of llik lam'))
        plot(x$llik_y, type = 'l', ylab = 'llik y', xlab = 'Iteration',
             main = paste0('Trace Plot of llik y'))
      }
    }
    
    if (predict) {
      if (ncol(x$x) == 1) {
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
        if (is.null(x$Sigma) && is.null(x$Sigma_ci)) {
          s2 <- ifel(is.null(x$s2_y), x$s2_y_ci, x$s2_y)
          q1 <- x$mean + qnorm(0.05, 0, sqrt(s2))
          q3 <- x$mean + qnorm(0.95, 0, sqrt(s2))
        } else {
          # Sigma_smooth <- x$Sigma - diag(exp(x$mean_lnugs) * mean(x$tau2), nrow(x$x_new))
          # y_samples <- t(mvtnorm::rmvnorm(50, x$mean, x$Sigma_ci))
          # y_samples <- x$mean
          sig <- ifel(is.null(x$Sigma), x$Sigma_ci, x$Sigma)
          q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(sig)))
          q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(sig)))
        }
        o <- order(x$x_new)
        plot(x$x_new[o], x$mean[o], type = 'l', xlab = 'X', ylab = 'Y', 
             ylim = c(min(q1), max(q3)),
             col = 'blue')
        # if (!is.null(x$Sigma_ci)) {
        #   matlines(x$x_new[o], y_samples[o,], col = 'lightblue', lty = 1)
        #   lines(x$x_new[o], x$mean[o], col = 'blue')
        # }
        lines(x$x_new[o], q1[o], col = 'blue', lty = 2)
        lines(x$x_new[o], q3[o], col = 'blue', lty = 2)
        points(rep(x$x, x$A), unlist(x$Ylist), pch = 20, col = "gray")
        points(x$x, x$y, pch = 20)
      } else if (ncol(x$x) == 2) {
        if (!requireNamespace("interp", quietly = TRUE)) {
          stop("Package \"interp\" needed for this plot. Please install it.",
               call. = FALSE)
        }
        cols <- heat.colors(128)
        i1 <- interp::interp(x$x_new[, 1], x$x_new[, 2], x$mean)
        if (is.null(x$Sigma) && is.null(x$Sigma_ci)) {
          s2 <- ifel(is.null(x$s2_y), x$s2_y_ci, x$s2_y)
          i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(s2))
        } else{
          sig <- ifel(is.null(x$Sigma), x$Sigma_ci, x$Sigma)
          i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(diag(sig)))
        } 
        
        par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
        image(i1, col = cols, main = 'Posterior Mean', xlab = 'X1', ylab = 'X2')
        contour(i1, add = TRUE)
        points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
        image(i2, col = cols, main = 'Posterior Variance', xlab = 'X1', 
              ylab = 'X2')
        contour(i2, add = TRUE)
        points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      } else cat('Dimension of X too large for default plotting')
    }
  }else warning("set verb = TRUE to view plots")
  
}

#' @rdname plot
#' @export
plot.bhomgp <- function(x, trace = NULL, predict = NULL, verb = TRUE,...) {
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if(verb){
    if (is.null(x$mean)) {
      if (is.null(trace)) trace <- TRUE
      if (is.null(predict)) predict <- FALSE
    } else {
      if (is.null(trace)) trace <- FALSE
      if (is.null(predict)) predict <- TRUE
    }
    
    if(trace){
      if(!x$settings$sep){
        par(mfrow = c(1, 4), mar = c(5, 4, 2, 2))
        plot(x$theta_y, type = 'l', lwd = 1, main = "Trace Plot (theta_y)")
        plot(x$g_y, type = 'l', lwd = 1, main = "Trace Plot (g_y)")
        plot(x$tau2, type = 'l', lwd = 1, main = "Trace Plot (tau2)")
        plot(x$llik, type = 'l', lwd = 1, main = "Trace Plot (llik)")
      }else{
        par(mfrow = c(1, ncol(x$theta_y)), mar = c(5, 4, 2, 2))
        for (i in 1:ncol(x$theta_y))
          plot(x$theta_y[, i], type = 'l', ylab = 'theta_y', xlab = 'Iteration',
               main = paste0('Trace Plot of theta_y[', i, ']'))
        par(mfrow = c(1, 4), mar = c(5, 4, 2, 2))
        plot(x$g_y, type = 'l', ylab = 'g', xlab = 'Iteration',
             main = paste0('Trace Plot of g'))
        plot(x$tau2, type = 'l', ylab = 'tau2', xlab = 'Iteration',
             main = paste0('Trace Plot of tau2'))
        plot(x$llik, type = 'l', ylab = 'llik', xlab = 'Iteration',
             main = paste0('Trace Plot of llik'))
      }
    }
    
    if (predict) {
      if (ncol(x$x) == 1) {
        par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
        if (is.null(x$Sigma) && is.null(x$Sigma_ci)) {
          s2 <- ifel(is.null(x$s2_y), x$s2_y_ci, x$s2_y)
          q1 <- x$mean + qnorm(0.05, 0, sqrt(s2))
          q3 <- x$mean + qnorm(0.95, 0, sqrt(s2))
        } else {
          # Sigma_smooth <- x$Sigma - diag(exp(x$mean_lnugs) * mean(x$tau2), nrow(x$x_new))
          # y_samples <- t(mvtnorm::rmvnorm(50, x$mean, x$Sigma_ci))
          # y_samples <- x$mean
          sig <- ifel(is.null(x$Sigma), x$Sigma_ci, x$Sigma)
          q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(sig)))
          q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(sig)))
        }
        o <- order(x$x_new)
        plot(x$x_new[o], x$mean[o], type = 'l', xlab = 'X', ylab = 'Y', 
             ylim = c(min(q1), max(q3)),
             col = 'blue')
        # if (!is.null(x$Sigma_ci)) {
        #   matlines(x$x_new[o], y_samples[o,], col = 'lightblue', lty = 1)
        #   lines(x$x_new[o], x$mean[o], col = 'blue')
        # }
        lines(x$x_new[o], q1[o], col = 'blue', lty = 2)
        lines(x$x_new[o], q3[o], col = 'blue', lty = 2)
        points(rep(x$x, x$A), unlist(x$Ylist), pch = 20, col = "gray")
        points(x$x, x$y, pch = 20)
      } else if (ncol(x$x) == 2) {
        if (!requireNamespace("interp", quietly = TRUE)) {
          stop("Package \"interp\" needed for this plot. Please install it.",
               call. = FALSE)
        }
        cols <- heat.colors(128)
        i1 <- interp::interp(x$x_new[, 1], x$x_new[, 2], x$mean)
        if (is.null(x$Sigma) && is.null(x$Sigma_ci)) {
          s2 <- ifel(is.null(x$s2_y), x$s2_y_ci, x$s2_y)
          i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(s2))
        } else{
          sig <- ifel(is.null(x$Sigma), x$Sigma_ci, x$Sigma)
          i2 <- interp::interp(x$x_new[, 1], x$x_new[, 2], sqrt(diag(sig)))
        } 
        
        par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
        image(i1, col = cols, main = 'Posterior Mean', xlab = 'X1', ylab = 'X2')
        contour(i1, add = TRUE)
        points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
        image(i2, col = cols, main = 'Posterior Variance', xlab = 'X1', 
              ylab = 'X2')
        contour(i2, add = TRUE)
        points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      } else cat('Dimension of X too large for default plotting')
    }
  } else warning("set verb = TRUE to view plot")
  
}

#' @rdname plot
#' @export
plot.bhetgp_vec <- plot.bhetgp

#' @rdname plot
#' @export
plot.bhomgp_vec <- plot.bhomgp
