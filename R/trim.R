#' @title Trim/Thin MCMC iterations
#' @description Acts on a \code{bhetgp}, \code{bhetgp_vec}, \code{bhomgp}, or 
#'    \code{bhomgp_vec} object.
#'    Removes the specified number of MCMC iterations (starting at the first 
#'    iteration).  After these samples are removed, the remaining samples are
#'    optionally thinned.
#' 
#' @details The resulting object will have \code{nmcmc} equal to the previous 
#'     \code{nmcmc} minus \code{burn} divided by \code{thin}. Removing burn-ins 
#'     are necessary following convergence. Thinning is recommended as it can 
#'     eliminate highly correlated consecutive samples. Additionally, the size of
#'     the object reduces and ensures faster prediction.
#'     
#' @param object object from \code{bhetGP}, or \code{bhomGP}
#' @param burn integer specifying number of iterations to cut off as burn-in
#' @param thin integer specifying amount of thinning (\code{thin = 1} keeps all 
#'        iterations, \code{thin = 2} keeps every other iteration, 
#'        \code{thin = 10} keeps every tenth iteration, etc.)
#'        
#' @return object of the same class with the selected iterations removed
#' 
#' @rdname trim
#' @export
trim <- function(object, burn, thin){
  
  tic <- proc.time()[3]
  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')
  nmcmc <- object$nmcmc
  its <- seq((burn + 1), nmcmc, by= thin)
  object$nmcmc <- length(its)
  
  if(inherits(object, "bhetgp") | inherits(object, "bhetgp_vec")){
    object$llam_samples <- object$llam_samples[its, , drop = FALSE]
    object$theta_lam <- object$theta_lam [its, , drop = FALSE]
    object$theta_y <- object$theta_y[its, , drop = FALSE]
    object$tau2 <- object$tau2[its, drop = FALSE]
    object$llik_lam <- object$llik_lam[its, drop = FALSE]
    object$llik_y <- object$llik_y[its, drop = FALSE]
    object$g <- object$g[its, drop = FALSE]
    if(!is.null(object$tau2_lam)) object$tau2_lam <- object$tau2_lam[its, drop = FALSE]
  }
  else{
    object$g_y <- object$g_y [its, , drop = FALSE]
    object$theta_y <- object$theta_y[its, , drop = FALSE]
    object$tau2 <- object$tau2[its, drop = FALSE]
    object$llik_y <- object$llik_y[its, drop = FALSE]
  }
  
  toc <- proc.time()[3]
  object$time <- object$time + unname(toc - tic)
  
  return(object)
}

