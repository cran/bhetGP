# ----- scaled Vecchia approximation (Katzfuss, Guinness, Lawrence)  ----------
# This code was appropriated from Katzfuss, and is basically a wrapper around 
# GpGp and GPVeccia packages, and is used to set good initial values for the 
# variance process, lambda. Please navigate to their git repo: 
# \url{https://github.com/katzfuss-group/scaledVecchia/blob/master/vecchia_scaled.R}
# for the script/details about the functions and examples

fit_scaled <- function(y, inputs, ms = 30, nu = 4.5, nug = NULL) {

  ## dimensions
  n <- nrow(inputs)
  d <- ncol(inputs)
  
  X <- as.matrix(sample(c(-1,1), n, replace=TRUE))
  beta <- mean(y)
  y <- y - beta
  cur.var <- summary(stats::lm(y ~ X-1))$sigma^2
  
  ## default range parameters
  # stop("here")
  input.ranges <- apply(inputs, 2, function(x) diff(range(x)))
  cur.ranges = 0.2*input.ranges
  active = rep(TRUE, d)

  # params
  if(is.null(nug)){
    fix.nug=FALSE; nug=0.01*var(y)
  } else fix.nug=TRUE
  
  # nug <- 0.1 * var(y)
  covfun=paste0("matern",nu*10,"_scaledim")
  
  n.est = min(5000, nrow(inputs))
  
  ## only use subsample for estimation.
  if(n.est < n){
    ind.est <- 1:n.est ## sample(1:n, n.est)
    y.full <- y; inputs.full <- inputs; X.full <- X
    y <- y[ind.est]; inputs <- inputs[ind.est, ,drop=FALSE]; X <- X[ind.est, ,drop=FALSE]
  }
  
  ### increase maxit until convergence
  conv = FALSE; maxit = 2
  while(conv == FALSE & maxit <= 30){
    
    ## check for inactive input dims (large range params)
    active = (cur.ranges < input.ranges * Inf)
    if(sum(active,na.rm=TRUE) == 0) stop('all inputs inactive. increase select?')
    
    ## specify how to scale input dimensions
    cur.ranges[!active] <- Inf
    scales <- 1/cur.ranges
    
    # stop("check here")
    ## order and condition based on current params
    # ord <- GPvecchia::order_maxmin_exact(t(t(inputs)*scales))
    ord <- 1:nrow(inputs)
    inputs.ord <- inputs[ord, ,drop=FALSE]
    y.ord <- y[ord]
    X.ord <- X[ord, ,drop=FALSE]
    NNarray <- GpGp::find_ordered_nn(t(t(inputs.ord)*scales), ms)
    
    ## starting and fixed parameters
    cur.parms <- c(cur.var, cur.ranges[active], nug)
    fixed=NULL      
    if(fix.nug) fixed=c(fixed,length(cur.parms))
    
    ## fisher scoring
    fit <- GpGp::fit_model(y.ord, inputs.ord[ , active, drop=FALSE], X.ord,
                           NNarray = NNarray, m_seq = ms,
                           start_parms = cur.parms, max_iter = maxit,
                           covfun_name = covfun, silent=TRUE,
                           reorder=FALSE, fixed_parms=fixed)
    
    cur.var <- fit$covparms[1]
    cur.ranges[active] <- fit$covparms[1+(1:sum(active))]
    nug <- fit$covparms[-(1:(1+sum(active)))]
    conv <- fit$conv
    maxit <- maxit * 2
  }

  fit$covparms = c(cur.var, cur.ranges, nug)
  if(n.est < n){
    fit$y <- y.full
    fit$locs <- inputs.full
    fit$X <- X.full
  } else {
    fit$locs <- inputs.ord
  }
  
  fit$betahat <- beta
  fit$y <- fit$y + beta
  fit$X <- as.matrix(rep(1,n))
  
  return(fit)
  
}

predictions_scaled <- function(fit, locs_pred, m = 100){
  
  y_obs <- fit$y
  locs_obs <- fit$locs
  X_obs <- fit$X
  n_pred <- nrow(locs_pred)
  
  X_pred <- as.matrix(rep(1, n_pred))
  scales <- 1/apply(locs_obs, 2, function(x) diff(range(x)))
  
  y  = y_obs - X_obs %*% fit$betahat
  
  # find the NNs 
  m <- min(m, nrow(locs_obs))
  NNarray <- FNN::get.knnx(t(t(locs_obs)*scales), t(t(locs_pred)*scales), m)$nn.index
  
  means <- numeric(length = n_pred)
  for(i in 1:n_pred){
    NN <- NNarray[i,]
    K <- get(fit$covfun_name)(fit$covparms, rbind(locs_obs[NN, , drop = FALSE], locs_pred[i, , drop = FALSE]))
    cl <- t(chol(K))
    means[i] <- cl[m+1,1:m] %*% forwardsolve(cl[1:m,1:m], y[NN])
  }
  
  means <- means + c(X_pred %*% fit$betahat)
  return(means)
}

