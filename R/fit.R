# --------- Functions ----------------------------------------------------------
# External 
# bhetGP: fits hetGP with vecchia/no vec, sep/iso and reps/no reps w/ vdims
# bhomGP: fits homGP with vec/no vec, sep/iso and reps/no reps
# ------------------------------------------------------------------------------

# bhetGP  ---------------------------------------------------------------
#' @title MCMC sampling for Heteroskedastic GP
#' @description Conducts MCMC sampling of hyperparameters and latent noise 
#'     process \code{llam} for a hetGP.  Separate length scale 
#'     parameters \code{theta_lam} and \code{theta_y} govern the correlation 
#'     strength of the hidden layer and outer layer respectively.  
#'     \code{lam} layer may have a non-zero nugget \code{g} which governs 
#'     noise for the latent noise layer. \code{tau2_y} and \code{tau2_lam}
#'     control the amplitude of the mean and noise process respectively.
#'    In Matern covariance, \code{v} governs smoothness.
#'     
#' @details Maps inputs \code{x} to mean response \code{y} and noise levels
#'     \code{llam}. Conducts sampling of the latent noise process using Elliptical 
#'     Slice sampling.  Utilizes Metropolis Hastings sampling of the length 
#'     scale and nugget parameters with proposals and priors controlled by 
#'     \code{priors}. \code{g} for the noise process is set to a specific 
#'     value, and by default, is not estimated.  When \code{vecchia = TRUE}, 
#'     all calculations leverage the Vecchia approximation with 
#'     specified conditioning set size \code{m}. \code{tau2_y} is always 
#'     inferred from likelihood; \code{tau2_lam} is inferred by default but 
#'     may be pre-specified and fixed.
#'     
#'     NOTE on OpenMP: The Vecchia implementation relies on OpenMP parallelization
#'     for efficient computation.  This function will produce a warning message 
#'     if the package was installed without OpenMP (this is the default for 
#'     CRAN packages installed on Apple machines).  To set up OpenMP 
#'     parallelization, download the package source code and install 
#'     using the gcc/g++ compiler.    
#'     
#'     Proposals for \code{g} and \code{theta} follow a uniform sliding window 
#'     scheme, e.g. 
#'     
#'     \code{theta_star <- runif(1, l * theta_t / u, u * theta_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{priors}.
#'     To adjust these, set \code{priors = list(l = new_l, u = new_u)}.    
#'     Priors on \code{g}, \code{theta_y}, and \code{theta_lam} follow Gamma 
#'     distributions with shape parameters (\code{alpha}) and rate parameters 
#'     (\code{beta}) controlled within the \code{priors} list object.  
#'     Defaults are
#'     \itemize{
#'         \item \code{priors$alpha$theta_lam <- 1.5}
#'         \item \code{priors$beta$theta_lam <- 4}
#'         \item \code{priors$alpha$theta_y <- 1.5}
#'         \item \code{priors$beta$theta_y <- 4}
#'         \item \code{priors$alpha$g <- 1.5}
#'         \item \code{priors$beta$g <- 4}
#'     }
#'     
#'     \code{tau2_y} and \code{tau2_lam} are not sampled; rather directly inferred
#'     under conjugate Inverse Gamma prior with shape (\code{alpha}) and scale 
#'     parameters (\code{beta}) within the \code{priors} list object
#'     \itemize{       
#'         \item \code{priors$a$tau2_y <- 10}
#'         \item \code{priors$a$tau2_y <- 4}
#'         \item \code{priors$a$tau2_lam <- 10}
#'         \item \code{priors$a$tau2_lam <- 4}
#'     }
#'     These priors are designed for \code{x} scaled to 
#'     [0, 1] and \code{y} having mean \code{mean_y}.  These may be 
#'     adjusted using the \code{priors} input.
#'     
#'     Initial values for \code{theta_y}, \code{theta_lam}, \code{llam} may be
#'     specified by the user. If no initial values are specified, \code{stratergy}
#'     will determine the initialization method. \code{stratergy = "default"} 
#'     leverages mleHetGP for initial values of hyper-parameters if 
#'     \code{vecchia = FALSE} and Scaled-Vecchia with Stochastic Kriging (Sk-Vec)
#'     hybrid approach if \code{vecchia = TRUE}.
#'     
#'     For SK-Vec hybrid approach, scaled Vecchia code from 
#'     \url{https://github.com/katzfuss-group/scaledVecchia/blob/master/vecchia_scaled.R}
#'     is used to fit two GPs using the Vecchia approximation. The first for (x, y) pairs,
#'     which result in estimated residual sums of squares 
#'     based on predicted y values. Another GP on (x, s) to obtain
#'     latent noise estimates which are smoothed.  A script is leveraged
#'     internally within this package that fits this method. 
#'     
#'     Optionally, choose stratergy = "flat" which which will start at 
#'     uninformative initial values; \code{llam} = log(var(y) * 0.1) or
#'     specify initial values. 
#'     
#'     The output object of class \code{bhetgp} or \code{bhetgp_vec} is designed for 
#'     use with \code{trim}, \code{predict}, and \code{plot}.   
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param reps_list list object from hetGP::find_reps 
#' @param nmcmc number of MCMC iterations
#' @param sep logical indicating whether to fit isotropic GP (\code{sep = FALSE})
#'            or seperable GP (\code{sep = TRUE})
#' @param inits set initial values for hyparameters: \code{llam}, \code{theta_y}, 
#'              \code{theta_lam}, \code{g}, \code{mean_y}, \code{mean_lam},
#'              \code{scale_y}, \code{scale_lam}.
#'              Additionally, set initial conditions for tuning:
#'              \itemize{ 
#'                \item \code{theta_check}: logical; if \code{theta_check = TRUE},
#'                then ensures that theta_lam > theta_y i.e., decay of correlation
#'                for noise process is slower than mean process.
#'                \item \code{prof_ll_lam}: logical; if \code{prof_ll_lam = TRUE},
#'                infers tau2_lam i.e., scale parameter for latent noise process
#'                \item \code{noise}: logical; if \code{noise = TRUE}, infers nugget
#'                \code{g} throught M-H for latent noise process.
#'              }
#' @param priors hyperparameters for priors and proposals (see details)
#' @param reps logical; if \code{reps = TRUE} uses Woodbury inference adjusting for
#'             replication of design points and \code{reps = FALSE} does not
#'             use Woodbury inference
#' @param cov covariance kernel, either Matern, ARD Matern 
#'        or squared exponential (\code{"exp2"})
#' @param stratergy choose initialization stratergy; "default" uses hetGP for 
#'        \code{vecchia = FALSE} settings and sVecchia for \code{vecchia = TRUE}. 
#'        See details.
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @param vecchia logical indicating whether to use Vecchia approximation
#' @param m size of Vecchia conditioning sets (only used if 
#'        \code{vecchia = TRUE})
#' @param ordering optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x}, defaults to random, is applied to \code{x} (only used if 
#'        \code{vecchia = TRUE})
#' @param verb logical indicating whether to print progress
#' @param omp_cores logical; if \code{vecchia = TRUE}, user may specify the number of cores to
#'        use for OpenMP parallelization. Uses min(4, limit) where limit is max openMP 
#'        cores available on the machine.
#' @return a list of the S3 class \code{bhetgp} or \code{bhetgp_vec} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response mean at inputs (x)
#'   \item \code{Ylist}: list of all responses observed per location (x)
#'   \item \code{A}: number of replicates at each location
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{priors}: copy of proposal/priors
#'   \item \code{settings}: setting for predictions using \code{bhetgp} or \code{bhetgp_vec}
#'   object
#'   \item \code{theta_y}: vector of MCMC samples for \code{theta_y} (length
#'         scale of mean process)
#'   \item \code{theta_lam}: matrix of MCMC samples for \code{theta_lam} (length 
#'         scale of latent noise process)
#'   \item \code{llam_samples}: matrix of ESS samples for \code{log lambda} (latent 
#'   noise process samples)
#'   \item \code{g}: vector of MCMC samples for \code{g} if infered
#'   \item \code{tau2}: vector of MAP estimates for \code{tau2} (scale 
#'         parameter of mean process)
#'   \item \code{tau2_lam}: vector of MAP estimates for \code{tau2_lam} (scale 
#'         parameter of latent noise process)
#'   \item \code{llik_y}: vector of MVN log likelihood of Y for reach Gibbs iteration
#'   \item \code{llik_lam}: vector of MVN log likelihood of \code{llam} i.e.
#'              the latent noise process for reach Gibbs iteration
#'   \item \code{x_approx}: conditioning set, NN and ordering for \code{vecchia = TRUE}
#'   \item \code{m}: copy of size of conditioning set for \code{vecchia = TRUE}            
#'   \item \code{time}: computation time in seconds
#'   
#' }
#' 
#' @references
#' 
#' Binois, Mickael, Robert B. Gramacy, and Mike Ludkovski. "Practical heteroscedastic Gaussian process 
#' modeling for large simulation experiments." Journal of Computational and Graphical 
#' Statistics 27.4 (2018): 808-821.
#' 
#' Katzfuss, Matthias, Joseph Guinness, and Earl Lawrence. "Scaled Vecchia approximation for 
#' fast computer-model emulation." SIAM/ASA Journal on Uncertainty Quantification 10.2 (2022): 537-554.
#' 
#' Sauer, Annie Elizabeth. "Deep Gaussian process surrogates for computer experiments." (2023).
#' 
#' @examples 
#' 
#' # 1D function with 1D noise 
#' 
#' # Truth
#' fx <- function(x){
#' result <- (6 * x - 2)^2* sin(12 * x - 4)
#' }
#' 
#' # Noise
#' rx <- function(x){
#' result <- (1.1 + sin(2 * pi * x))^2
#' return(result)
#' }
#' 
#' # Training data
#' r <- 10 # replicates
#' xn <- seq(0, 1, length = 25)
#' x <- rep(xn, r)
#' 
#' rn <- drop(rx(x))
#' noise <- as.numeric(t(mvtnorm::rmvnorm(r, sigma = diag(rn, length(xn)))))
#' 
#' f <- fx(x) 
#' y <- f + noise
#' 
#' # Testing data
#' xx <- seq(0, 1, length = 100)
#' yy <- fx(xx)
#' 
#' #--------------------------------------------------------------------------- 
#' # Example 1: Full model, no Vecchia 
#' #---------------------------------------------------------------------------
#' 
#' # Fitting a bhetGP model using all the data
#' fit <- bhetGP(x, y, nmcmc = 100, verb = FALSE)
#' 
#' # Trimming the object to remove burn in and thin samples
#' fit <- trim(fit, 50, 10)
#' 
#' # Predition using the bhetGP object (indepedent predictions)
#' fit <- predict(fit, xx, cores = 2) 
#' 
#' # Visualizing the mean predictive surface. 
#' # Can run plot(fit, trace = TRUE) to view trace plots
#' plot(fit) 
#' 
#' 
#' #---------------------------------------------------------------------------
#' # Example 2: Vecchia approximated model
#' #---------------------------------------------------------------------------
#' 
#' # Fitting a bhetGP model with vecchia approximation. Two cores for OpenMP
#' fit <- bhetGP(x, y, nmcmc = 100, vecchia = TRUE, m = 5, omp_cores = 2, verb = FALSE)
#' 
#' # Trimming the object to remove burn in and thin samples
#' fit <- trim(fit, 50, 10)
#' 
#' # Predition using the bhetGP_vec object with joint predictions (lite = FALSE)
#' # Two cores for OpenMP, default setting (omp_cores = 2). No SNOW
#' fit <- predict(fit, xx, lite = FALSE, vecchia = TRUE) 
#'
#' # Visualizing the mean predictive surface
#' plot(fit)
#' 
#' \donttest{
#' #--------------------------------------------------------------------------- 
#' # Example 3: Vecchia inference, non-vecchia predictions
#' #---------------------------------------------------------------------------
#' 
#' # Fitting a bhetGP model with vecchia approximation. Two cores for OpenMP
#' fit <- bhetGP(x, y, nmcmc = 200, vecchia = TRUE, m = 5, omp_cores = 2)
#' 
#' # Trimming the object to remove burn in and thin samples
#' fit <- trim(fit, 100, 10)
#'
#' # Predition using the bhetGP object with joint predictions (lite = FALSE)
#' # Two cores for OpenMP which is default setting (omp_cores = 2)
#' # Two cores for SNOW (cores = 2)
#' fit <- predict(fit, xx, vecchia = FALSE, cores = 2, lite = FALSE)
#' 
#' # Visualizing the mean predictive surface
#' plot(fit)
#' }
#' 
#' @export

bhetGP <- function(x = NULL, y = NULL, reps_list = NULL, nmcmc = 1000, 
                   sep = TRUE, inits = NULL, # ls_check = FALSE,
                   priors = NULL, reps = TRUE,
                   cov = c("exp2", "matern", "ARD matern"), v = 2.5, 
                   stratergy = c("default", "flat"),
                   vecchia = FALSE, m = min(25, length(y) - 1), ordering = NULL, 
                   verb = TRUE, omp_cores = 4){
  
  tic <- proc.time()[[3]]
  cov <- match.arg(cov)
  stratergy <- match.arg(stratergy)
  
  D = ifelse(is.matrix(x), ncol(x), 1)
  
  if (cov == "exp2") v <- 999 # indicator
  if (!vecchia & length(y) > 300)
    message("'vecchia = TRUE' is recommended for faster computation.")
  if (nmcmc <= 0) stop("nmcmc must be at least 1")
  if(!is.null(reps_list)){x <- reps_list$X0; y <- reps_list$Z0}
  if (is.numeric(x)) x <- as.matrix(x)
  
  if(reps){ # woodbury version
    if(verb) print("obtaining replication")
    reps <- ifel(is.null(reps_list), hetGP::find_reps(x, y), reps_list)
    Yn <- reps$Z0
    Xn <- reps$X0
    # Ylist <- reps$Zlist # for vecchia
    YNs2 <- unlist(lapply(reps$Zlist, function(x) sum((x - mean(x))^2)))
    A <- reps$mult
  } else{ # no reps
    reps <- NULL # no mapping
    Yn <- y
    Xn <- x
    YNs2 <- NULL # indicator
    A <- NULL # indicator
  }

  #check matern 0.5 -> weird results
  if (m >= length(y)) stop("m must be less than the length of y")
  if (cov == "matern" | cov == "ARD matern")
    if(!(v %in% c(1.5, 2.5)))
      stop("v must be one of 1.5, or 2.5")
  if(cov == "matern")
    v <- 1000 + v # indicator + v
  if (!is.null(ordering)) {
    if (!vecchia) message("ordering only used when vecchia = TRUE")
    test <- check_ordering(ordering, nrow(Xn)) # returns NULL if all checks pass
  }
  
  # check prior settings
  priors <- check_settings(priors, gp_type = "hetgp", x = Xn, y = Yn)

  if(is.null(inits$mean_y)) inits$mean_y <- mean(y)
  initial <- list(theta_y = inits$theta_y_0, theta_lam = inits$theta_lam_0, 
                  g = inits$g_0, llam = inits$llam_0, 
                  mean_y = inits$mean_y , mean_lam = inits$mean_lam,
                  scale_y = 1, scale_lam = inits$scale_lam,
                  theta_check = inits$ls_check, prof_ll_lam = inits$prof_ll_lam, 
                  noise = inits$noise)

  initial <- check_initialization(reps = ifel(is.null(reps), hetGP::find_reps(x, y), reps), 
                                  initial, n = nrow(Xn), D = D, sep = sep, v = v, 
                                  vec = vecchia, verb = verb, stratergy = stratergy)

  settings <- list(v = v, sep = sep, mean_y = initial$mean_y, mean_lam = initial$mean_lam,
                   verb = verb)

  out <- list(x = Xn, y = Yn, Ylist = reps$Zlist, A = reps$mult, 
              nmcmc = nmcmc, priors = priors, settings = settings)

  if(!vecchia){ # non vecchia

    initial <- check_scale(x = Xn, v = v, vec = FALSE, priors, initial)
    calc <- precalc(Xn, out = Yn, outs2 = YNs2, A = A, v, vec = FALSE, priors, initial, 
                    latent = FALSE)
    initial$llik_y <- calc$llik
    initial$tau2_y <- calc$tau2
    
    test <- check_inputs(x, y) # returns NULL if all checks pass
    
    
    if (nmcmc == 1) {
      out <- c(out, initial)
      out$time <- proc.time()[[3]] - tic
      return(out)
    }
    
    if(verb) print("starting mcmc")
    # non vdims, YN = NULL will do no reps case
    # if(adj == T) obj <- gibbs_sep_N(YNs2 = YNs2, Yn = Yn, Xn = Xn, A= A, mcmc = nmcmc,
    #                                   initial = initial, priors = priors, v = v) # this is for ess fig illustration 
    # else{
    if(sep) obj <- gibbs_sep(YNs2 = YNs2, Yn = Yn, Xn = Xn, A= A, mcmc = nmcmc,
                                 initial = initial, priors = priors, v = v, verb = verb)
    else obj <- gibbs_iso(YNs2 = YNs2, Yn = Yn,  Xn = Xn, A= A, mcmc = nmcmc,
                              initial = initial, priors = priors, v = v, verb = verb)
  } else{ # vecchia
    out$m <- m # neighbours
    out$ordering <- ordering # order
    x_approx <- create_approx(Xn, m, ordering, cores = omp_cores) # create approximation

    if(is.null(reps))
      YNs2_ord <- NULL # If no reps, use Xn and Yn (no YN)
    else
      YNs2_ord <- YNs2[x_approx$ord]
      # YN_ord <- unlist(Ylist[x_approx$ord]) # order Ylists by NN and then unlist for YN

    Yn_ord <- Yn[x_approx$ord] # order Yn
    out$x_approx <- x_approx # store approximation

    initial <- check_scale(x = x_approx, v = v, vec = TRUE, priors, initial)
    calc <- precalc(x_approx, out = Yn_ord, outs2 = YNs2_ord, 
                    A, v, vec = TRUE, priors, initial, latent = FALSE)
    initial$llik_y <- calc$llik
    initial$tau2_y <- calc$tau2
    
    test <- check_inputs(x, y) # returns NULL if all checks pass
    
    if (nmcmc == 1) {
      out <- c(out, initial)
      out$time <- proc.time()[[3]] - tic
      return(out)
    }
    
    if(verb) print("starting mcmc")
    
    if(sep) obj <- gibbs_sep_vec(YNs2 = YNs2_ord, Yn = Yn_ord,  x_approx, A= A, mcmc = nmcmc,
                                 initial = initial, priors = priors, v = v, verb = verb)
    else obj <- gibbs_iso_vec(YNs2 = YNs2_ord, Yn = Yn_ord,  x_approx, A= A, mcmc = nmcmc,
                              initial = initial, priors = priors, v = v, verb = verb)
  }

  out <- c(out, obj)
  toc <- proc.time()[[3]]
  out$time <- unname(toc - tic)

  if (vecchia) class(out) <- "bhetgp_vec" else class(out) <- "bhetgp"

  return(out)
}

# bhomGP  ---------------------------------------------------------------
#' @title MCMC sampling for Homoskedastic GP
#' @description Conducts MCMC sampling of hyperparameters for a homGP.
#'     Separate length scale parameters \code{theta_y} 
#'     govern the correlation strength of the response. \code{g} governs 
#'     noise for the noise. \code{tau2_y} control the amplitude of the mean process.
#'    In Matern covariance, \code{v} governs smoothness.
#'     
#' @details Maps inputs \code{x} to mean response \code{y}. Utilizes Metropolis Hastings 
#'     sampling of the length scale and nugget parameters with proposals and priors 
#'     controlled by \code{priors}. \code{g} is estimated by default but may be specified and fixed.
#'     When \code{vecchia = TRUE}, all calculations leverage the Vecchia approximation with 
#'     specified conditioning set size \code{m}. \code{tau2_y} is inferred by default but 
#'     may be pre-specified and fixed.
#'     
#'     NOTE on OpenMP: The Vecchia implementation relies on OpenMP parallelization
#'     for efficient computation.  This function will produce a warning message 
#'     if the package was installed without OpenMP (this is the default for 
#'     CRAN packages installed on Apple machines).  To set up OpenMP 
#'     parallelization, download the package source code and install 
#'     using the gcc/g++ compiler.    
#'     
#'     Proposals for \code{g} and \code{theta} follow a uniform sliding window 
#'     scheme, e.g. 
#'     
#'     \code{theta_star <- runif(1, l * theta_t / u, u * theta_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{priors}.
#'     To adjust these, set \code{priors = list(l = new_l, u = new_u)}.    
#'     Priors on \code{g}, \code{theta_y} follow Gamma 
#'     distributions with shape parameters (\code{alpha}) and rate parameters 
#'     (\code{beta}) controlled within the \code{priors} list object.  
#'     Defaults are
#'     \itemize{
#'         \item \code{priors$alpha$theta_lam <- 1.5}
#'         \item \code{priors$beta$theta_lam <- 4}
#'         \item \code{priors$alpha$theta_y <- 1.5}
#'         \item \code{priors$beta$theta_y <- 4}
#'         \item \code{priors$alpha$g <- 1.5}
#'         \item \code{priors$beta$g <- 4}
#'     }
#'     
#'     \code{tau2_y} is not sampled; rather directly inferred
#'     under conjugate Inverse Gamma prior with shape (\code{alpha}) and scale 
#'     parameters (\code{beta}) within the \code{priors} list object
#'     \itemize{       
#'         \item \code{priors$a$tau2_y <- 10}
#'         \item \code{priors$a$tau2_y <- 4}
#'     }
#'     These priors are designed for \code{x} scaled to 
#'     [0, 1] and \code{y} having mean \code{mean_y}.  These may be 
#'     adjusted using the \code{priors} input.
#'     
#'     Initial values for \code{theta_y}, and \code{g} may be
#'     specified by the user. If no initial values are specified, \code{stratergy}
#'     will determine the initialization method. \code{stratergy = "default"} 
#'     leverages mleHomGP for initial values of hyper-parameters if 
#'     \code{vecchia = FALSE} and scaled vecchia approach if \code{vecchia = TRUE}.
#'     
#'     For scaled Vecchia code from 
#'     \url{https://github.com/katzfuss-group/scaledVecchia/blob/master/vecchia_scaled.R}
#'     is used to fit a vecchia approximated GP to (x, y). A script is leveraged
#'     internally within this package that fits this method. 
#'     
#'     Optionally, choose stratergy = "flat" which will start at 
#'     uninformative initial values or specify initial values. 
#'     
#'     The output object of class \code{bhomgp} or \code{bhomgp_vec} is designed for 
#'     use with \code{trim}, \code{predict}, and \code{plot}.   
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param reps_list list object from hetGP::find_reps 
#' @param nmcmc number of MCMC iterations
#' @param sep logical indicating whether to fit isotropic GP (\code{sep = FALSE})
#'            or seperable GP (\code{sep = TRUE})
#' @param inits set initial values for hyparameters: \code{theta_y}, \code{g}, \code{mean_y},
#'              \code{scale_y},
#'              Additionally, set initial conditions for tuning:
#'              \itemize{ 
#'                \item \code{prof_ll}: logical; if \code{prof_ll = TRUE},
#'                infers tau2_y i.e., scale parameter for homGP.
#'                \item \code{noise}: logical; if \code{noise = TRUE}, infers nugget
#'                \code{g} throught M-H for latent noise process.
#'              }
#' @param priors hyperparameters for priors and proposals (see details)
#' @param reps logical; if \code{reps = TRUE} uses Woodbury inference adjusting for
#'             replication of design points and \code{reps = FALSE} does not
#'             use Woodbury inference
#' @param cov covariance kernel, either Matern, ARD Matern 
#'        or squared exponential (\code{"exp2"})
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @param stratergy choose initialization stratergy; "default" uses hetGP for 
#'        \code{vecchia = FALSE} settings and sVecchia for \code{vecchia = TRUE}. 
#'        See details.
#' @param vecchia logical indicating whether to use Vecchia approximation
#' @param m size of Vecchia conditioning sets (only used if 
#'        \code{vecchia = TRUE})
#' @param ordering optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x}, defaults to random, is applied to \code{x}
#' @param verb logical indicating whether to print progress
#' @param omp_cores if \code{vecchia = TRUE}, user may specify the number of cores to
#'        use for OpenMP parallelization. Uses min(4, limit) where limit is max openMP 
#'        cores available on the machine.
#' @return a list of the S3 class \code{bhomgp} or \code{bhomgp_vec} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{Ylist}: list of all responses observed per location (x)
#'   \item \code{A}: number of replicates at each location
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{priors}: copy of proposal/prior settings
#'   \item \code{settings}: setting for predictions using \code{bhetgp} or \code{bhetgp_vec}
#'   object
#'   \item \code{g_y}: vector of MCMC samples for \code{g_y}
#'   \item \code{theta_y}: vector of MCMC samples for \code{theta_y} (length
#'         scale of mean process)
#'   \item \code{tau2}: vector of MAP estimates for \code{tau2} (scale 
#'         parameter of outer layer)
#'   \item \code{llik_y}: vector of MVN log likelihood of Y for reach Gibbs iteration
#'   \item \code{time}: computation time in seconds
#'   \item \code{x_approx}: conditioning set, NN and ordering for \code{vecchia = TRUE}
#'   \item \code{m}: copy of size of conditioning set for \code{vecchia = TRUE}      
#' }
#' 
#' @references 
#' Binois, Mickael, Robert B. Gramacy, and Mike Ludkovski. "Practical heteroscedastic Gaussian process 
#' modeling for large simulation experiments." Journal of Computational and Graphical 
#' Statistics 27.4 (2018): 808-821.
#' 
#' Katzfuss, Matthias, Joseph Guinness, and Earl Lawrence. "Scaled Vecchia approximation for 
#' fast computer-model emulation." SIAM/ASA Journal on Uncertainty Quantification 10.2 (2022): 537-554.
#' 
#' Sauer, Annie Elizabeth. "Deep Gaussian process surrogates for computer experiments." (2023).
#' 
#' @examples 
#' 
#' # 1D example with constant noise
#' 
#' # Truth
#' fx <- function(x){
#' result <- (6 * x - 2)^2* sin(12 * x - 4)
#' }
#' 
#' # Training data
#' r <- 10
#' xn <- seq(0, 1, length = 25)
#' x <- rep(xn, r)
#'
#' f <- fx(x) 
#' y <- f + rnorm(length(x)) # adding constant noise
#' 
#' # Testing data
#' xx <- seq(0, 1, length = 100)
#' yy <- fx(xx)
#' 
#---------------------------------------------------------------------------
#' # Example 1: Full model, no Vecchia
#---------------------------------------------------------------------------
#' 
#' # Fitting a bhomGP model using all the data
#' fit <- bhomGP(x, y, nmcmc = 100, verb = FALSE)
#' 
#' # Trimming the object to remove burn in and thin samples
#' fit <- trim(fit, 50, 10)
#' 
#' # Predition using the bhomGP object (indepedent predictions)
#' fit <- predict(fit, xx, lite = TRUE, cores = 2)
#' 
#' #' # Visualizing the mean predictive surface. 
#' # Can run plot(fit, trace = TRUE) to view trace plots
#' plot(fit) 
#' 
#---------------------------------------------------------------------------
#' # Example 2: Vecchia approximated model
#---------------------------------------------------------------------------
#' 
#' # Fitting a bhomGP model using vecchia approximation
#' fit <- bhomGP(x, y, nmcmc = 100, vecchia = TRUE, m = 5, omp_cores = 2, verb = FALSE)
#' 
#' # Trimming the object to remove burn in and thin samples
#' fit <- trim(fit, 50, 10)
#' 
#' # Predition using the bhomGP_vec object with Vecchia (indepedent predictions)
#' fit <- predict(fit, xx, vecchia = TRUE, cores = 2)
#' 
#' # Visualizing the mean predictive surface.
#' plot(fit) 
#'
#' @export
bhomGP <- function(x = NULL, y = NULL, reps_list = NULL, 
                   nmcmc = 1000,
                   sep = TRUE, inits = NULL, priors = NULL,
                   cov = c("exp2", "matern", "ARD matern"), v = 2.5,
                   stratergy = c("default", "flat"),
                   vecchia = FALSE, m = min(25, length(y) - 1),
                   ordering = NULL, reps = TRUE, verb = TRUE, omp_cores = 4){

  tic <- proc.time()[[3]]
  cov <- match.arg(cov)
  stratergy <- match.arg(stratergy)
  D = ifelse(is.matrix(x), ncol(x), 1)

  if (cov == "exp2") v <- 999 # indicator
  if (!vecchia & length(y) > 300)
    message("'vecchia = TRUE' is recommended for faster computation.")
  if (nmcmc <= 0) stop("nmcmc must be greater than 1")
  if (is.numeric(x)) x <- as.matrix(x)
  if (sep & ncol(x) == 1) sep <- FALSE # no need for separable theta

  # check reps for homGP case - faster computation
  if(reps){ # store mappings
    reps <- hetGP::find_reps(x, y)
    # YN <- reps$Z
    YNs2 <- unlist(lapply(reps$Zlist, function(x) sum((x - mean(x))^2)))
    Yn <- reps$Z0
    Xn <- reps$X0
    # Ylist <- reps$Zlist
    A <- reps$mult
  } else{
    reps <- NULL # indicator
    YN <- NULL # indicator
    YNs2 <- NULL
    Yn <- y
    Xn <- x
    # Ylist <- y
    A <- NULL
  }

  #check matern 0.5 -> weird results
  if (m >= length(y)) stop("m must be less than the length of y")
  if (cov == "matern" | cov == "ARD matern")
    if(!(v %in% c(1.5, 2.5)))
      stop("v must be one of 1.5, or 2.5")
  if(cov == "matern")
    v <- 1000 + v # indicator + v
  if (!is.null(ordering)) {
    if (!vecchia) message("ordering only used when vecchia = TRUE")
    test <- check_ordering(ordering, nrow(Xn)) # returns NULL if all checks pass
  }
  if(is.null(inits$mean_y)) inits$mean_y <- mean(y)
  
  # check priors and initialize
  
  priors <- check_settings(priors, gp_type = "homgp", x = Xn, y = Yn)
  # initial <- list(theta_y = theta_y_0, g_y = g_y_0, scale = 1)

  initial <- list(theta_y = inits$theta_y_0, g = inits$g_y_0,
                  mean_y = inits$mean_y, scale_y = inits$scale, prof_ll = inits$prof_ll, noise= inits$noise)

  initial <- check_initialization_hom(reps = ifel(is.null(reps), hetGP::find_reps(x, y), reps), initial, 
                                      n = nrow(Xn), D = D, sep = sep, v = v, vec = vecchia, 
                                      stratergy = stratergy, verb = verb)
  
  settings <- list(v = v, sep = sep, mean_y = initial$mean_y, verb = verb)
  
  out <- list(x = Xn, y = Yn, Ylist = reps$Zlist, A = reps$mult, 
              nmcmc = nmcmc, priors = priors, settings = settings)

  # priors = settings
  if(!vecchia){ # non vec homGP; YN = NULL will do no reps case
    
    calc <- precalc(Xn, out = Yn, outs2 = YNs2, A = A, v, vec = FALSE, priors, initial, 
                    latent = FALSE)
    initial$llik_y <- calc$llik
    initial$tau2_y <- calc$tau2
    
    test <- check_inputs(x, y) # returns NULL if all checks pass
    
    if (nmcmc == 1) {
      out <- c(out, initial)
      out$time <- proc.time()[[3]] - tic
      return(out)
    }
    
    if(sep) obj <- gibbs_sep_hom(YNs2 = YNs2, Yn = Yn,  Xn = Xn, A= A, mcmc = nmcmc,
                                 initial = initial, priors = priors, v = v, verb = verb)
    else obj <- gibbs_iso_hom(YNs2 = YNs2, Yn = Yn,  Xn = Xn, A= A, mcmc = nmcmc,
                              initial = initial, priors = priors, v = v, verb = verb)
  }
  else{ # vecchia
    out$m <- m
    out$ordering <- ordering
    x_approx <- create_approx(Xn, m, ordering, cores = omp_cores)

    if(!is.null(reps)){ # woodbury
      # YN_ord <- unlist(Ylist[x_approx$ord])
      YNs2_ord <- YNs2[x_approx$ord]
      Yn_ord <- Yn[x_approx$ord]
    } else{ # no reps i.e. full N
      YNs2_ord <- NULL
      Yn_ord <- Yn[x_approx$ord]
    }

    out$x_approx <- x_approx
    calc <- precalc(x_approx, out = Yn_ord, outs2 = YNs2_ord, A = A, v, 
                    vec = TRUE, priors, initial, latent = FALSE)
    initial$llik_y <- calc$llik
    initial$tau2_y <- calc$tau2
    
    test <- check_inputs(x, y) # returns NULL if all checks pass
    
    if (nmcmc == 1) {
      out <- c(out, initial)
      out$time <- proc.time()[[3]] - tic
      return(out)
    }
    
    # YN_ord = NULL will do no reps case
    if(sep) obj <- gibbs_sep_vec_hom(YNs2 = YNs2_ord, Yn = Yn_ord,  x_approx, A= A, mcmc = nmcmc,
                                     initial = initial, priors = priors, v = v, verb = verb)
    else obj <- gibbs_iso_vec_hom(YNs2 = YNs2_ord, Yn = Yn_ord,  x_approx, A= A, mcmc = nmcmc,
                                  initial = initial, priors = priors, v = v, verb = verb)
  }

  out <- c(out, obj)
  toc <- proc.time()[[3]]
  out$time <- unname(toc - tic)

  if (vecchia) class(out) <- "bhomgp_vec" else class(out) <- "bhomgp"
  
  return(out)
}


bhetGP_vdims <- function(x = NULL, y = NULL, reps_list = NULL, nmcmc = 1000, D = ifelse(is.matrix(x), ncol(x), 1),
                   sep = TRUE, inits = NULL, # ls_check = FALSE,
                   priors = NULL, reps = TRUE,
                   vdims = NULL, # does variance change in all dims
                   cov = c("exp2", "matern", "ARD matern"), v = 2.5, 
                   stratergy = c("default", "flat"), vecchia = FALSE,
                   m = min(25, length(y) - 1), ordering = NULL, verb = TRUE, omp_cores = 4){
  
  tic <- proc.time()[[3]]
  cov <- match.arg(cov)

  if (cov == "exp2") v <- 999 # indicator
  if (!vecchia & length(y) > 300)
    message("'vecchia = TRUE' is recommended for faster computation.")
  if (nmcmc <= 0) stop("nmcmc must be at least 1")
  if(!is.null(reps_list)){x <- reps_list$X0; y <- reps_list$Z0}
  
  if (is.numeric(x)) x <- as.matrix(x)
  test <- check_inputs(x, y) # returns NULL if all checks pass
  
  if(reps){ # woodbury version
    if(verb) print("obtaining replication")
    reps <- ifel(is.null(reps_list), hetGP::find_reps(x, y), reps_list)
    Yn <- reps$Z0
    Xn <- reps$X0
    xv <- reps$X0[, vdims, drop = FALSE]
    # Ylist <- reps$Zlist # for vecchia
    YNs2 <- unlist(lapply(reps$Zlist, function(x) sum((x - mean(x))^2)))
    A <- reps$mult
  } else{ # no reps
    reps <- NULL # no mapping
    Yn <- y
    Xn <- x
    xv <- x[, vdims, drop = FALSE]
    YNs2 <- NULL # indicator
    A <- NULL # indicator
  }
  
  # if variance does not change in every dimension
  if(!is.null(reps)){ # find mapping between X[, vdims] and Xn
    reps_vdims <- hetGP::find_reps(xv, 1:nrow(Xn))
    xv <- reps_vdims$X0
  }else{  # no reps version
    reps_vdims <- NULL # no mapping needed
  }  
  
  #check matern 0.5 -> weird results
  if (m >= length(y)) stop("m must be less than the length of y")
  if (cov == "matern" | cov == "ARD matern")
    if(!(v %in% c(1.5, 2.5)))
      stop("v must be one of 1.5, or 2.5")
  if(cov == "matern")
    v <- 1000 + v # indicator + v
  if (!is.null(ordering)) {
    if (!vecchia) message("ordering only used when vecchia = TRUE")
    test <- check_ordering(ordering, nrow(Xn)) # returns NULL if all checks pass
  }
  
  # check prior settings
  priors <- check_settings(priors, gp_type = "hetgp", x = Xn, y = Yn)
  
  if(is.null(inits$mean_y)) inits$mean_y <- mean(y)
  initial <- list(theta_y = inits$theta_y_0, theta_lam = inits$theta_lam_0, g = inits$g_0,
                  llam = inits$llam_0, mean_y =inits$mean_y , mean_lam = inits$mean_lam,
                  scale_lam = inits$scale_lam, scale_y = inits$scale_y, theta_check = inits$ls_check,
                  prof_ll_lam = inits$prof_ll_lam, noise = inits$noise)
  
  initial <- check_inits_vdims(reps = ifel(is.null(reps), hetGP::find_reps(x, y), reps), 
                               initial, n = nrow(Xn), D = D, sep = sep,
                               vdims = vdims, nlam = nrow(xv), 
                               reps_vdims = ifel(is.null(reps_vdims), hetGP::find_reps(xv, y), reps_vdims), 
                               v = v, stratergy = stratergy, vec = vecchia, verb = verb)
  
  settings <- list(v = v, vdims = vdims, sep = sep, mean_y = initial$mean_y,
                   mean_lam = initial$mean_lam)
  
  mappings <- list(reps_vdims = reps_vdims) # store mappings for predictions
  
  out <- list(x = Xn, y = Yn, Ylist = reps$Zlist, A = reps$mult, nmcmc = nmcmc, 
              priors = priors, settings = settings)
  
  if(!vecchia){ # non vecchia
    initial <- check_scale_vdims(x = xv, v = v, init = initial, vec = FALSE, sep = sep)
    if (nmcmc == 1) {
      out <- c(out, initial)
      return(out)
    }
    if(verb) print("starting mcmc")
    if(is.null(reps)){ # no reps
      if(sep) obj <- gibbs_sep_vdims_N(Yn = Yn,  Xn, nmcmc,  initial, 
                                       priors = priors, v, vdims = vdims, verb = verb)
      else obj <- gibbs_iso_vdims_N(Yn = Yn,  Xn, nmcmc, initial, 
                                    priors = priors, v, vdims = vdims, verb = verb)
    }
    else{ # reps (woodbury mtd)
      if(sep) obj <- gibbs_sep_vdims(YNs2= YNs2, Yn = Yn,  Xn,  A= A, nmcmc,
                                     initial, priors = priors, v, vdims = vdims, reps_vdims, verb = verb)
      else obj <- gibbs_iso_vdims(YNs2 = YNs2, Yn = Yn,  Xn,  A= A, nmcmc,
                                  initial, priors = priors, v, vdims = vdims, reps_vdims, verb = verb)
    }
  } else{ # vecchia
    out$m <- m # neighbours
    out$ordering <- ordering # order
    x_approx <- create_approx(Xn, m, ordering, cores = omp_cores) # create approximation
    
    if(is.null(reps))
      YNs2_ord <- NULL # If no reps, use Xn and Yn (no YN)
    else
      YNs2_ord <- YNs2[x_approx$ord]
    # YN_ord <- unlist(Ylist[x_approx$ord]) # order Ylists by NN and then unlist for YN
    
    Yn_ord <- Yn[x_approx$ord] # order Yn
    out$x_approx <- x_approx # store approximation
    
    initial <- check_scale_vdims(x = x_approx, v = v, init = initial, vec = TRUE, sep = sep)
    if (nmcmc == 1) {
      out <- c(out, initial)
      return(out)
    }
    if(verb) print("starting mcmc")
    if(is.null(reps)){ # if no reps
      if(sep) obj <- gibbs_sep_vdims_N_vec(Yn = Yn_ord,  x_approx, nmcmc,  initial, 
                                           priors = priors, v, vdims = vdims, verb = verb)
      else obj <- gibbs_iso_vdims_N_vec(Yn = Yn_ord,  x_approx, nmcmc, initial, 
                                        priors = priors, v, vdims = vdims, verb = verb)
    }
    else{ # if reps (woodbury)
      if(sep) obj <- gibbs_sep_vdims_vec(YNs2_ord, Yn_ord, x_approx, xv, reps$mult, vdims = vdims, nmcmc,
                                         initial, priors = priors, v = v, reps_vdims, verb = verb)
      else obj <- gibbs_iso_vdims_vec(YNs2_ord, Yn_ord, x_approx, xv, reps$mult, vdims = vdims, nmcmc,
                                      initial, priors = priors, v = v, reps_vdims, verb = verb)
    }
  }
  
  out <- c(out, obj)
  toc <- proc.time()[[3]]
  out$time <- unname(toc - tic)
  
  if (vecchia) class(out) <- "bhetgp_vec" else class(out) <- "bhetgp"
  
  return(out)
}


