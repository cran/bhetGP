# Imported Functions ----------------------------------------------------------
#' @importFrom Matrix t solve
#' @importFrom grDevices heat.colors
#' @importFrom graphics image lines matlines par plot points contour abline matplot
#' @importFrom stats cov dgamma dnorm pnorm qnorm rnorm runif var loess predict optimize sd
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm rmvnorm
#' @importFrom FNN get.knnx
#' @importFrom GpGp find_ordered_nn matern35_scaledim matern45_scaledim matern15_scaledim exponential_scaledim matern25_scaledim
#' @importFrom GPvecchia order_maxmin_exact order_middleout order_coordinate order_dist_to_point
#' @importFrom laGP darg
#' @importFrom hetGP find_reps mleHetGP mleHomGP 

# Package Documentation -------------------------------------------------------
#' @useDynLib bhetGP, .registration = TRUE
#' @title Package bhetGP
#' @author Parul V. Patil \email{parulvijay@vt.edu}
#' @name bhetGP-package
#'
#' @description Performs Bayesian posterior inference for heteroskedastic Gaussian processes.
#' Models are trained through MCMC including elliptical slice sampling (ESS) of 
#' latent noise processes and Metropolis-Hastings sampling of 
#' kernel hyperparameters. Replicates are handled efficientyly through a
#' Woodbury formulation of the joint likelihood for the mean and noise process 
#' (Binois, M., Gramacy, R., Ludkovski, M. (2018) <doi:10.1080/10618600.2018.1458625>)
#' For large data, Vecchia-approximation for faster 
#' computation is leveraged (Sauer, A., Cooper, A., and Gramacy, R.,
#' (2023), <doi:10.1080/10618600.2022.2129662>). Incorporates 'OpenMP' and 
#' SNOW parallelization and utilizes 'C'/'C++' under the hood.
#' 
#' @section Important Functions:
#' \itemize{
#'   \item \code{\link[bhetGP]{bhetGP}}: conducts MCMC sampling of 
#'   hyperparameters and latent noise layer for a heteroskedatic GP.
#'   \item \code{\link[bhetGP]{bhomGP}}: conducts MCMC sampling of 
#'   hyperparameters for a homoskedastic GP.
#'   \item \code{\link[bhetGP]{trim}}: cuts off burn-in and optionally thins 
#'   samples
#'   \item \code{\link[bhetGP]{predict}}: calculates posterior mean and 
#'   variance over a set of input locations (optionally calculates EI or entropy)
#'   \item \code{\link[bhetGP]{plot}}: produces trace plots, hidden layer 
#'   plots, and posterior predictive plots
#' }
#' 
#' @references 
#' M. Binois, Robert B. Gramacy, M. Ludkovski (2018), Practical heteroskedastic 
#' Gaussian process modeling for large simulation experiments,
#' Journal of Computational and Graphical Statistics, 27(4), 808--821.
#' 
#' Katzfuss, Matthias, Joseph Guinness, and Earl Lawrence. 
#' Scaled Vecchia approximation for fast computer-model emulation. 
#' SIAM/ASA Journal on Uncertainty Quantification 10.2 (2022): 537-554.
#' 
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
#' processes for computer experiments. 
#' Journal of Computational and Graphical Statistics, 32(3), 824-837.  
#' 
#' @examples
#' # See ?bhetGP, or ?bhomGP for examples
#' 
"_PACKAGE"