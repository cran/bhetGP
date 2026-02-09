This repo consists of an in progress interface for the implementation of a Bayesian HetGP

Run R CMD INSTALL bhGP from outside this folder. If already installed and package is updated: remove package, restart R, and then run this install command. 


R files include:

- **functions.R** : all functions for log likelihood and Sigma contruction

- **mcmc.R** : Does mcmc and ess sampling for theta, nugget vector, etc.

- **gibbs.R** : Gibbs sampler for MCMC for seperable, iso-tropic, and combination

- **fit.R** : Fits bhGP to the data (X,Y)

- **predict.R**: Predicts at new locations Xp

- **mcmc_vecchia.R** : Does mcmc and ess sampling for theta, nugget vector, w/ vecchia

- **gibbs_vecchia.R** : Gibbs sampler for MCMC for seperable, iso-tropic, and combination w/ vecchia

- **gibbs_vdims_N.R**: Gibbs for full N version w/ variance changing in some dimenstions

- **predict_vecchia.R**: Predicts at new locations Xp w/ vecchia

- **trim.R**: trims mcmc object

- **plots.R**: Trace plots and prediction plots

- **vecchia.R**: Helper vecchia functions

- **general.R**: Helper functions

- **checks.R** : checks initilizations

- **inits.R**: initialization stratergy 

src files include:

- **vecchia.cpp**: Vecchia related fucntions

- **cov.cpp**: Builds cov matrices

- **invdet.c**: inv and det fns

- **general.c**: used for invdet 

- **RcppExports.cpp**: created for package

## Updates

- vdims is working now

- removed sq from matern for SV

- predict lam_ub is now by default FALSE

- stratergy = "flat" now working

- priors on theta_lam adjusted ---> ocean example working out better

- ls_check = TRUE by default??? maybe not. 


- checks.R --> check_scale() ---> remove print statement
- check all prints before reupload.


------------------------- Next to Dos --------------------------------------

- for nmcmc = 1, make predictions also possible (to be done)

- will need to make joint U preds faster.

## Version 1.0.2

- corrects span error
- fixed stratergy for initialization


