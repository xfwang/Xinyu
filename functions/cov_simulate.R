#####---------------------------------------------------------------------------------------------------------#####
##### R script to simulate a covariate matrix for group LASSO                                                 #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       n           :=  the number of samples                                                             #####
#####       p           :=  the number of covariates                                                          #####
#####       cor         :=  p*p correlation matrix                                                            #####
#####       sd          :=  the standard deviation of covariates                                              #####
#####       random.seed :=  random seed                                                                       #####
##### Output:                                                                                                 #####
#####       x           :=  simulated covariance matrix, n*p                                                  #####
#####---------------------------------------------------------------------------------------------------------#####
cov_groupLasso <- function(n, p, cor = diag(p), sd = 1, random.seed = NULL) {
  if(!is.null(random.seed)) set.seed(random.seed)
  CovX = cor*(sd)^2
  x = rmvnorm(n, sigma=CovX, method="chol")
  return(x)
}

















#####       eff.size    :=  the effect size of active covariates, a number or a vector                        #####
#####       G           :=  the number of groups of covariates, p/G covariates in each group                  #####
