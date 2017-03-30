#####---------------------------------------------------------------------------------------------------------#####
##### R script to simulate a predictor matrix for group LASSO (X)                                             #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       n           :=  the number of samples                                                             #####
#####       p           :=  the number of covariates                                                          #####
#####       cor         :=  p*p correlation matrix                                                            #####
#####       sd          :=  the standard deviation of covariates, a number rather than a vector               #####
#####       random.seed :=  random seed                                                                       #####
##### Output:                                                                                                 #####
#####       x           :=  simulated covariance matrix, n*p                                                  #####
#####---------------------------------------------------------------------------------------------------------#####
X_groupLasso <- function(n, p, cor = diag(p), sd = 1, random.seed = NULL) {
  if(!is.null(random.seed)) set.seed(random.seed)
  if(!(ncol(cor) == nrow(cor) & ncol(cor) == p)) stop("the dimension of correlation matrix does not match")
  CovX = cor*(sd)^2
  x = rmvnorm(n, sigma=CovX, method="chol")
  return(x)
}






#####---------------------------------------------------------------------------------------------------------#####
##### R script to create a correlation matrix for group LASSO (corr)                                          #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       K           :=  the number of groups                                                              #####
#####       p0          :=  the number of predictor variables in each group                                   #####
#####       rho         :=  the correlation coefficient within each group                                     #####
##### Output:                                                                                                 #####
#####       corr        :=  correlation matrix, (K*p0) * (K*p0)                                               #####
#####---------------------------------------------------------------------------------------------------------#####
corr_create <- function(K, p0, rho=0) {
  corr <- matrix(0, nrow = K*p0, ncol = K*p0)
  for(k in seq(K)) {
    corr[((k-1)*p0+1):(k*p0), ((k-1)*p0+1):(k*p0)] = rho
  }
  diag(corr) = 1
  return(corr)
}










#####       eff.size    :=  the effect size of active covariates, a number or a vector                        #####
#####       G           :=  the number of groups of covariates, p/G covariates in each group                  #####
