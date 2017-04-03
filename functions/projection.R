# 
# functions for performing the projection predictive variable subset selection
# for linear Gaussian model.
#
#####---------------------------------------------------------------------------------------------------------#####
##### R script to calculate KL for a submodel                                                                 #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       w                 :=  coefficient matrix: (P+1) * niter(posterior draws)                          #####
#####       sigma2            :=  sigma square vector: length = niter                                         #####
#####       x                 :=  predictor matrix containing all predictors: N * (P+1)                       #####
#####       indproj           :=  indicators of the predictors selected in the model, a subset of 1:(P+1)     #####
##### Output:                                                                                                 #####
#####       w                 :=  LSE of coefficients in the submodel                                         #####
#####       sigma2            :=  mse of submodel                                                             #####
#####       kl                :=  a vector of KL (Kullback-Leibler divergence) of all submodels of length P+1 #####
#####---------------------------------------------------------------------------------------------------------#####
lm_proj <- function(w,sigma2,x,indproj) {
  
  # assume the intercept term is stacked into w, and x contains 
  # a corresponding vector of ones. returns the projected samples
  # and estimated kl-divergence.
  
  # pick the columns of x that form the projection subspace
  n <- dim(x)[1]
  xp <- x[,indproj]
  
  # solve the projection equations
  fit <- x %*% w # fit of the full model
  wp <- solve(t(xp) %*% xp, t(xp) %*% fit,tol=1e-20)
  sigma2p <- sigma2 + colMeans((fit - xp %*% wp)^2)
  
  # this is the estimated kl-divergence between the full and projected model
  kl <- mean(0.5*log(sigma2p/sigma2))
  
  # reshape wp so that it has same dimensionality as x, and zeros for
  # those variables that are not included in the projected model
  d <- dim(w)[1]
  S <- dim(w)[2]
  wptemp <- matrix(0, d, S)
  wptemp[indproj,] <- wp
  wp <- wptemp
  
  return(list(w=wp, sigma2=sigma2p, kl=kl))
}


#####---------------------------------------------------------------------------------------------------------#####
##### R script to generate all models selected from a stan object                                             #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       w                 :=  coefficient matrix: (P+1) * niter(posterior draws)                          #####
#####       sigma2            :=  sigma square vector: length = niter                                         #####
#####       x                 :=  predictor matrix containing all predictors: N * (P+1)                       #####
##### Output:                                                                                                 #####
#####       chosen            :=  a vector of indices of predictors ordered by the priority of entering models#####
#####                             of length P+1                                                               #####
#####       kl                :=  a vector of KL (Kullback-Leibler divergence) of all submodels of length P+1 #####
#####---------------------------------------------------------------------------------------------------------#####
lm_fprojsel <- function(w, sigma2, x) {
  
  # forward variable selection using the projection
  
  d = dim(x)[2]
  chosen <- c() # chosen variables, start from the model with the intercept only
  notchosen <- setdiff(1:d, chosen)
  
  # start from the model having only the intercept term
  kl <- rep(0,d)
  # kl[1] <- lm_proj(w,sigma2,x,1)$kl
  
  # start adding variables one at a time
  for (k in 1:d) {
    
    nleft <- length(notchosen)
    val <- rep(0, nleft)
    
    for (i in 1:nleft) {
      ind <- sort( c(chosen, notchosen[i]) )
      proj <- lm_proj(w,sigma2,x,ind)
      val[i] <- proj$kl
    }
    
    # find the variable that minimizes the kl
    imin <- which.min(val)
    chosen <- c(chosen, notchosen[imin])
    notchosen <- setdiff(1:d, chosen)
    
    kl[k] <- val[imin]
  }
  return(list(chosen=chosen, kl=kl))
}


#####---------------------------------------------------------------------------------------------------------#####
##### R script to generate all models selected at group level from a stan object                              #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       w                 :=  coefficient matrix: (P+1) * niter(posterior draws)                          #####
#####       sigma2            :=  sigma square vector: length = niter                                         #####
#####       x                 :=  predictor matrix containing all predictors: N * (P+1)                       #####
#####       group_index       :=  a vector of length P indicating the group each predictor belongs to         #####
##### Output:                                                                                                 #####
#####       chosen            :=  a vector of indices of groups ordered by the priority of entering models    #####
#####                             of length G (the number of groups)                                          #####
#####       kl                :=  a vector of KL (Kullback-Leibler divergence) of all submodels of length G   #####
#####---------------------------------------------------------------------------------------------------------#####


lm_fprojsel_grp <- function(w, sigma2, x, group_index) {
  
  group_names <- unique(group_index)
  G = length(group_names)
  
  # forward group selection using the projection
  d = dim(x)[2]
  chosen <- c() # chosen groups, start from the model with the intercept only
  notchosen <- setdiff(1:G, chosen)
  
  # start from the model having only the intercept term
  kl <- rep(0,G)
  # kl[1] <- lm_proj(w,sigma2,x,1)$kl
  
  # start adding variables one at a time
  for (g in 1:G) {
    
    nleft <- length(notchosen)
    val <- rep(0, nleft)
    
    for (i in 1:nleft) {
      ind <- sort(which(group_index %in% c(chosen, notchosen[i])))
      proj <- lm_proj(w,sigma2,x,ind)
      val[i] <- proj$kl
    }
    
    # find the variable that minimizes the kl
    imin <- which.min(val)
    chosen <- c(chosen, notchosen[imin])
    notchosen <- setdiff(1:d, chosen)
    
    kl[g] <- val[imin]
  }
  return(list(chosen=chosen, kl=kl))
}

