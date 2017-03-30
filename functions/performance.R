#####---------------------------------------------------------------------------------------------------------#####
##### R script to calculate mean squared error of an estimated beta given x and y                             #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       beta              :=  the esimated coefficients                                                   #####
#####       x                 :=  the covariates (usually in test dataset )                                   #####
#####       y                 :=  the response   (usually in test dataset )                                   #####
##### Output:                                                                                                 #####
#####       mse               :=  a real number, mean((y - yhat)^2)                                           #####
#####---------------------------------------------------------------------------------------------------------#####
# mse.fnc <- function(beta, x, y) {
#   return(mean((y - x %*% beta) ^ 2))
# }

#####---------------------------------------------------------------------------------------------------------#####
##### R script to calculate mean squared error of an estimated beta given x and y                             #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       beta.list         :=  all iterations of beta values (number of parameters * iterations)           #####
#####       x                 :=  the covariates (usually in test dataset )                                   #####
#####       y                 :=  the response   (usually in test dataset )                                   #####
##### Output:                                                                                                 #####
#####       mse               :=  a real number, mean((y - yhat)^2)                                           #####
#####---------------------------------------------------------------------------------------------------------#####
mse.fnc <- function(beta.list, x, y) {
  if(!is.matrix(beta.list)) beta.list = matrix(beta.list, ncol = 1)
  N = length(y)
  K = ncol(x)
  if(nrow(x) != N) stop("dimension mismatch in x and y") else
    if(nrow(beta.list) != K) stop("dimension mismatch in x and beta")
  return(mean((y-rowMeans(x %*% beta.list))^2))
}

#####---------------------------------------------------------------------------------------------------------#####
##### R script to calculate mean log predictive density (mlpd)                                                #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       beta.list         :=  all iterations of beta values (number of parameters * iterations)           #####
#####       sigma.list        :=  all iterations of sigma values (length = number of iterations)              #####
#####       x                 :=  the covariates (usually in test dataset )                                   #####
#####       y                 :=  the response   (usually in test dataset )                                   #####
##### Output:                                                                                                 #####
#####       mlpd              :=  a real number, mean(log(p(y|x, beta, sigma)))                               #####
#####---------------------------------------------------------------------------------------------------------#####
mlpd.fnc <- function(beta.list, sigma.list, x, y) {
  if(!is.matrix(beta.list)) beta.list = matrix(beta.list, ncol = 1)
  N = length(y)
  K = ncol(x)
  if(nrow(x) != N) stop("dimension mismatch in x and y") else
    if(nrow(beta.list) != K) stop("dimension mismatch in x and beta") else
      if(length(sigma.list) != ncol(beta.list)) stop("length of iterations does not match in beta and sigma")
  sigma2p <- rep(sigma.list, each=N)
  pd <- dnorm(y, x %*% beta.list, sqrt(sigma2p))
  return(mean(log(rowMeans(pd))))
}












