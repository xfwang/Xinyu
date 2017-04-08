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
  # if(!is.matrix(beta.list)) beta.list = matrix(beta.list, ncol = 1)
  # N = length(y)
  # K = ncol(x)
  # if(nrow(x) != N) stop("dimension mismatch in x and y") else
  #   if(nrow(beta.list) != K) stop("dimension mismatch in x and beta")
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
  # if(!is.matrix(beta.list)) beta.list = matrix(beta.list, ncol = 1)
  N = length(y)
  # K = ncol(x)
  # if(nrow(x) != N) stop("dimension mismatch in x and y") else
  #   if(nrow(beta.list) != K) stop("dimension mismatch in x and beta") else
  #     if(length(sigma.list) != ncol(beta.list)) stop("length of iterations does not match in beta and sigma")
  sigma2p <- rep(sigma.list, each=N)
  pd <- dnorm(y, x %*% beta.list, sqrt(sigma2p))
  return(mean(log(rowMeans(pd))))
}

#####---------------------------------------------------------------------------------------------------------#####
##### R script for genereating various predictive performance measures for binary classifiers                 #####
#####                                                                                                         #####
##### Input:                                                                                                  #####
#####       y    := vector of true binary labels; 1/0                                                         #####
#####       yhat := vector of predicted binary labels; 1/0                                                    #####
##### Output:                                                                                                 #####
#####       output := a 1x5 vector for ppv, npv, sen, spe, acc, auc, (mcc is not included in this simulation) #####
#####---------------------------------------------------------------------------------------------------------#####

binary.classification.evaluation  = function(y, yhat) {
  
  ppv = 0
  npv = 0
  sen = 0
  spe = 0
  acc = 0
  ### use the pROC package, we could calculate the auc ###
  
  ### 2013-06-07 tp=true positive (y=1 & yhat=1)  ###
  tp = sum((y + yhat) == 2)/length(y)
  ### 2013-06-07 tn=true negative (y=0 & yhat=0)  ###
  tn = sum((y + yhat) == 0)/length(y)
  ### 2013-06-07 fp=false positive (y=0 & yhat=1) ###
  fp = sum((y - yhat) == -1)/length(y)
  ### 2013-06-07 fn=fase negative (y=1 & yhat=0)  ###
  fn = sum((y - yhat) == 1)/length(y)
  
  ### 2013-06-07 ppv=positive predictive value ###
  ### 2013-06-07 For biomarkers to predict good response: expect HIGH PPV ###
  ppv = tp/(tp+fp)
  ### 2013-06-07 npv=negative predictive value ###
  npv = tn/(tn+fn)
  ### 2013-06-07 sen=sensitivity=number of true positive/total number of positive ###
  ### = number of true positive/(number of true positive+number of false negative)###
  sen = tp/(tp+fn)
  ### 2013-06-07 spe=specificity=number of true negative/total number of negative ###
  ### = number of true negative/(number of true negative+number of false positive)###
  spe = tn/(tn+fp)
  
  ### 2013-06-07 acc=accuracy=(number of tp+number of tn)/(number of tp +number of tn +number of fp +number of fn) ###
  acc = (tn+tp)/(tp+tn+fp+fn)
  
  ### 2013-06-07 mcc=Matthews correlation coefficient=correlation coefficient between the observed and predicted binary classification ###
  #  mcc = (tp*tn - fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fn)*(tn+fp))
  
  output= matrix(c(ppv, npv, sen, spe, acc), ncol = 5, nrow = 1)
  
  return(output)
}











