library(CoxBoost); library(dplyr)
source("./functions/result_process.R")
sapply(paste("./functions", dir(paste0("./functions")),sep="\\"),source)

# Generate some survival data with 10 informative covariates
n <- 200; K <- 50
beta <- c(rep(1.5,10),rep(0,K-10))
x <- matrix(rnorm(n*K),n,K)
real.time <- -(log(runif(n)))/(10*exp(drop(x %*% beta)))
cens.time <- rexp(n,rate=1/10)
status <- ifelse(real.time <= cens.time,1,0)
obs.time <- ifelse(real.time <= cens.time,real.time,cens.time)

# Fit a Cox proportional hazards model by CoxBoost
cbfit <- CoxBoost(time=obs.time,status=status,x=x,stepno=100,penalty=100)
summary(cbfit)

# Fit a stan model
stan_file <- "test4_survival/coxph.stan"
stan_data <- list(
  Nobs = sum(status == 1),
  Ncen = sum(status == 0),
  yobs = obs.time[status == 1],
  ycen = obs.time[status == 0],
  K = K,
  xobs = x[status == 1,],
  xcen = x[status == 0,]
)

## MODEL1 random initial values
stan.fit <- stan(stan_file,
                 data = stan_data,
                 chains = 1,
                 pars = c("beta","alpha","lambda"), 
                 iter = 1000,
                 seed = 1328025050
)
crossprod(param.stan(stan.fit)-beta)
which(param.stan(stan.fit, param = "beta") != 0)


library(psbcGroup)
# generate some survival data
set.seed(204542)
p = 20
n = 100
beta.true <- c(rep(4, 10), rep(0, (p-10)))
CovX<-matrix(0,p,p)
for(i in 1:10){
  for(j in 1:10){
    CovX[i,j] <- 0.5^abs(i-j)
  }
}
diag(CovX) <- 1
survObj <- list()
survObj$x <- apply(rmvnorm(n, sigma=CovX, method="chol"), 2, scale)
pred <- as.vector(exp(rowSums(scale(survObj$x, center = FALSE, scale = 1/beta.true))))
t <- rexp(n, rate = pred)
cen <- runif(n, 0, 8)
survObj$t <- pmin(t, cen)
survObj$di <- as.numeric(t <= cen)
priorPara <- list()
priorPara$eta0 <- 1
priorPara$kappa0 <- 1
priorPara$c0 <- 2
priorPara$r <- 0.5
priorPara$delta <- 0.0001
priorPara$s <- sort(survObj$t[survObj$di == 1])
priorPara$s <- c(priorPara$s, 2*max(survObj$t)
                 -max(survObj$t[-which(survObj$t==max(survObj$t))]))
priorPara$J <- length(priorPara$s)
priorPara$groupInd <- c(rep(1,10),2:11)
mcmcPara <- list()
mcmcPara$numBeta <- p
mcmcPara$beta.prop.var <- 1
initial <- list()
initial$beta.ini <- rep(0.5, p)
initial$lambdaSq <- 1
initial$sigmaSq <- runif(1, 0.1, 10)
initial$tauSq <- rexp(length(unique(priorPara$groupInd)),
                      rate = initial$lambdaSq/2)
initial$h <- rgamma(priorPara$J, 1, 1)
rw = FALSE
num.reps = 20000
chain = 1
thin = 5
save = 5
fitGL <- psbcGL(survObj, priorPara, initial, rw=FALSE, mcmcPara,
                num.reps, thin, chain, save)
## End(Not run)