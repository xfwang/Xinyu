/*  Variable naming: 
 obs       = observed 
 cen       = (right) censored 
 N         = number of samples 
 tau       = scale parameter 
*/ 
data { 
  int<lower=0> Nobs; 
  int<lower=0> Ncen; 
  vector[Nobs] yobs; 
  vector[Ncen] ycen; 
  int<lower=1> K;  // number of covariates
  matrix[Nobs,K] xobs;
  matrix[Ncen,K] xcen;
} 

transformed data { 
  real<lower=0> tau_al; 
  tau_al = 10.0; 
} 
 
parameters { 
  vector[K] beta;
  real alpha_raw; 
    real<lower=0> lambda;
  real mu;
  real<lower=0> sigma;
} 
 
transformed parameters { 
  real alpha; 
  vector[Nobs] mu_obs; 
  vector[Ncen] mu_cen; 
  alpha = exp(tau_al * alpha_raw); 
  mu_obs = xobs * beta;
  mu_cen = xcen * beta;
} 
 
model { 
  yobs ~ weibull(alpha, exp(-(mu_obs)/alpha)); 
  target += weibull_lccdf(ycen | alpha, exp(-(mu_cen)/alpha)); 
    for (k in 1:K)
      target += - lambda * fabs(beta[k]); 
  beta ~ double_exponential(mu, sigma);
  alpha_raw ~ normal(0.0, 1.0); 
} 
