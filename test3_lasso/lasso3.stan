## Laplace prior, lambda as a parameter

data {
  int<lower=0> N;
  int<lower=1> K;
  vector[N] y;
  matrix[N,K] x;
}

parameters {
  vector[K] beta;
  real mu;
  real<lower=0> sigma;
  real<lower=0> lambda;
  real<lower=0> gamma;
}

model {
  y ~ normal(x * beta, sigma);
  sigma ~ normal(0,1);
  beta ~ double_exponential(mu, gamma);
  for (k in 1:K)
    target += - lambda * N * fabs(beta[k]); 
}
