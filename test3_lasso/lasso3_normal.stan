## Gasussian prior, lambda as a parameter

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
  real<lower=0> gamma;
  real<lower=0> lambda;
}

model {
  y ~ normal(x *beta, sigma);
  beta ~ normal(mu, gamma);
  for (k in 1:K)
    target += - lambda * N * fabs(beta[k]); 
}
