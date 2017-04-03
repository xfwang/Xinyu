## Laplace prior, lambda as a constant

data {
  int<lower=0> N;
  int<lower=1> K;
  vector[N] y;
  matrix[N,K] x;
  real<lower=0> lambda;
}

parameters {
  vector[K] beta;
  real mu;
  real<lower=0> sigma;
  real<lower=0> gamma;
}

model {
  y ~ normal(x * beta, sigma);
  beta ~ double_exponential(mu, gamma);
  for (k in 1:K)
    target += - lambda * fabs(beta[k]); 
}
