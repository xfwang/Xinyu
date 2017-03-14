## Laplace prior, lambda as a constant

data {
  int<lower=0> N;
  int<lower=1> K;
  int y[N];
  matrix[N,K] x;
  real<lower=0> lambda;
}

parameters {
  vector[K] beta;
  real m;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] prob;
  vector[N] mu;
  mu = x * beta;
  for (n in 1:N)
    prob[n] = 1 / (1 + exp(- mu[n]));
}

model {
  beta ~ double_exponential(m, sigma);
  y ~ bernoulli(prob);
  for (k in 1:K)
    target += - lambda * fabs(beta[k]); 
}
