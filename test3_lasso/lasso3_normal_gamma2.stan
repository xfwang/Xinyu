## Gasussian_Gamma prior, lambda as a parameter

data {
  int<lower=0> N;
  int<lower=1> K;
  vector[N] y;
  matrix[N,K] x;
}

parameters {
  vector[K] beta;
  real<lower=0> tau[K];
  real mu;
  real<lower=0> sigma;
  real<lower=0> lambda;
  // real<lower=1> a;
  // real<lower=0> b;
}

transformed parameters {
  real<lower=0> lambda_sq;
  real<lower=0> tau_sq[K];
  lambda_sq = lambda ^2;
  for (j in 1:K)
    tau_sq[j] = tau[j] ^ 2;
}

model {
  y ~ normal(x *beta, sigma);
  for (j in 1:K)
    beta[j] ~ normal(mu, sigma*tau[j]);
  tau_sq ~ exponential(lambda_sq * N^2 / 8);
  // lambda_sq ~ gamma(a,b);
  // a ~ gamma(1,1);
  for (k in 1:K)
    target += - lambda * N * fabs(beta[k]); 
  target += sum(log(tau));
}
