data {
  int<lower=0> N;
  int<lower=1> K;
  vector[N] y;
  matrix[N,K] x;
  real<lower=0> lambda;
}

parameters {
  vector[K] beta;
}

transformed parameters {
  real<lower=0> squared_error;
  real<lower=0> S1;
  real<lower=0> S2;
  real<lower=0> S3;

  squared_error = dot_self(y - x * beta);
  S1 = dot_self(beta[1:4]);
  S2 = dot_self(beta[5:7]);
  S3 = dot_self(beta[8:10]);
}

model {
  target += -squared_error;
  target += - lambda * (S1+S2+S3); 
}

generated quantities {
  real<lower=0> sigma_squared;
  sigma_squared = squared_error / N;
}
