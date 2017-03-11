## strong normal prior ridge regression, lambda as a parameter

data {
  int<lower=0> N;
  int<lower=1> K;
  vector[N] y;
  matrix[N,K] x;
}

parameters {
  vector[K] beta;
  real<lower=0> lambda;
}

transformed parameters {
  real<lower=0> squared_error;
  squared_error = dot_self(y - x * beta);
}

model {
  target += -squared_error;
  target += - lambda * dot_self(beta); 
  beta ~ normal(0, 1);
}
