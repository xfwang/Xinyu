## non-informative prior  

data {
  int<lower=0> N;
  int<lower=1> K;
  vector[N] y;
  matrix[N,K] x;
}

parameters {
  vector[K] beta;
}

transformed parameters {
  real<lower=0> squared_error;
  squared_error = dot_self(y - x * beta);
}

model {
  target += -squared_error;
}
