## group lasso, Laplace prior, lambda as a PARAMETER

data {
    int<lower=0> N;
    int<lower=1> K;
    vector[N] y;
    matrix[N,K] x;
    }
    
    parameters {
    real<lower=0> lambda;
    vector[K] beta;
    }
    
    transformed parameters {
    real<lower=0> squared_error;
    real S1;
    real S2;
    real S3;
    vector[10] beta1;
    vector[10] beta2;
    vector[80] beta3;

    squared_error = dot_self(y - x * beta);
    for(i in 1:10) beta1[i] = beta[i];
    for(i in 11:20) beta2[i-10] = beta[i];
    for(i in 21:100) beta3[i-20] = beta[i];
    S1 = dot_self(beta1);
    S2 = dot_self(beta2);
    S3 = dot_self(beta3);
    }

model {
  beta ~ double_exponential(0,1);
  target += -squared_error;
  target += - lambda * (S1+S2+S3); 
}

generated quantities {
  real<lower=0> sigma_squared;
  sigma_squared = squared_error / N;
}
