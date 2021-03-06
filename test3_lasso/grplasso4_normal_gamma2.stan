## group lasso, Laplace prior, lambda as a PARAMETER, with hyperparameter

data {
    int<lower=0> N; // number of samples
    int<lower=1> K; // number of covariates in each group 
    int<lower=1> G; // number of groups
    vector[N] y;
    matrix[N,K*G] x;
    }
    
parameters {
    real<lower=0> lambda;
    vector[K*G] beta;
    vector[G] mu;
    real<lower=0> sigma[G];
    }
    
transformed parameters {
    real<lower=0> squared_error;
    vector[G] SS;
    real<lower=0> sigma_sq[G];

    squared_error = dot_self(y - x * beta);
    for(i in 1:G) {
      SS[i] = dot_self(beta[((i-1)*K+1) : (i*K)]);
      sigma_sq[i] = sigma[i];
    }
    }

model {
    for(i in 1:G) 
      beta[((i-1)*K+1) : (i*K)] ~ normal(mu[i], sigma[i]);
    mu ~ normal(0, 1);
    sigma_sq ~ exponential(lambda^2 * N^2);
    target += -squared_error;
    target += - lambda * N * sum(sqrt(SS));
    target += sum(log(sigma));
}

