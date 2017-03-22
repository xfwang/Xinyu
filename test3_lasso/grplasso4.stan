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
    real mu;
    real<lower=0> sigma;
    }
    
transformed parameters {
    real<lower=0> squared_error;
    vector[G] SS;

    squared_error = dot_self(y - x * beta);
    for(i in 1:G) 
      SS[i] = dot_self(beta[((i-1)*K+1) : (i*K)]);
    }

model {
    for(j in 1:(K*G)) 
      beta[j] ~ double_exponential(mu, sigma);
    target += -squared_error;
    target += - lambda * sum(sqrt(SS)); 
}

