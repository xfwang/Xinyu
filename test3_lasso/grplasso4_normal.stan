## group lasso, Gasussian prior, lambda as a PARAMETER, with hyperparameter

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
    real<lower=0> gamma;
    real<lower=0> sigma;
    }
    
transformed parameters {
    vector[G] SS;

    for(i in 1:G) 
      SS[i] = sqrt(dot_self(beta[((i-1)*K+1) : (i*K)]));
    }

model {
    y ~ normal(x * beta, sigma);
    for(j in 1:(K*G)) 
      beta[j] ~ normal(mu, gamma);
    sigma ~ normal(0,1);
    target += - lambda * N * sum(SS); 
}

