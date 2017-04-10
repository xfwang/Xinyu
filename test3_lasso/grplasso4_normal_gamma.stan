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
    real<lower=0> tau_sq[G];
    real mu;
    real<lower=0> sigma;
    }
    
transformed parameters {
    vector[G] SS;
    real<lower=0> lambda_sq;
    lambda_sq = lambda ^2;
    for(i in 1:G) 
      SS[i] = sqrt(dot_self(beta[((i-1)*K+1) : (i*K)]));
    }

model {
    y ~ normal(x * beta, sigma);
    sigma ~ normal(0,1);
    for(i in 1:G) 
      beta[((i-1)*K+1) : (i*K)] ~ normal(mu, sigma*sqrt(tau_sq[i]));
    tau_sq ~ exponential(lambda_sq * N^2 / 8);
    target += - lambda * N * sum(SS); 
}


    