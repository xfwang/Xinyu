## Sparse group lasso, Laplace prior, lambda as a PARAMETER

data {
    int<lower=0> N; // number of samples
    int<lower=1> K; // number of covariates in each group 
    int<lower=1> G; // number of groups
    vector[N] y;
    matrix[N,K*G] x;
    }
    
parameters {
    real<lower=0, upper=1> p;
    real<lower=0> lambda;
    vector[K*G] beta;
    }
    
transformed parameters {
    real<lower=0> squared_error;
    vector[G] SS;

    squared_error = dot_self(y - x * beta);
    for(i in 1:G) 
      SS[i] = dot_self(beta[((i-1)*K+1) : (i*K)]);
    }

model {
    beta ~ double_exponential(0,1);
    target += -squared_error;
    target += - p * lambda * sqrt(K) * sum(sqrt(SS)); 
    for (k in 1:K)
      target += - (1-p) * lambda * fabs(beta[k]); 
    }

