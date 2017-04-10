## group lasso, without any prior

data {
    int<lower=0> N; // number of samples
    int<lower=1> K; // number of covariates in each group 
    int<lower=1> G; // number of groups
    vector[N] y;
    matrix[N,K*G] x;
    real<lower=0> lambda;
    }
    
parameters {
    vector[K*G] beta;
  real<lower=0> sigma;
    }
    
transformed parameters {
    vector[G] SS;

    for(i in 1:G) 
      SS[i] = sqrt(dot_self(beta[((i-1)*K+1) : (i*K)]));
    }

model {
  // beta ~ double_exponential(0,1);
  y ~ normal(x * beta, sigma);
  sigma ~ normal(0,1);
  target += - lambda * sum(SS); 
}
