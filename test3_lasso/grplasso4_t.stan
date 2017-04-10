## group lasso, Student-t prior 
functions {
	// square root of a vector (elementwise)
	vector sqrt_vec(vector x) {
		vector[dims(x)[1]] res;

		for (m in 1:dims(x)[1]){
			res[m] = sqrt(x[m]);
		}
		return res;
	}
}

data {
    int<lower=0> N; // number of samples
    int<lower=1> K; // number of covariates in each group 
    int<lower=1> G; // number of groups
    vector[N] y;
    matrix[N,K*G] x;
	real<lower=1> nu; // degrees of freedom for the half t-priors
    }
    
parameters {
    real<lower=0> sigma;
	vector[K*G] z;
	real<lower=0> r1_global;
	real<lower=0> r2_global;
	vector<lower=0>[K*G] r1_local;
	vector<lower=0>[K*G] r2_local;
    }
    
transformed parameters {
    vector[G] SS;
	real<lower=0> tau;
	vector<lower=0>[K*G] lambda;
    vector[K*G] beta;
	
    for(i in 1:G) 
      SS[i] = sqrt(dot_self(beta[((i-1)*K+1) : (i*K)]));
	tau = r1_global * sqrt(r2_global);
	lambda = r1_local .* sqrt(r2_local);
	beta = z .* lambda*tau;
    }

model {
    y ~ normal(x * beta, sigma);
	z ~ normal(0, 1);
	r1_local ~ normal(0.0, 1.0);
	r2_local ~ inv_gamma(0.5*nu, 0.5*nu);
	
	// half cauchy for tau
	r1_global ~ normal(0.0, 1.0);
	r2_global ~ inv_gamma(0.5, 0.5);
    target += - lambda * N * sum(SS); 
}


