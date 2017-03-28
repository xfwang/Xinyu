## Student-t prior plus

/* lg_t.stan */

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
  int<lower=0> N;
  int<lower=1> K;
  vector[N] y;
  matrix[N,K] x;
	real<lower=1> nu; // degrees of freedom for the half t-priors
}

parameters {

	// noise std
	real<lower=0> sigma;
  
	// auxiliary variables for the variance parameters
	vector[K] z;
	real<lower=0> r1_global;
	real<lower=0> r2_global;
	vector<lower=0>[K] r1_local;
	vector<lower=0>[K] r2_local;
	vector<lower=0>[K] r1_localplus;
	vector<lower=0>[K] r2_localplus;
}

transformed parameters {
	
	// global and local variance parameters, and the input weights
	real<lower=0> tau;
	vector<lower=0>[K] lambda;
	vector<lower=0>[K] lambdaplus;
	vector[K] beta;
	
	tau = r1_global * sqrt(r2_global);
	lambda = r1_local .* sqrt(r2_local);
	lambdaplus = r1_localplus .* sqrt_vec(r2_localplus);
	beta = z .* lambda .*lambdaplus*tau;
}

model {
	
	// observation model
	y ~ normal(x*beta, sigma);
	
	// half t-priors for lambdas (nu = 1 corresponds to horseshoe)
	z ~ normal(0, 1);
	r1_local ~ normal(0.0, 1.0);
	r2_local ~ inv_gamma(0.5*nu, 0.5*nu);
	r1_localplus ~ normal(0.0, 1.0);
	r2_localplus ~ inv_gamma(0.5*nu, 0.5*nu);
	
	// half cauchy for tau
	r1_global ~ normal(0.0, 1.0);
	r2_global ~ inv_gamma(0.5, 0.5);
	

	// using uniform prior on the noise variance
}
