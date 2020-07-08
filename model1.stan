data {
  int<lower=0> N;                         // number of areas
  int<lower=1> T;                         // number of times
  int<lower=0> y[N,T];                    // count outcomes of chikungunya
  vector<lower=0>[N] E;                   // expected number of cases
  int<lower=1> K;                         // num covariates
  matrix[N, K] x;                         // design matrix
  real m0;                                
  real<lower=0> C0;
  
}

transformed data {
  
  vector[N] log_E = log(E);
  
}


parameters {
  
  real beta0;                        // intercept
  matrix[K,T] betas;                  // covariates' coefficients 
  real<lower=0> W;                    // standard deviation of the coefficients' random walk

}


transformed parameters{
  
  matrix<lower=0> [N,T] mu;   // poisson parameter

  for(i in 1:N){
      mu[i,1]=exp(log_E[i]+beta0+x[i,]*betas[,1]);
    for(t in 2:T){
        mu[i,t]=exp(log_E[i]+beta0+x[i,]*betas[,t]);
        }
  }
  
}


model {
  
  matrix[N,T] lps;
  
  beta0~normal(m0,C0);
  
  // likelihood function
  for(i in 1:N){
    for(t in 1:T){
      y[i,t] ~ poisson(mu[i,t]);
      }
}
 
  betas[,1]~normal(0.0,5);
  for(t in 2:T){
   betas[,t]~normal(betas[,t-1],W);
  }  
  W~normal(0,1);
  
}


generated quantities {
  
  int<lower=0> yfit[N,T];
  matrix[N,T] log_lik;

  
  for(i in 1:N){
    for(t in 1:T){
       yfit[i,t]=poisson_rng(mu[i,t]);
       log_lik[i,t]=poisson_lpmf(y[i,t] | mu[i,t]);
    }
  }
}


