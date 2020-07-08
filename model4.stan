data {
  int<lower=0> N;                         // number of areas
  int<lower=1> T;                         // number of times
  int<lower=0> N_edges;                   // number of edges
  int<lower=1, upper=N> node1[N_edges];   // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];   // and node1[i] < node2[i]
  int<lower=0> y[N,T];                    // count outcomes of chikungunya
  vector<lower=0>[N] E;                   // expected number of cases
  int<lower=1> K;                         // num covariates
  matrix[N, K] x;                         // design matrix
  matrix[N, T] tmincenter;                // standardized max temperature
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
  vector[N] zeta;                     // instantaneous effect of the temperature per neighbourhood
  real zetacent;                      // centered instantaneous effect of the temperature
  real <lower=0> sigmacent;           // standard deviation of the centered zeta
  vector<lower=0,upper=1> [N] rho;    // memory effect of the temperature
  real<lower=0> sigma;                // standard deviation of the spatial random effect
  vector[N] phi;                      // spatial effects

  
  
}


transformed parameters{
  
  matrix<lower=0> [N,T] mu;   // poisson parameter
  matrix[N,T] U;              // transfer function

  U[,1]=tmincenter[,1].*zeta; 
  for(i in 1:N){
      mu[i,1]=exp(log_E[i]+beta0+x[i,]*betas[,1]+U[i,1]+sigma*phi[i]);
    for(t in 2:T){
        U[i,t]=U[i,(t-1)]*rho[i]+tmincenter[i,t]*zeta[i];
        mu[i,t]=exp(log_E[i]+beta0+x[i,]*betas[,t]+U[i,t]+sigma*phi[i]);
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
 
  // This is the prior for phi! (up to proportionality)
  target += -0.5 * dot_self(phi[node1] - phi[node2]);

  rho~uniform(0,1);
  betas[,1]~normal(0.0,5);
  for(t in 2:T){
   betas[,t]~normal(betas[,t-1],W);
  }  
  W~normal(0,1);
  sigma~gamma(1,1);
  zeta~normal(zetacent,sigmacent);
  zetacent~normal(0,5);
  sigmacent~normal(0,1);  

// soft sum-to-zero constraint on phi
  sum(phi) ~ normal(0, 0.001 * N);  // equivalent to mean(phi) ~ normal(0,0.001)
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
