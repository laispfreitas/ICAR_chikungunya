data {
  int<lower=0> N;                         // number of areas
  int<lower=1> T;                         // number of times
  int<lower=0> N_edges;                   // number of edges
  int<lower=1, upper=N> node1[N_edges];   // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];   // and node1[i] < node2[i]
  int<lower=0> y[N,T];                    // count outcomes of chikungunya
  vector<lower=0>[N] E;                   // expected number of cases
  real m0;                                
  real<lower=0> C0;
  
}

transformed data {
  
  vector[N] log_E = log(E);
  
}


parameters {
  
  real beta0;                        // intercept
  real<lower=0> sigma;                // standard deviation of the spatial random effect
  vector[N] phi;                      // spatial effects

  
  
}


transformed parameters{
  
  matrix<lower=0> [N,T] mu;   // poisson parameter

  for(i in 1:N){
      mu[i,1]=exp(log_E[i]+beta0+sigma*phi[i]);
    for(t in 2:T){
        mu[i,t]=exp(log_E[i]+beta0+sigma*phi[i]);
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

  sigma~gamma(1,1);

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
