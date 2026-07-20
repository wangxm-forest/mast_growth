# May 18, 2026
# started by Mao, a simpler trade-off model linking seed and growth directly, the seed production is directly linked to tree size, growth, and climate variables.

data {
  int<lower=1> N;
  vector[N] BAI;
  int<lower=1> year[N];
  int<lower=1> N_years;

  int<lower=1> N_sc;
  vector[N_sc] sc;
  int<lower=2> year_sc[N_sc];
}

parameters {
  // seed
  real alpha_sc;
  real<lower=0> sigma_sc;
  
  //growth
  vector[N_years] G;
  real<lower=0> sigma_BAI;
  
  
  //trade-off
  real gamma_current;
  real gamma_lag;
  
}

model {
  G ~ normal(0, 5);
  sigma_BAI ~ exponential(1);

  alpha_sc ~ normal(0, 2);
  gamma_current ~ normal(0, 1);
  gamma_lag ~ normal(0, 1);
  sigma_sc ~ exponential(1);

  BAI ~ lognormal(G[year], sigma_BAI);

  vector[N_sc] log_mu_sc;
  
  for (n in 1:N_sc){
    log_mu_sc[n] = alpha_sc
                 + gamma_current * G[year_sc[n]]
                 + gamma_lag * G[year_sc[n] - 1];}

  sc ~ lognormal(log_mu_sc, sigma_sc);
}
