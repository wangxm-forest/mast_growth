# May 18, 2026
# started by Mao, a simpler trade-off model linking seed and growth directly, the seed production is directly linked to tree size, growth, and climate variables.

data {
  int<lower=1> N;
  vector[N] BAI;
  int<lower=1> year[N];
  int<lower=1> site[N];
  
  int<lower=1> N_years;
  int<lower=1> N_sites;
  
  matrix[N_sites, N_years] Temp;
  matrix[N_sites, N_years] Prep;

  int<lower=1> N_sc;
  vector[N_sc] sc;
  int<lower=2> year_sc[N_sc];
  int<lower=1> site_sc[N_sc];
}

parameters {
  // seed
  real alpha_sc;
  real<lower=0> sigma_sc;
  real beta_sc_temp;
  
  //growth
  real alpha_BAI;
  real<lower=0> sigma_BAI;
  real beta_BAI_temp;
  real beta_BAI_prep;
  
  
  //trade-off
  real gamma_current;
  real gamma_lag;
  
}

transformed parameters {
  vector[N] G;

  matrix[N_sites, N_years] Gbar;
  matrix[N_sites, N_years] sum_G = rep_matrix(0, N_sites, N_years);
  matrix[N_sites, N_years] cnt   = rep_matrix(0, N_sites, N_years);
    for (n in 1:N) {
      G[n] = BAI[n] - alpha_BAI - beta_BAI_temp * Temp[site[n], year[n]] - beta_BAI_prep * Prep[site[n], year[n]];
      
      sum_G[site[n], year[n]] += G[n];
      cnt[site[n], year[n]]   += 1;
    }
    
    for (s in 1:N_sites){
      for (y in 1:N_years) {
    Gbar[s, y] = sum_G[s, y] / cnt[s, y];
      }
    }
}

model {
  alpha_BAI ~ normal(7, 3);
  sigma_BAI ~ normal(0, 1);
  beta_BAI_temp ~ normal(0, 1);
  beta_BAI_prep ~ normal(0, 1);

  alpha_sc ~ normal(0, 1);
  gamma_current ~ normal(0, 1);
  gamma_lag ~ normal(0, 1);
  sigma_sc ~ normal(0, 1);
  beta_sc_temp ~ normal(0, 1);

  for (n in 1:N){
  BAI[n] ~ normal(alpha_BAI + beta_BAI_temp * Temp[site[n], year[n]] + beta_BAI_temp * Prep[site[n], year[n]], sigma_BAI);
  }

  vector[N_sc] mu_sc;
  
  for (n in 1:N_sc){
    mu_sc[n] = alpha_sc
                 + gamma_current * Gbar[site_sc[n], year_sc[n]]
                 + gamma_lag * Gbar[site_sc[n], year_sc[n] - 1]
                 + beta_sc_temp * Temp[site_sc[n], year_sc[n]];
                 }
 
  sc ~ lognormal(mu_sc, sigma_sc);
}
