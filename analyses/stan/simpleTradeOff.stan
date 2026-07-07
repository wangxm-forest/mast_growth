# May 18, 2026
# started by Mao, a simpler trade-off model linking seed and growth directly, the seed production is directly linked to tree size, growth, and climate variables.

data {
  int<lower=1> N;
  vector[N] BAI;
  vector[N] BAI_lag; 
  int<lower=1> year[N];
  int<lower=1> N_years;
  
  int<lower=1> N_t;
  int<lower=0> sc[N_t];
  int<lower=0> sc_lag[N_t];
  int<lower=1> year_t[N_t];
  int<lower=1> year_t_lag[N_t];
}

parameters {
  // seed
  real alpha_sc;
  real beta_sc;
  real<lower=0> phi_sc;
  
  //growth
  real alpha_BAI;
  real beta_BAI;
  real<lower=0> sigma_BAI;
  
  
  //trade-off
  real gamma1;
  real gamma2;
  
}

model {
  alpha_BAI    ~ normal(0, 1);
  beta_BAI    ~ normal(0, 1);
  sigma_BAI ~ normal(0, 1);

  alpha_sc ~ normal(0, 2);
  beta_sc  ~ normal(0, 1);
  gamma1   ~ normal(0, 1);
  gamma2   ~ normal(0, 1);
  phi_sc   ~ exponential(1);

BAI ~ normal(alpha_BAI + beta_BAI * BAI_lag, sigma_BAI);

  vector[N] G = alpha_BAI + beta_BAI * BAI_lag;
  vector[N_years] sum_G = rep_vector(0, N_years);
  vector[N_years] cnt   = rep_vector(0, N_years);
  for (n in 1:N) {
    sum_G[year[n]] += G[n];
    cnt[year[n]]   += 1;
  }
  vector[N_years] Gbar = sum_G ./ cnt;
  
  vector[N_t] log_sc_mu;
  for (t in 1:N_t)
    log_sc_mu[t] = alpha_sc
              + gamma1 * Gbar[year_t_lag[t]]
              + gamma2 * Gbar[year_t[t]];

  sc ~ neg_binomial_2_log(log_sc_mu, phi_sc);
  
  }


