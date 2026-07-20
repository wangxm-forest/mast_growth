# May 18, 2026
# started by Mao, a simpler trade-off model linking seed and growth directly, the seed production is directly linked to tree size, growth, and climate variables.

data {
  int<lower=1> N;
  vector[N] BAI;
  vector[N] BAI_lag; 
  int<lower=1> year[N];
  int<lower=1> N_years;
  
  int<lower=1> N_t;
  vector[N_t] repro;
  vector[N_t] repro_lag;
  int<lower=1> year_t[N_t];
  int<lower=1> year_t_lag[N_t];
}

parameters {
  real alpha1;
  real<lower=0> sigma_BAI;

  real beta;
  real gamma1;
  real gamma2;
  real<lower=0> sigma_sc;
  
}

model {
  alpha1    ~ normal(0, 5);
  sigma_BAI ~ normal(0, 5);
  beta      ~ normal(0, 5);
  gamma1    ~ normal(0, 1);
  gamma2    ~ normal(0, 5);
  sigma_sc  ~ normal(0, 5);
  
  BAI ~ normal(alpha1 * BAI_lag, sigma_BAI);

  vector[N] G = alpha1 * BAI_lag;
  vector[N_years] sum_G = rep_vector(0, N_years);
  vector[N_years] cnt   = rep_vector(0, N_years);
  for (n in 1:N) {
    sum_G[year[n]] += G[n];
    cnt[year[n]]   += 1;
  }
  vector[N_years] Gbar = sum_G ./ cnt;


  vector[N_t] R;
  for (t in 1:N_t){
    R[t] = beta * repro_lag[t]
         + gamma1 * Gbar[year_t_lag[t]]
         + gamma2 * Gbar[year_t[t]];
  }
  repro ~ normal(R, sigma_sc);
  
  

}
