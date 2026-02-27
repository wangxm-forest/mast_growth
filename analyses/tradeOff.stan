data {
  int<lower=1> N; // Total observations (i * t)
  int<lower=1> N_sp; // Number of species
  int<lower=1> N_site; // Number of sites
  
  int<lower=1, upper=N_sp> sp[N];
  int<lower=1, upper=N_site> site[N];
  vector[N] sc;
  vector[N] rw;
  vector[N] dbh;
  vector[N_site] climate;

}

parameters {
  vector[N_sp] gamma_sp;
  vector[N_sp] alpha_sp;
  vector[N_sp] beta1_sp;           // Effect of DBH
  vector[N_sp] beta2_sp;           // Species-specific response to climate
  vector<lower=0>[N_sp] beta3_sp;
  vector[N_sp] beta4_sp;
  
  vector[N_site] alpha_site;       // Site random effect
  
  real<lower=0> sigma;             // Residual error
}

model {
  // Priors
  gamma_sp ~ normal(0, 1);
  alpha_sp ~ normal(0, 5);
  beta1_sp ~ normal(0, 1);
  beta2_sp ~ normal(0, 1);
  beta3_sp ~ lognormal(0, 1);
  beta4_sp ~ normal(0, 1);
  alpha_site ~ normal(0, 1);
  sigma ~ exponential(1);

for (n in 1:N) {
    real C_it = alpha_sp[sp[n]] + 
                    alpha_site[site[n]] + 
                    beta1_sp[sp[n]] * dbh[n] + 
                    beta2_sp[sp[n]] * climate_site[site[n]];
    
    real R_it = gamma_sp[sp[n]] * sc[n];
    real G_it = beta3_sp[sp[n]] * pow(rw[n], beta4_sp[sp[n]]) * 0.5;
    

    C_it ~ normal(R_it + G_it, sigma);
  }
}