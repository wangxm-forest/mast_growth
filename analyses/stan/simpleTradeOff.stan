# May 18, 2026
# started by Mao, a simpler trade-off model linking seed and growth directly, the seed production is directly linked to tree size, growth, and climate variables.

data {
  int<lower=1> I;
  int<lower=2> T;
  int<lower=0> sc[I,T]; 
  real<lower=0> BAI[I,T];
  vector[T] GST;
}

parameters {
  // seed
  real alpha_sc;
  real beta_GST1;
  real<lower=0> phi_sc;
  
  //growth
  real alpha_BAI;
  real beta_dbh;
  real beta_GST2;
  real<lower=0> sigma_BAI;
  
  //trade-off
  real gamma_current;
  real gamma_lag;
}

model {
  matrix[I,T] G;
  
  alpha_BAI ~ normal(0, 5);
  beta_GST1 ~ normal(0, 5);
  sigma_BAI ~ normal(0, 1);
  
  alpha_sc ~ normal(0, 5);
  beta_dbh ~ normal(0, 5);
  beta_GST2 ~ normal(0, 5);
  phi_sc ~ gamma(2, 0.1);
  
  gamma_current ~ normal(0, 5);
  gamma_lag ~ normal(0, 5);

  for (i in 1:I) {

  G[i, 1] = alpha_BAI + beta_GST2 * GST[1];
  BAI[i, 1] ~ lognormal(G[i, 1], sigma_BAI);
  
    for (t in 2:T){  
      G[i, t] = alpha_BAI + beta_GST2 * GST[t];
      BAI[i, t] ~ lognormal(G[i, t], sigma_BAI);
      
  real log_mu_sc = alpha_sc + beta_GST2 * GST[t] + gamma_current * G[i, t] + gamma_lag * G[i, t-1];
                       
      sc[i, t] ~ neg_binomial_2_log(log_mu_sc, phi_sc);

  }
}
}


