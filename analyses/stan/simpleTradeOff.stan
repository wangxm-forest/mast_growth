# May 18, 2026
# started by Mao, a simpler trade-off model linking seed and growth directly, the seed production is directly linked to tree size, growth, and climate variables.

data {
  int<lower=1> N;
  int<lower=0> sc[N]; 
  real<lower=0> rw[N];
  real<lower=0> DBH[N];
  real GST[N];
}

parameters {
  real alpha;
  real beta_rw;
  real beta_dbh;
  real beta_GST;
  real<lower=0> phi_sc;
}

model {
  alpha ~ normal(0, 5);
  beta_rw ~ normal(0, 5); // If beta_rw is negative, there's trade-off
  beta_dbh ~ normal(0, 5);
  beta_GST ~ normal(0, 5);
  phi_sc ~ gamma(2, 0.1);

  for (n in 1:N) {
    real log_mu = alpha + beta_rw * rw[n] + beta_dbh * DBH[n] + beta_GST * GST[n];
    
    sc[n] ~ neg_binomial_2_log(log_mu, phi_sc);
  }
}