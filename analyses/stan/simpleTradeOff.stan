# May 18, 2026
# started by Mao, a simpler trade-off model linking seed and growth directly, the seed production is directly linked to tree size, growth, and climate variables.

data {
  int<lower=1> I;
  int<lower=2> T;
  real<lower=0> sc[I,T]; 
  real<lower=0> BAI[I,T];
  vector[T] GST;
}

parameters {
  // seed
  real alpha_sc;
  real beta_GST1;
  real<lower=0> sigma_sc;
  
  //growth
  real alpha_BAI;
  real beta_GST2;
  real<lower=0> sigma_BAI;
  matrix[I,T] G; 
  real<lower=0> sigma_G;
  
  //trade-off
  real gamma_current;
  real gamma_lag;
  
  
}

model {
  alpha_BAI ~ normal(5, 2);
  beta_GST1 ~ normal(0, 1);
  sigma_BAI ~ normal(0, 5);
  
  alpha_sc ~ normal(0.5, 1);
  beta_GST2 ~ normal(0, 1);
  sigma_sc ~ normal(0, 1);
  
  gamma_current ~ normal(0, 1);
  gamma_lag ~ normal(0, 1);

  for (i in 1:I) {

    real G_mu_1 = alpha_BAI + beta_GST2 * (GST[1]-15);

        G[i,1] ~ normal(G_mu_1, sigma_G);

        BAI[i,1] ~ lognormal(G[i,1], sigma_BAI);
    
    real log_mu_sc_1 = alpha_sc + beta_GST1 * (GST[1] - 15) + gamma_current * G[i, 1];
    
      sc[i, 1] ~ lognormal(log_mu_sc_1, sigma_sc);
  
    for (t in 2:T){
      
    real G_mu = alpha_BAI + beta_GST2 * (GST[t]-15);

        G[i,t] ~ normal(G_mu, sigma_G);

        BAI[i,t] ~ lognormal(G[i,t], sigma_BAI);
        
    real log_mu_sc = alpha_sc + beta_GST1 * (GST[t]-15) + gamma_current * G[i, t] + gamma_lag * G[i, t-1];
                       
      sc[i, t] ~ lognormal(log_mu_sc, sigma_sc);

  }
}
}


