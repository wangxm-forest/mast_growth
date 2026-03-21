data {
  int<lower=1> N;
  int<lower=1> NSite;
  int<lower=1, upper=NSite> site[N];
  int<lower=0> sc[N]; 
  real<lower=0> rw[N];
  real<lower=0> DBH[N];
  real GST[N];
}

parameters {
  // carbon
  real<lower=0> C[N];
  real alpha;
  vector[NSite] alpha_site;
  real beta_dbh;
  real beta_GST;
  real<lower=0> sigma_c;

  // allocation
  vector<lower=0,upper=1>[N] gamma;
  // growth
  real<lower=0> beta_growth1;
  real<lower=0> beta_growth2;
  real<lower=0> sigma_rw;

  // reproduction
  real<lower=0> phi_sc;
  real<lower=0> theta;

}


model {
  alpha ~ lognormal(log(90),0.5);
  alpha_site ~ normal(0, 10);

  beta_dbh ~ normal(0, 1);
  beta_GST ~ normal(0, 5);

  sigma_c ~ normal(0, 1);

  gamma ~ normal(0, 1);

  beta_growth1 ~ lognormal(1, 0.5);
  beta_growth2 ~ lognormal(1, 0.5);

  sigma_rw ~ normal(0, 1);

  phi_sc ~ gamma(2, 0.1);
  theta ~ normal(0, 5);

  for (n in 1:N) {

    // carbon
    target += normal_lpdf(
      C[n] |
      alpha
      + alpha_site[site[n]]
      + beta_dbh * DBH[n]
      + beta_GST * GST[n],
      sigma_c
    );

    // growth
    {
      real G_n = (1 - gamma[n]) * C[n];

      target += normal_lpdf(rw[n] | (pow(G_n / beta_growth1 + pow(DBH[n], beta_growth2),
            1 / beta_growth2)
        - DBH[n]
      ) / 2, sigma_rw);
    }

    // reproducton
    {
      real R_n = gamma[n] * C[n];
      
    target += neg_binomial_2_lpmf(
      sc[n] | R_n / theta, phi_sc);
    }
  }
}

generated quantities {

  vector[N] C_pred;
  vector[N] gamma_pred;
  vector[N] G_pred;
  vector[N] R_pred;

  vector[N] rw_pred;
  array[N] int sc_pred;

  for (n in 1:N) {
  // carbon
  C_pred[n] = normal_rng(alpha + alpha_site[site[n]] + beta_dbh * DBH[n] + beta_GST * GST[n],sigma_c);

  gamma_pred[n] = gamma[n];
  // allocation
  G_pred[n] = (1 - gamma_pred[n]) * C_pred[n];
  R_pred[n] = gamma_pred[n] * C_pred[n];
  // growth
  rw_pred[n] = normal_rng((pow(G_pred[n] / beta_growth1 + pow(DBH[n], beta_growth2),
            1 / beta_growth2)
        - DBH[n]
      ) / 2, sigma_rw);
  // reproduction
  sc_pred[n] = neg_binomial_2_rng(R_pred[n] / theta, phi_sc);

  }
}
