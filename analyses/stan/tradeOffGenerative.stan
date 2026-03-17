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

      real rw_pred = (
        pow(G_n / beta_growth1 + pow(DBH[n], beta_growth2),
            1.0 / beta_growth2)
        - DBH[n]
      ) / 2.0;

      target += normal_lpdf(rw[n] | rw_pred, sigma_rw);
    }

    // reproducton
    target += neg_binomial_2_lpmf(
      sc[n] | gamma[n] * C[n], phi_sc
    );
  }
}

#generated quantities {

#  real C_pred;
#  real gamma_pred;
#  real G_pred;
#  real R_pred;

#  real rw_pred;
#  int sc_pred;

  // carbon
#  C_pred = lognormal_rng(alpha + beta_dbh * DBH[N] + beta_GST * GST[N],sigma_c);

#  gamma_pred = beta_rng(2,2);
  // allocation
#  G_pred = (1 - gamma_pred) * C_pred;
#  R_pred = gamma_pred * C_pred;
  // growth
#  rw_pred = normal_rng(theta * G_pred, sigma_rw);

#  sc_pred = poisson_rng(eta * R_pred);

#}
