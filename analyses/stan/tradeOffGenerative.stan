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
  theta ~ normal(0, 1);

  for (n in 1:N) {

    // carbon
    target += lognormal_lpdf(
      C[n] |
      log(alpha
      + alpha_site[site[n]]
      + beta_dbh * (DBH[n]-20)
      + beta_GST * GST[n]),
      sigma_c
    );

    // growth
    {
      real G_n = (1 - gamma[n]) * C[n];

      target += lognormal_lpdf(rw[n] | log((pow(G_n / (beta_growth1/100) + pow(DBH[n], beta_growth2),
            1 / beta_growth2)
        - DBH[n]
      ) / 2), sigma_rw);
    }

    // reproducton
    {
      real R_n = gamma[n] * C[n];
      
    target += neg_binomial_2_lpmf(
      sc[n] | R_n / (theta/100), phi_sc);
    }
  }
}


