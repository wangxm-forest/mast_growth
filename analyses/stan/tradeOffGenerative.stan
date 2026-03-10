data {
  int<lower=1> N;
  real<lower=0> rw[N];
  int<lower=0> sc[N];
  real<lower=0> DBH[N];
  real GST[N];
}

parameters {
  // carbon
  vector<lower=0>[N] C;
  vector<lower=0,upper=0.5>[N] gamma;
  real alpha;
  real beta_dbh;
  real beta_GST;
  real<lower=0> sigma_c;

  // growth
  real theta;
  real<lower=0> sigma_rw;

  // reproduction
  real eta;

}

model {
  alpha ~ lognormal(log(90),0.5);
  beta_dbh ~ normal(0, 5);
  beta_GST ~ normal(0, 5);
  sigma_c ~ normal(0, 5);
  sigma_rw ~ normal(0, 0.5);
  theta ~ normal(0,0.05);
  eta ~ normal(0,50);
  
  gamma ~ beta(2,2);

  for (n in 1:N) {

    // cabon availability
    C[n] ~ normal(alpha + beta_dbh * DBH[n] + beta_GST * GST[n], sigma_c);

    // carbon allocations
    real G;
    real R;

    G = (1 - gamma[n]) * C[n];
    R = gamma[n] * C[n];

    // growth
    rw[n] ~ normal(theta * G, sigma_rw);

    // reproduction
    sc[n] ~ poisson(eta * R);

  }

}

generated quantities {

  real C_pred;
  real gamma_pred;
  real G_pred;
  real R_pred;

  real rw_pred;
  int sc_pred;

  // carbon
  C_pred = lognormal_rng(alpha + beta_dbh * DBH[N] + beta_GST * GST[N],sigma_c);

  gamma_pred = beta_rng(2,2);
  // allocation
  G_pred = (1 - gamma_pred) * C_pred;
  R_pred = gamma_pred * C_pred;
  // growth
  rw_pred = normal_rng(theta * G_pred, sigma_rw);

  sc_pred = poisson_rng(eta * R_pred);

}
