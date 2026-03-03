data {
  int<lower=1> N;
  int<lower=1> NSite;
  int<lower=1, upper=NSite> site[N];
  int<lower=0> sc[N]; 
  real<lower=0> rw[N];
  real<lower=0> DBH[N];
  real GST[N];
  real GSP[N];
}

parameters {
  // Carbon availability
  real alpha;
  vector[NSite] alpha_site;
  real beta_dbh;
  real beta_GST;
  real beta_GSP;

  // Reproduction
  real<lower=0> gamma;
  vector<lower=0>[N] R;
  
  // Growth
  real<lower=0> beta_growth1;
  real<lower=0> beta_growth2;

  real<lower=0> phi;
  real<lower=0> sigma;
}

transformed parameters {
  vector[N] C;
  vector[N] G;
  vector[N] sc_pred;
  vector[N] rw_pred;

  for (n in 1:N) {
    C[n] = alpha + alpha_site[site[n]] + beta_dbh * DBH[n]
           + beta_GST * GST[n] + beta_GSP * GSP[n];

    sc_pred[n] = R[n] / gamma;

    G[n] = C[n] - R[n];

    rw_pred[n] = pow(G[n] / (0.5 * beta_growth1), 1 / beta_growth2);
  }
}

model {
  // Priors (adjust as appropriate)
  alpha ~ normal(0, 10);
  alpha_site ~ normal(0, 5);
  beta_dbh ~ normal(0, 5);
  beta_GST ~ normal(0, 5);
  beta_GSP ~ normal(0, 5);
  beta_growth1 ~ normal(0, 10);
  beta_growth2 ~ normal(0, 5);
  gamma ~ normal(0, 5);
  phi ~ normal(0, 5);
  sigma ~ normal(0, 5);
  
  R ~ uniform(0, C);

  // Likelihoods
  for (n in 1:N) {
    sc[n] ~ neg_binomial_2(sc_pred[n], phi);
    rw[n] ~ lognormal(log(rw_pred[n]), sigma);
  }
}

generated quantities {
  vector[N] sc_rep;
  vector[N] rw_rep;
  vector[N] frac_R;

  for (n in 1:N) {

    sc_rep[n] = neg_binomial_2_rng(R[n] / gamma, phi);

    rw_rep[n] = lognormal_rng(log(pow(G[n] / (0.5 * beta_growth1), 1 / beta_growth2)), sigma);
  }
}
