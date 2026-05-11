data {
  int<lower=1> N;
  //int<lower=1> NSite;
  //int<lower=1, upper=NSite> site[N];
  int<lower=0,upper=1> mast[N];
  int<lower=0> sc[N]; 
  real<lower=0> rw[N];
  real<lower=0> DBH[N];
  real GST[N];
  
  //Fix intercept for total carbon
  real<lower=0> alpha;
  //Fix growth parameters
  real<lower=0> beta_growth1;
  real<lower=0> beta_growth2;
  
}

parameters {
  // carbon
  // real<lower=0> C[N];
  // real alpha;
  // vector[NSite] alpha_site;
  real beta_dbh;
  real beta_GST;

  // allocation
  real<lower=0,upper=1> gamma_mast;
  real<lower=0,upper=1> gamma_nonmast;
  //real<lower=0, upper=1> mu_gamma; 
  //real<lower=0.1> kappa_gamma;
  //vector<lower=0,upper=1>[N] gamma;
  // growth
  //real<lower=0> beta_growth1;
  //real<lower=0> beta_growth2;
  real<lower=0> sigma_rw;

  // reproduction
  real<lower=0> phi_sc;
  real<lower=0> theta;

}

transformed parameters {
  vector[N] C;

// Version with the site specific intercept
  //for (n in 1:N) {
   // C[n] = alpha + alpha_site[site[n]] + beta_dbh * DBH[n]
         //  + beta_GST * GST[n];
           
  for (n in 1:N) {
    C[n] = alpha + beta_dbh * DBH[n]
           + beta_GST * GST[n];
  }
}
model {
// Version with undefined intercept
  //alpha ~ lognormal(3,0.5);
  //alpha_site ~ normal(0, 10);

  beta_dbh ~ normal(0, 1);
  beta_GST ~ normal(0, 1);

// Version adding variance to the total carbon
//  sigma_c ~ normal(0, 1);

// Version hierchical gamma
  gamma_mast ~ beta(2,2);
  gamma_nonmast ~ beta(2,2);
//mu_gamma ~ beta(2, 2); 
//kappa_gamma ~ exponential(0.1); 
//gamma ~ beta(mu_gamma * kappa_gamma, (1 - mu_gamma) * kappa_gamma);

//putting a very strong priors on allometric parameters
  //beta_growth1 ~ normal(3, 0.5);
  //beta_growth2 ~ normal(3, 0.5);

  sigma_rw ~ normal(0, 1);

  phi_sc ~ gamma(2, 0.1);
  theta ~ normal(0, 1);

  for (n in 1:N) {
    real gamma_n;

if (mast[n] == 1)
  gamma_n = gamma_mast;
else
  gamma_n = gamma_nonmast;

    // carbon
//    target += lognormal_lpdf(
//      C[n] |
//      log(alpha
//      + alpha_site[site[n]]
//      + beta_dbh * DBH[n]
//      + beta_GST * GST[n]),
//      sigma_c
//    );

    // growth
    {
      real G_n = (1 - gamma_n) * C[n];

      target += lognormal_lpdf(rw[n] | log((pow(G_n / (beta_growth1/100) + pow(DBH[n], beta_growth2),
            1 / beta_growth2)
        - DBH[n]
      ) / 2), sigma_rw);
    }

    // reproducton
    {
      real R_n = gamma_n * C[n];
      
    target += neg_binomial_2_lpmf(
      sc[n] | R_n / (theta/100), phi_sc);
    }
  }
}


