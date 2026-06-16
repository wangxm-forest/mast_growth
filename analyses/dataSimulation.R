#### Growth-reproduction trade-off ####
## Started by Mao ##
## Feb-28-2026 ##

library(MASS) 
library(rstan)

rm(list = ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

setwd("C:/PhD/Project/PhD_thesis/mast_growth/analyses")
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

set.seed(123)
###for tradeoff model###
###Set values for parameters###
N <- 50
NSite <- 3

DBH <- runif(N, 20, 40)
GST <- runif(N, 15, 25)
#GSP <- runif(N, 50, 150)
site <- sample(1:NSite, N, replace = TRUE)

alpha <- 50
alpha_site <- c(0, 20, -10)
beta_dbh <- 0.5
beta_GST <- 0.5
#beta_GSP <- 0.3
beta_growth1 <- 4
beta_growth2 <- 3.2
theta <- 1
phi_sc <- 2
sigma_rw <- 0.5


C <- alpha + alpha_site[site] + beta_dbh * DBH + beta_GST * GST

frac_R <- runif(N, 0.1, 0.5)
R <- C * frac_R
G <- C - R 

sc_pred <- R / theta
sc <- rnegbin(N, mu = sc_pred, theta = phi_sc)

rw_pred <- (((G / beta_growth1 + DBH^beta_growth2)^(1 / beta_growth2)) - DBH) / 2
rw <- rlnorm(N, meanlog = log(rw_pred), sdlog = sigma_rw)

# Data
stanData <- list(
  N = N,
  NSite = NSite,
  site = site,
  sc = sc,
  rw = rw,
  DBH = DBH,
  GST = GST
)

fit <- stan(file='stan/tradeOff.stan', data=stanData, seed=112234, control=list(adapt_delta=0.99))
diagnostics <- util$extract_hmc_diagnostics(fit)

print(util$check_all_hmc_diagnostics(diagnostics))


samples <- util$extract_expectand_vals(fit)
names <- c(grep('alpha', names(samples), value = TRUE),
           grep('alpha_site', names(samples), value = TRUE),
           grep('beta_dbh', names(samples), value = TRUE),
           grep('beta_GST', names(samples), value = TRUE),
           grep('beta_growth1', names(samples), value = TRUE),
           grep('beta_growth2', names(samples), value = TRUE),
           grep('gamma', names(samples), value = TRUE),
           grep('phi', names(samples), value = TRUE),
           grep('sigma', names(samples), value = TRUE),
           grep('theta', names(samples), value = TRUE))

base_samples <- util$filter_expectands(samples,names)
print(util$check_all_expectand_diagnostics(base_samples))
pdf("pairsPlot.pdf", height = 9, width = 9)
par(mfrow = c(3,3))
util$plot_div_pairs(names, names, samples, diagnostics)
dev.off()

pdf("figures/priorPosteriorPlotBasic.pdf", height = 9, width = 9)
par(mfrow = c(3,3))
util$plot_expectand_pushforward(samples[['alpha']], 50, display_name = "alpha")
curve(dlnorm(x, log(90), 0.5),
      add = TRUE,
      col = "blue",
      lwd = 2)
abline(v = 50, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['alpha_site[1]']], 50, display_name = "alpha_site[1]")
curve(dnorm(x, 0, 10),
      add = TRUE,
      col = "blue",
      lwd = 2)
abline(v = 0, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['beta_dbh']], 50, display_name = "beta_dbh")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.5, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['beta_GST']], 50, display_name = "beta_GST")
curve(dnorm(x, 0, 5),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.5, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['beta_growth1']], 50, display_name = "beta_growth1")
curve(dlnorm(x, 1, 0.5),
      add = TRUE,
      col = "blue",
      lwd = 2)
abline(v = 0.04, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['beta_growth2']], 50, display_name = "beta_growth2", flim=c(1.5,3.5))
curve(dlnorm(x, 1, 0.5),
      add=TRUE, col="blue", lwd=2, lty=2)
abline(v = 3.2, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['theta']], 50, display_name = "theta")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.01, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['phi_sc']], 50, display_name = "phi_sc")
curve(dgamma(x, 2, 0.1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 2, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['sigma_rw']], 50, display_name = "sigma_rw")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.5, col = "red", lwd = 2)
dev.off()

rm(list = ls())
###for generative model###
N <- 1000
mast <- rbinom(N, 1, 0.5)
#NSite <- 3

DBH <- runif(N, 20, 40)
GST <- runif(N, 15, 25)
#site <- sample(1:NSite, N, replace = TRUE)

alpha <- 30
#alpha_site <- c(0, 20, -10)
beta_dbh <- 0.5
beta_GST <- 0.5
#sigma_c <- 0.1
#beta_GSP <- 0.3
beta_growth1 <- 3
beta_growth2 <- 3
#mu_gamma <- 0.5
#kappa_gamma <- 50
gamma_mast <- 0.8
gamma_nonmast <- 0.2
gamma <- ifelse(mast == 1,
                gamma_mast,
                gamma_nonmast)
sigma_rw <- 0.1
phi_sc <- 2
theta <- 1


  #Carb <- alpha + alpha_site[site] + beta_dbh * DBH + beta_GST * GST
  Carb <- alpha + beta_dbh * DBH + beta_GST * GST
  # carbon allocation
  G_n <- (1-gamma) * Carb
  R_n <- gamma * Carb
  
  # ring width
  rw <- rlnorm(N, meanlog = log(((G_n / (beta_growth1/100) + DBH^beta_growth2)^(1 / beta_growth2) - DBH) / 2), sdlog = sigma_rw)
  
  # seed counts
  sc <- rnegbin(N, mu = R_n/(theta/100), theta = phi_sc)

stanData <- list(
  N = N,
  mast = mast,
  alpha = alpha,
#  NSite = NSite,
#  site = site,
  sc = sc,
  rw = rw,
  DBH = DBH,
  GST = GST,
  beta_growth1 = beta_growth1,
  beta_growth2 = beta_growth2
)

mod <- stan_model(file='stan/tradeOffGenerative.stan')

fit <- stan(file='stan/tradeOffGenerative.stan', data=stanData, seed=112234, control=list(adapt_delta=0.99))
diagnostics <- util$extract_hmc_diagnostics(fit)

print(util$check_all_hmc_diagnostics(diagnostics))


samples <- util$extract_expectand_vals(fit)
names <- c(grep('alpha', names(samples), value = TRUE),
           #grep('alpha_site', names(samples), value = TRUE),
           grep('beta_dbh', names(samples), value = TRUE),
           grep('beta_GST', names(samples), value = TRUE),
           grep('beta_growth1', names(samples), value = TRUE),
           grep('beta_growth2', names(samples), value = TRUE),
           grep('gamma', names(samples), value = TRUE),
           grep('phi', names(samples), value = TRUE),
           grep('sigma', names(samples), value = TRUE),
           grep('theta', names(samples), value = TRUE))

base_samples <- util$filter_expectands(samples,names)
print(util$check_all_expectand_diagnostics(base_samples))
print(fit, pars = names)

pdf("pairsPlotGeneChangePriors.pdf", height = 9, width = 9)
par(mfrow = c(3,3))
util$plot_div_pairs(names, names, samples, diagnostics)
dev.off()

pdf("figures/priorPosteriorPlot.pdf", height = 9, width = 9)
par(mfrow = c(3,3))
#util$plot_expectand_pushforward(samples[['alpha']], 50, display_name = "alpha",flim = c(10,55))
#curve(dlnorm(x, 3,0.5),
#      add = TRUE,
#      col = "blue",
#      lwd = 2,
#      xlim = c(10,55))
#abline(v = 30, col = "red", lwd = 2)

#util$plot_expectand_pushforward(samples[['alpha_site[1]']], 50, display_name = "alpha_site[1]")
#curve(dnorm(x, 0, 10),
#      add = TRUE,
#      col = "blue",
#      lwd = 2)
#abline(v = 0, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['beta_dbh']], 50, display_name = "beta_dbh")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.5, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['beta_GST']], 50, display_name = "beta_GST")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.5, col = "red", lwd = 2)

#util$plot_expectand_pushforward(samples[['beta_growth1']], 50, display_name = "beta_growth1")
#curve(dlnorm(x, 3, 0.01),
#      add = TRUE,
#      col = "blue",
#      lwd = 2)
#abline(v = 3, col = "red", lwd = 2)

#util$plot_expectand_pushforward(samples[['beta_growth2']], 50, display_name = "beta_growth2", flim=c(1.5,3.5))
#curve(dlnorm(x, 3, 0.01),
#      add=TRUE, col="blue", lwd=2, lty=2)
#abline(v = 3, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['theta']], 50, display_name = "theta")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 1, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['phi_sc']], 50, display_name = "phi_sc")
curve(dgamma(x, 2, 0.1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 2, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['sigma_rw']], 50, display_name = "sigma_rw",flim = c(0,1))
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.1, col = "red", lwd = 2)

#util$plot_expectand_pushforward(samples[['sigma_c']], 50, display_name = "sigma_c")
#curve(dnorm(x, 0, 1),
#      add = TRUE,
#      col = "blue",
#      lwd = 2,
#      lty = 2)
#abline(v = 0.1, col = "red", lwd = 2)

#util$plot_expectand_pushforward(samples[['mu_gamma']], 50, display_name = "mu_gamma",flim = c(0.25,0.55))
#curve(dbeta(x, 2, 2),
#      add = TRUE,
#      col = "blue",
#      lwd = 2,
#      lty = 2)
#abline(v = 0.5, col = "red", lwd = 2)

#util$plot_expectand_pushforward(samples[['kappa_gamma']], 50, display_name = "kappa_gamma", flim = c(10,100))
#curve(dexp(x, 0.1),
#      add = TRUE,
#      col = "blue",
#      lwd = 2,
#      lty = 2)
#abline(v = 20, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['gamma_mast']], 50, display_name = "gamma_mast",flim = c(0,1))
curve(dbeta(x, 2, 2),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.8, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['gamma_nonmast']], 50, display_name = "gamma_nonmast",flim = c(0,1))
curve(dbeta(x, 2, 2),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.2, col = "red", lwd = 2)



dev.off()

###for simple model###
set.seed(2266)
I <- 10
T <- 30

# Simulate predictors
# Growing season temperature
GST <- rnorm(T, mean = 15, sd = 1)
# Define True Parameters
alpha_BAI <- 2
beta_GST2 <- 0.05
sigma_BAI <- 0.5

alpha_sc <- 0.5
beta_GST1 <- 0.3
gamma_current <- -0.2
gamma_lag <- -0.1

sigma_sc <- 1
sigma_G <- 0.5


BAI <- matrix(NA, nrow = I, ncol = T)
sc  <- matrix(0,  nrow = I, ncol = T)
G <- matrix(NA, nrow = I, ncol =T)

for (i in 1:I) {
  for (t in 1:T) {

    G_mu <- alpha_BAI + beta_GST2 * (GST[t]-15)
    
    G[i, t] <- rnorm(1, G_mu, sigma_G)

    BAI[i, t] <- rlnorm(1, meanlog = G[i, t], sdlog = sigma_BAI)
    
    if (t > 1) {
      log_mu_sc <- alpha_sc + 
        beta_GST1 * (GST[t]-15) + 
        gamma_current * G[i, t] + 
        gamma_lag * G[i, t-1]
      
      sc[i, t] <- rlnorm(1, meanlog = log_mu_sc, sdlog = sigma_sc)
      
    } else {
      G_mu_1 <- alpha_BAI + beta_GST2 * (GST[1]-15)
      
      G[i, 1] <- rnorm(1, G_mu_1, sigma_G)
      
      BAI[i, 1] <- rlnorm(1, meanlog = G[i, 1], sdlog = sigma_BAI)
      
      log_mu_sc_year1 <- alpha_sc + beta_GST1 * (GST[t]-15) + gamma_current * log(BAI[i, 1])
      
      sc[i, 1] <- rlnorm(1, meanlog = log_mu_sc_year1, sdlog = sigma_sc)
    }
  }
}

stanData <- list(
  I = I,
  T = T,
  sc = sc,
  BAI = BAI,
  GST = GST
)

mod <- stan_model(file='stan/simpleTradeOff.stan')

fit <- stan(file='stan/simpleTradeOff.stan', data=stanData, seed=112234, control=list(adapt_delta=0.99))
diagnostics <- util$extract_hmc_diagnostics(fit)

print(util$check_all_hmc_diagnostics(diagnostics))

samples <- util$extract_expectand_vals(fit)
names <- c(grep('alpha_BAI', names(samples), value = TRUE),
           grep('beta_GST2', names(samples), value = TRUE),
           grep('sigma_BAI', names(samples), value = TRUE),
           grep('alpha_sc', names(samples), value = TRUE),
           grep('beta_GST1', names(samples), value = TRUE),
           grep('gamma_current', names(samples), value = TRUE),
           grep('gamma_lag', names(samples), value = TRUE),
           grep('sigma_sc', names(samples), value = TRUE))

base_samples <- util$filter_expectands(samples,names)
print(util$check_all_expectand_diagnostics(base_samples))
print(fit, pars = names)

pdf("figures/priorPosteriorPlotSimple.pdf", height = 9, width = 9)
par(mfrow = c(3,3))
util$plot_expectand_pushforward(samples[['alpha_BAI']], 50, display_name = "alpha_BAI")
curve(dlnorm(x, 5,2),
      add = TRUE,
      col = "blue",
      lwd = 2)
abline(v = 2, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['beta_GST2']], 50, display_name = "beta_GST2")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.05, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['sigma_BAI']], 50, display_name = "sigma_BAI")
curve(dnorm(x, 0, 5),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.5, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['alpha_sc']], 50, display_name = "alpha_sc")
curve(dnorm(x, 0.5, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.5, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['beta_GST1']], 50, display_name = "beta_GST1")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 0.3, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['gamma_current']], 50, display_name = "gamma_current")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = -0.2, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['gamma_lag']], 50, display_name = "gamma_lag")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = -0.1, col = "red", lwd = 2)

util$plot_expectand_pushforward(samples[['sigma_sc']], 50, display_name = "sigma_sc")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
abline(v = 1, col = "red", lwd = 2)

dev.off()
