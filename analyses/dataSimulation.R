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
beta_growth1 <- 0.04
beta_growth2 <- 3.2
gamma <- 0.01
phi <- 2
sigma <- 0.2

C <- alpha + alpha_site[site] + beta_dbh * DBH + beta_GST * GST

frac_R <- runif(N, 0.1, 0.5)
R <- C * frac_R
G <- C - R 

sc_pred <- R / gamma
sc <- rnegbin(N, mu = sc_pred, theta = phi)

rw_pred <- (((G / beta_growth1 + DBH^beta_growth2)^(1 / beta_growth2)) - DBH) / 2
rw <- rlnorm(N, meanlog = log(rw_pred), sdlog = sigma)

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
           grep('sigma', names(samples), value = TRUE))

base_samples <- util$filter_expectands(samples,names)
print(util$check_all_expectand_diagnostics(base_samples))
pdf("pairsPlot.pdf", height = 9, width = 9)
par(mfrow = c(3,3))
util$plot_div_pairs(names, names, samples, diagnostics)
dev.off()

pdf("figures/priorPosteriorPlot.pdf", height = 9, width = 9)
par(mfrow = c(3,3))
util$plot_expectand_pushforward(samples[['alpha']], 50, display_name = "alpha")
curve(dlnorm(x, log(90), 0.5),
      add = TRUE,
      col = "blue",
      lwd = 2)

util$plot_expectand_pushforward(samples[['alpha_site[1]']], 50, display_name = "alpha_site[1]")
curve(dnorm(x, 0, 10),
      add = TRUE,
      col = "blue",
      lwd = 2)

util$plot_expectand_pushforward(samples[['beta_dbh']], 50, display_name = "beta_dbh")
curve(dnorm(x, 0, 1),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)

util$plot_expectand_pushforward(samples[['beta_GST']], 50, display_name = "beta_GST")
curve(dnorm(x, 0, 5),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)

util$plot_expectand_pushforward(samples[['beta_growth1']], 50, display_name = "beta_growth1")
curve(dlnorm(x, 1, 0.5),
      add = TRUE,
      col = "blue",
      lwd = 2)

util$plot_expectand_pushforward(samples[['beta_growth2']], 50, display_name = "beta_growth2")
curve(dlnorm(x, 1, 0.5),
      add=TRUE, col="blue", lwd=2, lty=2)

util$plot_expectand_pushforward(samples[['gamma']], 50, display_name = "gamma")
curve(dnorm(x, 0, 5),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)

util$plot_expectand_pushforward(samples[['phi']], 50, display_name = "phi")
curve(dnorm(x, 0, 5),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)

util$plot_expectand_pushforward(samples[['sigma']], 50, display_name = "sigma")
curve(dnorm(x, 0, 5),
      add = TRUE,
      col = "blue",
      lwd = 2,
      lty = 2)
dev.off()

###for generative model###
N <- 20
DBH <- runif(N, 20, 40)
GST <- runif(N, 15, 25)

C <- numeric(N)
gamma <- runif(N, 0.1, 0.5)

G <- numeric(N)
G[1] <- 48
R <- numeric(N)
R[1] <- 32

rw <- numeric(N)
rw[1] <- 0.2
sc <- numeric(N)
sc[1] <- 1000

# true parameters
alpha <- 50
beta_dbh <- 0.5
beta_GST <- 0.5
sigma_C <- 5

theta <- 0.005
sigma_rw <- 0.1

eta <- 50


for(n in 1:N){

  muC <- alpha + beta_dbh*DBH[n] + beta_GST*GST[n]
  C[n] <- rnorm(1, muC, sigma_C)
  
  
  # carbon allocation
  G[n] <- (1-gamma[n]) * C[n]
  R[n] <- gamma[n] * C[n]
  
  # ring width
  rw[n] <- rnorm(1, theta*G[n], sigma_rw)
  
  # seed counts
  sc[n] <- rpois(1, eta*R[n])
  
  # update DBH
  DBH[n] <- DBH[n-1] + rw[n]
  
}

stanData <- list(
  N = N,
  sc = sc,
  rw = rw,
  DBH = DBH,
  GST = GST
)

mod <- stan_model(file='stan/tradeOffGenerative.stan')


fit <- stan(file='stan/tradeOffGenerative.stan', data=stanData, seed=112234, control=list(adapt_delta=0.99))
