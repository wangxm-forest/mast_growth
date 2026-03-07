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
source('stan_utility.R')
lsf.str()
set.seed(123)
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

