#### Growth-reproduction trade-off ####
## Started by Mao ##
## Feb-28-2026 ##

library(MASS) 
library(rstan)

setwd("C:/PhD/Project/PhD_thesis/mast_growth/analyses")
set.seed(123)

N <- 50
NSite <- 3

DBH <- runif(N, 20, 40)
GST <- runif(N, 15, 25)
GSP <- runif(N, 50, 150)
site <- sample(1:NSite, N, replace = TRUE)

alpha <- 10
alpha_site <- c(0, 2, -1)
beta_dbh <- 0.5
beta_GST <- 1
beta_GSP <- 0.3
beta_growth1 <- 5
beta_growth2 <- 2
gamma <- 0.01
phi <- 2
sigma <- 0.1

C <- alpha + alpha_site[site] + beta_dbh * DBH + beta_GST * GST + beta_GSP * GSP

frac_R <- runif(N, 0.1, 0.5)
R <- C * frac_R
G <- C - R 

sc_pred <- R / gamma
sc <- rnegbin(N, mu = sc_pred, theta = phi)

rw_pred <- (G / (0.5 * beta_growth1))^(1 / beta_growth2)
rw <- rlnorm(N, meanlog = log(rw_pred), sdlog = sigma)

stanData <- list(
  N = N,
  NSite = NSite,
  site = site,
  sc = sc,
  rw = rw,
  DBH = DBH,
  GST = GST,
  GSP = GSP
)

tradeOff <-stan_model("stan/tradeOff.stan")
fit <- sampling(tradeOff, stanData, 
                iter = 4000, warmup = 3000,
                chains = 4)
