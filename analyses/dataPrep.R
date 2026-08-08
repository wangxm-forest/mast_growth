###################Started by Mao####################
###################July 23,2026######################

rm(list = ls())
options(stringsAsFactors = FALSE)
options(mc.cores = parallel::detectCores())

library("dplR")
library("dplyr")
library("ggplot2")
library("rstan")
setwd("C:/PhD/Project/PhD_thesis/mast_growth/analyses")

# Read in the ring width data
AB08THSE <- read.rwl("C:/PhD/Project/PhD_thesis/mast_growth/data/measurement/crossDated/AB08_TSHE_dated.rwl", format = "auto")

treeID <- sub("\\(.*\\)$", "", colnames(AB08THSE))

# Average across cores from the same tree
ring <- sapply(split(seq_along(treeID), treeID), function(i) {
  if (length(i) == 1) {
    AB08THSE[, i]
  } else {
    rowMeans(AB08THSE[, i, drop = FALSE], na.rm = TRUE)
  }
})

rownames(ring) <- rownames(AB08THSE)
# Merge all available seed data
merge <- FALSE
if(merge){
seed <- read.csv("C:/PhD/Project/PhD_thesis/mast_growth/data/MORA_cleanseeds_2009-2017.csv")
y2018 <- read.csv("C:/PhD/Project/Masting/data/rawdata/sortedseeds/clean&notes/SortedSeeds_MORA_2018.csv",header = TRUE)
y2019 <- read.csv("C:/PhD/Project/Masting/data/rawdata/sortedseeds/clean&notes/SortedSeeds_MORA_2019.csv",header = TRUE)
y2020 <- read.csv("C:/PhD/Project/Masting/data/rawdata/sortedseeds/clean&notes/SortedSeeds_MORA_2020.csv",header = TRUE)
y2021 <- read.csv("C:/PhD/Project/Masting/data/rawdata/sortedseeds/clean&notes/SortedSeeds_MORA_2021.csv",header = TRUE)
y2022 <- read.csv("C:/PhD/Project/Masting/data/rawdata/sortedseeds/clean&notes/SortedSeeds_MORA_2022.csv",header = TRUE)
y2023 <- read.csv("C:/PhD/Project/Masting/data/rawdata/sortedseeds/clean&notes/SortedSeeds_MORA_2023.csv",header = TRUE)
y2024 <- read.csv("C:/PhD/Project/Masting/data/rawdata/sortedseeds/clean&notes/SortedSeeds_MORA_2024.csv",header = TRUE)

y1824 <- rbind(y2018,y2019,y2020,y2021,y2022,y2023,y2024)
colnames(y1824)
colnames(seed)
y1824 <- y1824[-11]
colnames(y1824) <- c('year','stand','trap','species','filledseeds','emptyseeds','cones','conefilledseeds','coneemptyseeds','totconeseeds')
fullSeed <- rbind(y1824,seed)
write.csv(fullSeed,"C:/PhD/Project/PhD_thesis/mast_growth/data/seedMORAFull.csv")
}

# Read in seed data
seed <- read.csv("C:/PhD/Project/PhD_thesis/mast_growth/data/seedMORAFull.csv",header = TRUE)

# Make new columns for total seeds and total filled seeds
# Since cone filled seeds data are not available during 2009-2012, I will make a dataframe separately for total seeds and total filled seeds
seed$totalseeds <- seed$filledseeds + seed$emptyseeds + seed$totconeseeds

seed$totalfilled <- seed$filledseeds + seed$conefilledseeds
seed_stand_total <- aggregate(
  totalseeds ~ year + stand + species,
  data = seed,
  FUN = sum,
  na.rm = TRUE
)

seed_stand_filled <- aggregate(
  totalfilled ~ year + stand + species,
  data = seed,
  FUN = sum,
  na.rm = TRUE
)
# subset seed data

seed_sub_total <- seed_stand_total[seed_stand_total$stand == "AB08" & seed_stand_total$species == "TSHE", ]
seed_sub_filled <- seed_stand_filled[seed_stand_filled$stand == "AB08" & seed_stand_filled$species == "TSHE", ]


# Read in the DBH data
dbh <- read.csv("C:/PhD/Project/PhD_thesis/mast_growth/data/dbhMORA.csv",header = TRUE)

mora_stand <- c("AB08", "AV06", "TO04", "TA01", "AO03", "AG05", "AE10", "AM16")
dbh <- dbh[dbh$STANDID %in% mora_stand, ]

latest_dbh <- do.call(rbind, lapply(split(dbh, list(dbh$STANDID, dbh$TAG)), function(x) {
  x[which.max(x$YEAR), ]
}))

# Subset for only trees we cored

trees <- colnames(ring)

dbh_trees <- latest_dbh[latest_dbh$TAG %in% trees, ]


# Make ring width data in the correct format
ringwidth_reshape <- data.frame(
  TAG = rep(colnames(ring), each = nrow(ring)),
  year    = rep(as.numeric(rownames(ring)), times = ncol(ring)),
  ringWidth = as.vector(as.matrix(ring))
)
ringwidth_reshape <- ringwidth_reshape[!is.na(ringwidth_reshape$ringWidth), ]

# Prepare BAI data
# Calculate radius per year, I used genAI to help me finish writing this function

reconstruct_tree <- function(tree_tag, dbh_trees, ringwidth_reshape) {
  dbh_row  <- dbh_trees[dbh_trees$TAG == tree_tag, ]
  ref_year <- dbh_row$YEAR
  r_ref    <- dbh_row$DBH / 2

  tree_rings <- ringwidth_reshape[ringwidth_reshape$TAG == tree_tag, ]
  tree_rings <- tree_rings[order(tree_rings$year), ]

  # Check the reference year actually falls within the tree's ring-width span
  if (ref_year < min(tree_rings$year) || ref_year > max(tree_rings$year)) {
    warning(paste("Tree", TAG, "- DBH year", ref_year,
                   "falls outside ring-width data range; skipping"))
    return(NULL)
  }

  tree_rings$cum_width <- cumsum(tree_rings$ringWidth)   # running total up to & including each year

  cum_at_ref <- tree_rings$cum_width[tree_rings$year == ref_year]
  if (length(cum_at_ref) == 0) {
    warning(paste("Tree", TAG, "- no ring width recorded exactly at DBH year", ref_year))
    return(NULL)
  }

  tree_rings$radius <- r_ref + tree_rings$cum_width - cum_at_ref
  tree_rings$TAG <- tree_tag
  tree_rings
}

all_trees <- unique(ringwidth_reshape$TAG)
recon_list <- lapply(all_trees, reconstruct_tree, dbh_trees = dbh_trees, ringwidth_reshape = ringwidth_reshape)
recon <- do.call(rbind, recon_list)

recon <- recon[order(recon$TAG, recon$year), ]
recon$BA <- pi * recon$radius^2
recon$BAI <- ave(recon$BA, recon$TAG, FUN = function(x) c(NA, diff(x)))



# Prepare seed data for filled seeds
bai_2013_2024 <- recon[recon$year >= 2013 & recon$year <= 2024 & !is.na(recon$BAI), ]
all_years <- 2013:2024
N_years <- length(all_years)
year_lookup <- data.frame(year = all_years, year_idx = 1:N_years)
seed_sub_filled <- merge(seed_sub_filled, year_lookup, by = "year")
seed_sub_filled <- seed_sub_filled[seed_sub_filled$year_idx >= 2, ]
bai_2013_2024 <- merge(bai_2013_2024, year_lookup, by = "year")

stan_data_growth <- list(
  N = nrow(bai_2013_2024),
  BAI = bai_2013_2024$BAI,
  year = bai_2013_2024$year_idx,
  N_years = N_years
)

stan_data_filled <- c(stan_data_growth, list(
  N_sc = nrow(seed_sub_filled),
  sc = seed_sub_filled$totalfilled,
  year_sc = seed_sub_filled$year_idx
))


# run stan model with filled seeds

mod <- stan_model(file='C:/PhD/Project/PhD_thesis/mast_growth/analyses/stan/simpleTradeOff.stan')

fit_filled <- stan(file='C:/PhD/Project/PhD_thesis/mast_growth/analyses/stan/simpleTradeOff.stan', data=stan_data_filled, seed=112234, control=list(adapt_delta=0.99))

# Plotting the parameters for model results with filled seeds
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

diagnostics <- util$extract_hmc_diagnostics(fit_filled)

print(util$check_all_hmc_diagnostics(diagnostics))

samples <- util$extract_expectand_vals(fit_filled)
names <- c(grep('alpha_BAI', names(samples), value = TRUE),
           grep('sigma_BAI', names(samples), value = TRUE),
           grep('alpha_sc', names(samples), value = TRUE),
           grep('gamma_current', names(samples), value = TRUE),
           grep('gamma_lag', names(samples), value = TRUE),
           grep('sigma_sc', names(samples), value = TRUE))

base_samples <- util$filter_expectands(samples,names)
print(util$check_all_expectand_diagnostics(base_samples))
print(fit_filled, pars = names)

post <- as.data.frame(fit_filled)

params_df <- data.frame(
  parameter = c("gamma_current", "gamma_lag"),
  mean = c(mean(post$gamma_current), mean(post$gamma_lag)),
  lower = c(quantile(post$gamma_current, 0.025), quantile(post$gamma_lag, 0.025)),
  upper = c(quantile(post$gamma_current, 0.975), quantile(post$gamma_lag, 0.975))
)

p <- ggplot(params_df, aes(x = parameter, y = mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 1, color = "#1D5145")  +
  scale_x_discrete(labels = c(
    "gamma_current" = expression(gamma[current]),
    "gamma_lag"     = expression(gamma[lag])
  )) +
  labs(y = "Effect on log(expected seed count)", x = NULL) +
  theme_minimal(base_size = 25) +
  coord_flip() + theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

p
ggsave("figures/parametersFilledSeeds.png", p, bg="transparent")

# Manually plot the model results with CI, adapted this code from egret budseed model visualization
Gbar_draws  <- as.matrix(fit_filled, pars = "Gbar")
Gbar_mean <- colMeans(Gbar_draws)

Gbar_seq <- seq(min(Gbar_mean), max(Gbar_mean), length.out = 50)

G_obs  <- Gbar_mean[stan_data_filled$year_sc]
sc_obs_log <- log(stan_data_filled$sc) 

pred_draws_log <- sapply(Gbar_seq, function(g) {
  post$alpha_sc + post$gamma_current * g + post$gamma_lag * mean(Gbar_mean)
})

# Make the dataframe so I can use ggplot2 to make a plot with transparent background
pred_df_log <- data.frame(
  G = Gbar_seq,
  mean = colMeans(pred_draws_log),
  low  = apply(pred_draws_log, 2, quantile, probs = 0.1),
  high = apply(pred_draws_log, 2, quantile, probs = 0.9)
)

obs_df_log <- data.frame(G = G_obs, sc = sc_obs_log)

p <- ggplot() +
  geom_ribbon(data = pred_df_log, aes(x = G, ymin = low, ymax = high),
              fill = "#9BAEAA", alpha = 0.3) +
  geom_line(data = pred_df_log, aes(x = G, y = mean),
            color = "#1D5145", size = 1.2) +
  geom_point(data = obs_df_log, aes(x = G, y = sc),
             size = 3, color = "black") +
  labs(x = "Growth deviation (G)", y = "log (seed count)")  +
  theme_minimal(base_size = 25) + theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

p

ggsave("figures/predictsFilledGvsR.png", p, bg="transparent")

# Prepare seed data for all seeds
bai_2010_2024 <- recon[recon$year >= 2010 & recon$year <= 2024 & !is.na(recon$BAI), ]
all_years <- 2010:2024
N_years <- length(all_years)
year_lookup <- data.frame(year = all_years, year_idx = 1:N_years)
seed_sub_total <- merge(seed_sub_total, year_lookup, by = "year")
seed_sub_total <- seed_sub_total[seed_sub_total$year_idx >= 2, ]
bai_2010_2024 <- merge(bai_2010_2024, year_lookup, by = "year")

stan_data_growth <- list(
  N = nrow(bai_2010_2024),
  BAI = bai_2010_2024$BAI,
  year = bai_2010_2024$year_idx,
  N_years = N_years
)

stan_data_total <- c(stan_data_growth, list(
  N_sc = nrow(seed_sub_total),
  sc = seed_sub_total$totalseeds,
  year_sc = seed_sub_total$year_idx
))

fit_all <- stan(file='C:/PhD/Project/PhD_thesis/mast_growth/analyses/stan/simpleTradeOff.stan', data=stan_data_total, seed=112234, control=list(adapt_delta=0.99))


# Plotting the parameters for model results with all seeds
util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

diagnostics <- util$extract_hmc_diagnostics(fit_all)

print(util$check_all_hmc_diagnostics(diagnostics))

samples <- util$extract_expectand_vals(fit_all)
names <- c(grep('alpha_BAI', names(samples), value = TRUE),
           grep('sigma_BAI', names(samples), value = TRUE),
           grep('alpha_sc', names(samples), value = TRUE),
           grep('gamma_current', names(samples), value = TRUE),
           grep('gamma_lag', names(samples), value = TRUE),
           grep('sigma_sc', names(samples), value = TRUE))

base_samples <- util$filter_expectands(samples,names)
print(util$check_all_expectand_diagnostics(base_samples))
print(fit_all, pars = names)

post <- as.data.frame(fit_all)

params_df <- data.frame(
  parameter = c("gamma_current", "gamma_lag"),
  mean = c(mean(post$gamma_current), mean(post$gamma_lag)),
  lower = c(quantile(post$gamma_current, 0.025), quantile(post$gamma_lag, 0.025)),
  upper = c(quantile(post$gamma_current, 0.975), quantile(post$gamma_lag, 0.975))
)

p <- ggplot(params_df, aes(x = parameter, y = mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 1, color = "#1D5145")  +
  scale_x_discrete(labels = c(
    "gamma_current" = expression(gamma[current]),
    "gamma_lag"     = expression(gamma[lag])
  )) +
  labs(y = "Effect on log(expected seed count)", x = NULL) +
  theme_minimal(base_size = 25) +
  coord_flip() + theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

p
ggsave("figures/parametersAllSeeds.png", p, bg="transparent")

# Manually plot the model results with CI, adapted this code from egret budseed model visualization
Gbar_draws  <- as.matrix(fit_all, pars = "Gbar")
Gbar_mean <- colMeans(Gbar_draws)

Gbar_seq <- seq(min(Gbar_mean), max(Gbar_mean), length.out = 50)

G_obs  <- Gbar_mean[stan_data_total$year_sc]
sc_obs_log <- log(stan_data_total$sc) 

pred_draws_log <- sapply(Gbar_seq, function(g) {
  post$alpha_sc + post$gamma_current * g + post$gamma_lag * mean(Gbar_mean)
})

# Make the dataframe so I can use ggplot2 to make a plot with transparent background
pred_df_log <- data.frame(
  G = Gbar_seq,
  mean = colMeans(pred_draws_log),
  low  = apply(pred_draws_log, 2, quantile, probs = 0.1),
  high = apply(pred_draws_log, 2, quantile, probs = 0.9)
)

obs_df_log <- data.frame(G = G_obs, sc = sc_obs_log)

p <- ggplot() +
  geom_ribbon(data = pred_df_log, aes(x = G, ymin = low, ymax = high),
              fill = "#9BAEAA", alpha = 0.3) +
  geom_line(data = pred_df_log, aes(x = G, y = mean),
            color = "#1D5145", size = 1.2) +
  geom_point(data = obs_df_log, aes(x = G, y = sc),
             size = 3, color = "black") +
  labs(x = "Growth deviation (G)", y = "log (seed count)")  +
  theme_minimal(base_size = 25) + theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

p

ggsave("figures/predictsTotalGvsR.png", p, bg="transparent")
