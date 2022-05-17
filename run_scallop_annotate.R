set.seed(42)
library(tidyverse)
library(tidybayes)
library(Cairo)
library(here)
library(magrittr)
library(rstan)
library(Matrix)
library(ggridges)
library(rstanarm)
library(geosphere)
library(ggridges)

stan_data <- readRDS(here("results", "scallop_stan_data_d0_grid_010522.rds"))

warmups <- 500
total_iterations <- 1000
max_treedepth <-  10
n_chains <-  1
n_cores <- 1
stan_model_fit <- stan(file = here::here("src","process_sdm_d0.stan"), # check that it's the right model!
                       data = stan_data,
                       chains = n_chains,
                       warmup = warmups,
                       init = list(list(log_mean_recruits = rep(log(1000), np),
                                        theta_d = 1)),
                       iter = total_iterations,
                       cores = n_cores,
                       refresh = 100,
                       control = list(max_treedepth = max_treedepth,
                                      adapt_delta = 0.85)
)

saveRDS(stan_model_fit, here("results", "scallop_stan_fit_grid_d0_010522.rds"))# stan_model_fit <- readRDS("results/scallop_stan_fit_122021.rds")
