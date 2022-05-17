#############
# load packages and data
#############

set.seed(42)
# library(tidyverse)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(tidybayes) # when using rstudio server (r3.5.2): devtools::install_version(package = "tidybayes", version = "1.0.3")
library(Cairo)
library(here)
library(magrittr)
library(rstan)
library(Matrix)
library(ggridges)
# library(rstanarm)
library(geosphere)
library(ggridges)

funs <- list.files("functions")
sapply(funs, function(x) source(file.path("functions",x)))

rstan_options(javascript=FALSE, auto_write =TRUE)

dat <- read_csv(here("processed-data","scallop_catch_at_length.csv")) %>%
  filter(year > 1967,
         year < 2015) %>%
  # add grid IDs for half-degree cells
  mutate(grid_lat = case_when(lat-trunc(lat)<0.5 ~ trunc(lat)+0.25,
                              lat-trunc(lat)>=0.5 ~ 
                                trunc(lat)+0.75),
         grid_lon = case_when(abs(lon-trunc(lon))<0.5 ~ 
                                trunc(lon)-0.25,
                              abs(lon-trunc(lon))>=0.5 ~ 
                                trunc(lon)-0.75),
         grid = paste0(grid_lat, grid_lon))

# identify cells that have scallop and sbt data for all years
# before filtering = 232 cells
# with scallop data all years = 71 cells
# AND with sea bottom temp all years = 40 cells
has_data_all_years <- dat %>%
  group_by(grid, year) %>%
  summarise(mean_btemp = mean(btemp, na.rm = T),
            .groups = "drop") %>%
  filter(is.na(mean_btemp) == F) %>%
  group_by(grid) %>%
  tally() %>%
  filter(n >= 47)

# how much variation is there in length frequency over time?
# see flounder code

##########
# prep data for fitting
##########

# identify top_n most abundant patches
top_n <- 100
top_patches <- dat %>%
  group_by(grid) %>% 
  summarise(total = sum(number_at_length)) %>% 
  arrange(-total) %>% 
  slice(1:top_n) %>% # CHECK THIS MANUALLY TO BE SURE IT'S SANE, AND PATCHES ARE CONTIGUOUS
  left_join(distinct(select(dat, grid_lat, grid_lon, grid)), by = "grid") %>%
  # remove patches with missing years (to allow dens matrix creation)
  filter(grid %in% has_data_all_years$grid)

# check to see which patches were selected
ggplot(top_patches) +
  geom_tile(aes(x = grid_lon, y = grid_lat, fill = total))

dat_train_lengths <- dat %>% 
  group_by(length, year, grid) %>% 
  summarise(sum_num_at_length = sum(number_at_length)) %>% 
  filter(grid %in% top_patches$grid)%>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(grid)))

dat_train_dens <- dat %>% 
  filter(grid %in% top_patches$grid) %>% 
  group_by(haulid) %>% 
  mutate(dens = sum(number_at_length)) %>% # get total no. fish in each haul, of any size
  group_by(year, grid) %>% 
  summarise(mean_dens = mean(dens)) %>%  # get mean density (all sizes) / haul for the patch*year combo 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(grid)))

dat_train_dens70 <- dat %>% 
  filter(grid %in% top_patches$grid,
         length >= 7) %>% 
  group_by(haulid) %>% 
  mutate(dens = sum(number_at_length)) %>% # get total no. fish in each haul, of any size
  group_by(year, grid) %>% 
  summarise(mean_dens = mean(dens)) %>%  # get mean density (all sizes) / haul for the patch*year combo 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(grid)))

# get patch area 
patchdat <- dat %>% 
  select(grid, grid_lat, grid_lon) %>%
  distinct() %>%
  mutate(max_lon = grid_lon+0.25,
            min_lon = grid_lon-0.25) %>% 
  rowwise() %>% 
  mutate(lon_dist = distGeo(p1=c(max_lon, grid_lat), p2=c(min_lon, grid_lat))/1000, # get distance between the furthest longitudes in km, at the midpoint of the lat band 
         patch_area_km2 = lon_dist * 111) %>%  # 1 degree latitude = 111 km 
  select(grid, patch_area_km2) %>% 
  filter(grid %in% top_patches$grid) %>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(grid)))

dat_train_sbt <- dat %>%   
  group_by(grid, year) %>% 
  summarise(sbt = mean(btemp, na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(grid %in% top_patches$grid)%>% 
  mutate(patch = as.integer(as.factor(grid)))

# set fixed parameters from stock assessment
loo.gb <- mean(c(143.9, 140.0, 148.6, 121.1, 144.9, 152.5, 143.6, 143.3, 145.1))
loo.ma <- mean(c(133.3, 151.8))
loo = mean(loo.gb, loo.ma) # mean from Hart & Chute: 142.6
k = 0.37 # means from Hart & Chute; GB: 0.33; MAB: 0.40
m = 0.35 # 2018 assessment: M = 0.35 for all years and ages.
f = 0.31 # average of 2008-2019 from 2018 Stock Assess & 2020 update
z = exp(-m-f)
age_at_maturity = 4
t0=-.2
cv= 0.2 # guess
min_age = 1
max_age = 14 # based on age data in scallop dredge survey


length_at_age_key <-
  generate_length_at_age_key(
    min_age = min_age,
    max_age = max_age,
    cv = cv,
    linf = loo,
    k = k,
    t0 = t0,
    time_step = 1,
    linf_buffer = 1.5
  )

length_at_age_key %>% 
  filter(age > 0) %>% 
  ggplot(aes(age, length_bin, fill = p_bin)) + 
  geom_hline(aes(yintercept = loo)) +
  geom_tile() + 
  scale_fill_viridis_c()

l_at_a_mat <- length_at_age_key %>% 
  select(age, length_bin, p_bin) %>% 
  pivot_wider(names_from = length_bin, values_from = p_bin) %>% 
  ungroup() %>% 
  select(-age) %>% 
  as.matrix()

# get time dimension
years <- sort(unique(dat_train_lengths$year)) 
ny <- length(years)
ny_proj <- 10

#get other dimensions
patches <- sort(unique(dat_train_lengths$grid))
np = length(patches) 

lbins <- unique(length_at_age_key$length_bin)
# lbins <- sort(unique(dat_train_lengths$length))
n_lbins <- length(lbins) 

n_ages <- nrow(l_at_a_mat)

# now that years are defined above, convert them into indices in the datasets
dat_train_dens$year = as.integer(as.factor(dat_train_dens$year))
dat_train_dens70$year = as.integer(as.factor(dat_train_dens70$year))
dat_train_lengths$year = as.integer(as.factor(dat_train_lengths$year))
dat_train_sbt$year= as.integer(as.factor(dat_train_sbt$year))

# make matrices/arrays from dfs
len <- array(0, dim = c(np, n_lbins, ny)) 
for(p in 1:np){
  print(p)
  for(l in 1:n_lbins){
    for(y in 1:ny){
      tmp <- dat_train_lengths %>% filter(patch==p, round(length)==lbins[l], year==y) 
      if (nrow(tmp) > 0){
      len[p,l,y] <- tmp$sum_num_at_length
      }
    }
  }
}

# save for later
# saveRDS(len, "results/len_100_top_patches_010522.rds")
# read in to save time
# len <- readRDS("results/len_100_top_patches_010522.rds")

plot(len[4,,20])

dens <- array(NA, dim=c(np, ny))
for(p in 1:np){
  print(p)
  for(y in 1:ny){
    print(paste0("doing patch ",p, ", year ", y))
    tmp2 <- dat_train_dens %>% filter(patch==p, year==y) %>% 
      left_join(patchdat, by = c("grid","patch"))%>% 
      mutate(mean_dens = mean_dens * patch_area_km2)
    dens[p,y] <- tmp2$mean_dens
  }
}

sbt <- array(NA, dim=c(np,ny))
for(p in 1:np){
  for(y in 1:ny){
    tmp3 <- dat_train_sbt %>% filter(patch==p, year==y) 
    sbt[p,y] <- tmp3$sbt
  }
}

a <- seq(min_age, max_age)

check <- a %*% l_at_a_mat

######
# fit model
######
saveRDS(stan_data, here("results", "scallop_stan_data_top100_minus_64_w_NA.rds"))
# stan_data <- readRDS(here("results", "scallop_stan_data_top100_minus_64_w_NA.rds")) # top 100 minus 64 with missing data

stan_data <- list(
  np=np,
  n_ages=n_ages,
  ny_train=ny,
  n_lbins=n_lbins,
  n_p_l_y = len,
  abund_p_y = dens,
  sbt = sbt,
  m=m,
  k=k,
  loo=loo,
  t0=t0,
  cv=cv,
  length_50_sel_guess=70, # THIS IS A RANDOM GUESS, I can't find help in the stock assessment
  n_lbins = n_lbins, 
  age_sel = 0,
  bin_mids=lbins+0.5, # also not sure if this is the right way to calculate the midpoints
  sel_100 = 3, # not sure if this should be 2 or 3. it's age 2, but it's the third age category because we start at 0, which I think Stan will classify as 3...?
  age_at_maturity = age_at_maturity,
  patcharea = patchdat$patch_area_km2,
  l_at_a_key = l_at_a_mat,
  do_dirichlet = 1
  
)

warmups <- 1000
total_iterations <- 2000
max_treedepth <-  10
n_chains <-  1
n_cores <- 1
thin <- 10 
np <- stan_data$np
stan_model_fit <- stan(file = here::here("src","process_sdm_d0_grid.stan"), # check that it's the right model!
                      data = stan_data,
                      chains = n_chains,
                      warmup = warmups,
                      init = list(list(log_mean_recruits = rep(log(1000), np),
                                       theta_d = 1)),
                      iter = total_iterations,
                      thin = thin,
                      cores = n_cores,
                      save_warmup = FALSE,
                      refresh = 100,
                      control = list(max_treedepth = max_treedepth,
                                     adapt_delta = 0.85)
)

saveRDS(stan_model_fit, here("results", "scallop_stan_fit_grid_d0_011122.rds"))
# stan_model_fit <- readRDS(here("results", "scallop_stan_fit_grid_d0_011122.rds"))

a = rstan::extract(stan_model_fit, "theta_d")
# write_rds(stan_model_fit,"sigh.rds")
# hist(a$sigma_obs)
# rstanarm::launch_shinystan(stan_model_fit)

# assess abundance fits

abund_p_y <- dat_train_dens %>%
  left_join(patchdat, by = c("grid", "patch")) %>% 
  group_by(patch, year) %>% 
  summarise(abundance = sum(mean_dens *patch_area_km2)) %>% 
  ungroup()

# just scallops >= 70 mm
abund_p_y70 <- dat_train_dens70 %>%
  left_join(patchdat, by = c("grid", "patch")) %>% 
  group_by(patch, year) %>% 
  summarise(abundance = sum(mean_dens *patch_area_km2)) %>% 
  ungroup()
  

abund_p_y_hat <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat[patch,year])
abund_p_y_hat70 <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat70[patch,year])


abund_p_y70_plot = aggregate(abund_p_y_hat70$dens_p_y_hat70, by = abund_p_y_hat70[,4:5], mean) %>%
  rename(dens_p_y_hat70 = x) 

abund_p_y_plot = aggregate(abund_p_y_hat$dens_p_y_hat, by = abund_p_y_hat[,4:5], median) %>%
  rename(dens_p_y_hat = x) %>%
  left_join(abund_p_y) %>%
  rename(dens_p_y = abundance) %>%
  left_join(abund_p_y70_plot) %>%
  left_join(abund_p_y70) %>%
  rename(dens_p_y70 = abundance)

# observed vs. predicted scatterplot
abund_p_y_plot %>%
  ggplot() + 
  geom_point(aes(dens_p_y+1, dens_p_y_hat+1), 
             alpha = 0.5, color = "steelblue") +
  scale_x_log10() + 
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw()
# ggsave("results/osb_vs_predicted_run011122.png", height = 4, width = 6, dpi = 400)
  
patch_filter = 31:36
abund_p_y_plot %>%
  filter(patch %in% patch_filter) %>%
  ggplot() + 
  stat_lineribbon(data = filter(abund_p_y_hat, patch %in% patch_filter), aes(year, dens_p_y_hat)) + 
  geom_point(aes(year, dens_p_y), color = "steelblue", size = 1) +
  facet_wrap(~patch, scales = "free_y")
ggsave("results/abund_p_y_hat_d0_grid31_36_run011122.png", dpi = 400, height = 10, width = 10)

patch_filter = 13:18
abund_p_y_hat %>%
  filter(patch %in% patch_filter) %>%
  ggplot(aes(year, dens_p_y_hat+1)) + 
  stat_lineribbon() + 
  geom_point(data = filter(abund_p_y, patch %in% patch_filter), aes(year, abundance+1), 
             color = "steelblue", size = 1) +
  facet_wrap(~patch, scales = "free_y") +
  scale_y_log10()
ggsave("results/abund_p_y_hat_d0_grid13_18_log_run011122.png", dpi = 400, height = 10, width = 10)

# assess length comp fits

n_p_l_y_hat <- tidybayes::gather_draws(stan_model_fit, n_p_l_y_hat[year,patch,length]) # originally
# saveRDS(n_p_l_y_hat, "results/n_p_l_y_hat_122021.rds")
# saveRDS(n_p_l_y_hat, "results/n_p_l_y_hat_010422.rds")
saveRDS(n_p_l_y_hat, "results/n_p_l_y_hat_011122.rds")
n_p_l_y_hat <- readRDS("results/n_p_l_y_hat_011122.rds")
# neff <- tidybayes::gather_draws(stan_model_fit, n_eff[patch,year], n = 500)

#neff <- tidybayes::gather_draws(stan_model_fit, n_eff[patch,year], n = 500)

p = 17

dat_train_lengths <- dat_train_lengths %>% 
  group_by(patch, year) %>% 
  mutate(p_length = sum_num_at_length / sum(sum_num_at_length))

dat_train_lengths %>% 
  group_by(year,patch) %>% 
  summarise(n = sum(sum_num_at_length)) %>% 
  ggplot(aes(year, n, color =factor(patch))) + 
  geom_point() +
  theme(legend.position = "none")
# ggsave("results/n_by_patch_year_run011122.png", dpi = 400, height = 4, width = 6)

n_p_l_y_hat %>% 
  ungroup() %>% 
  filter(patch == p) %>% 
  group_by(patch, year, .iteration) %>% 
  mutate(pvalue = .value / sum(.value)) %>% 
  ggplot(aes(length, pvalue)) + 
  stat_lineribbon() + 
  geom_point(data = dat_train_lengths %>% filter(patch == p), aes(length,p_length), color = "red", alpha = 0.2) +
  facet_wrap(~year, scales = "free_y")
 # ggsave("results/n_p_l_y_hat_d0.png", dpi = 400, height = 10, width = 10)

# plot important parameters 
plot(stan_model_fit, pars=c('sigma_r','sigma_obs','d','width','Topt','alpha','theta','theta_d','log_f'))
#Topt going to about 20, width to about 8; zoom in on the others
plot(stan_model_fit, pars=c('sigma_r','sigma_obs','d','alpha','theta','theta_d','log_f'))


# save model run if desired
saveRDS(stan_model_fit, here("results","scallop_stan_fit_d0_01042022.rds"))

###########
# calculate summary statistics and evaluate range shifts
###########

# centroid position by year 
dat_centroid <- abund_p_y %>% 
  group_by(year) %>% 
  summarise(centroid_lat = weighted.mean(x=patch, w=abundance))

# model fit centroid -- should eventually estimate in model for proper SE -- just exploring here
est_centroid <- abund_p_y_hat %>% 
  group_by(year, .draw) %>%  # IS THIS SUPPOSED TO BE .ITERATION? CHECK WHEN MODEL IS RUN FOR LONGER 
  summarise(centroid_lat = weighted.mean(x=patch, w=dens_p_y_hat)) %>% 
  ungroup()

est_centroid %>% 
  ggplot(aes(year, centroid_lat)) + 
  stat_lineribbon() + 
  scale_fill_brewer() +
  geom_point(data = dat_centroid, aes(year, centroid_lat), color = "red") 
ggsave("results/centroid_shift_run011122.png", dpi = 400, width = 6, height = 4)
# centroid didn't shift at all!

# patch abundance fraction every year 
est_patch_abund <- abund_p_y_hat %>% 
  group_by(year, patch) %>% 
  summarise(abundance = mean(dens_p_y_hat))

abund_p_y %>% 
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 47, 4)) +
  scale_y_continuous(breaks=seq(1, 36, 1)) +
  labs(title="Observed")
# ggsave("results/abund_p_y_run011122.png", dpi = 400, width = 6, height = 4)

est_patch_abund %>% 
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 47, 4)) +
  scale_y_continuous(breaks=seq(1, 36, 1)) +
  labs(title="Estimated")
# ggsave("results/abund_p_y_estimated_run011122.png", dpi = 400, width = 6, height = 4)

# who's doing the colonizing?
# dat_train_lengths %>% 
#   group_by(patch, length) %>% 
#   arrange(year) %>% 
#   mutate(logratio = log(sum_num_at_length / lag(sum_num_at_length))) %>% 
#   filter(logratio < Inf, logratio > -Inf) %>% 
#   ungroup() %>% 
#   ggplot(aes(x=year, y=logratio, group=length, color=length)) +
#   geom_point() + 
#   geom_line() +
#   scale_color_viridis_c() +
#   facet_wrap(~patch)
# # ggsave("results/who_is_colonizing_d0.png", dpi = 400, width = 9, height = 7)

patch_filter <- 25:36
dat_train_lengths %>% 
  filter(patch %in% patch_filter) %>%
  ggplot(aes(x=year, y=sum_num_at_length, fill=length)) + 
  geom_bar(stat="identity") +
  facet_wrap(~patch) +
  scale_fill_viridis_c() +
  theme_bw()
# ggsave("results/sum_num_at_length_grid25_36_run011122.png", dpi = 400, width = 9, height = 7)
