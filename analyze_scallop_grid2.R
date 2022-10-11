# this file is based on analyze_summer_flounder.R last changed ~1/14/22

#############
# load packages and data
#############

set.seed(42)
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
# library(ggridges)
# library(rstanarm)
library(geosphere)
# library(ggridges)
# library(Amelia) # seems to interfere with Stan
funs <- list.files("functions")
sapply(funs, function(x) source(file.path("functions",x)))

rstan_options(javascript=FALSE, auto_write =TRUE)

# set the 
plotsave <- "results/run20221010b"

# read in climate data
clim <- read.csv("processed-data/climate_formatted.csv")
clim_avg <- clim # rename so clim can be used again if 1x1 grid spacing is used

#############
# make model decisions
#############
train_years <- 1980:2004 # set training year range
test_years <- 2005:2014 # set testing year range
season <- "fall" # "spring", "fall", or "both
grid_1x1 <- 1 # 1 = 1x1 degree grid, 0 = 0.5 x 0.5 degree grid
max_NA_dens <- 4 # set max number of NA dens values for each cell
max_NA_sbt <- 4 # set max number of NA sbt values for each cell
lagged_sbt <- 0 # 1 = sbt lagged 1 year; 0 = not
impute_mean_dens <- 0 # impute missing dens values (1) or no (0)
topn <- NULL # NULL if not using "keep only top n cells" option
use_custom_grid = 1
custom_grid_vector = c("39.5-72.5", "39.5-73.5") # NULL # vector of grid names you want to include; NULL = none
manual_selectivity = 1
do_dirichlet = 1
eval_l_comps = 0 # evaluate length composition data? 0=no, 1=yes
T_dep_mortality = 0 # 
T_dep_recruitment = 1 #
spawner_recruit_relationship = 0
run_forecast=1
time_varying_f = TRUE
btemp_meas <- "mean" # "min", "mean", or "max"
wt_at_age <- rep(1, 14) # not used in scallop model so far

if(time_varying_f==TRUE){
# the f-at-age data starts in 1982; fill in the previous years with the earliest year of data
dat_f_age_prep <- read_csv(here("processed-data","scallop_F_by_age.csv")) %>%
  rename(year = Year, age = Age, f = F)
# f_early <- expand_grid(year=seq(1972, 1981, 1), age=unique(dat_f_age_prep$age)) %>% 
#   left_join(dat_f_age_prep %>% filter(year==1982) %>% select(age, f)) 
# dat_f_age_prep <- bind_rows(dat_f_age_prep, f_early)
}

##########
# prep data for fitting
##########

# read in data
dat <- read_csv(here("processed-data","scallop_catch_at_length_mo.csv")) %>%
  filter(year %in% c(unique(c(train_years, test_years)),min(train_years)-1)) %>% # -1 ensures an earlier year to get lagged sbt temps  
  # add grid IDs for half-degree cells
  mutate(grid_lat = case_when(lat-trunc(lat)<0.5 ~ trunc(lat)+0.25,
                              lat-trunc(lat)>=0.5 ~ 
                                trunc(lat)+0.75),
         grid_lon = case_when(abs(lon-trunc(lon))<0.5 ~ 
                                trunc(lon)-0.25,
                              abs(lon-trunc(lon))>=0.5 ~ 
                                trunc(lon)-0.75),
         grid = paste0(grid_lat, grid_lon))

# rename grid cells if cell size is set to 1x1 degree
if(grid_1x1 == 1){
  dat <- dat %>%
    mutate(grid_lat = case_when(substr(grid_lat,4,5)=="75" ~
                                  grid_lat-0.25,
                                substr(grid_lat,4,5)=="25" ~
                                  grid_lat+0.25),
           grid_lon = case_when(substr(grid_lon,5,6)=="75" ~
                                  grid_lon+0.25,
                                substr(grid_lon,5,6)=="25" ~
                                  grid_lon-0.25),
           grid = paste0(grid_lat,grid_lon))
  
  clim_avg <- clim %>%
    mutate(grid_lat = case_when(substr(grid_lat,4,5)=="75" ~
                                  grid_lat-0.25,
                                substr(grid_lat,4,5)=="25" ~
                                  grid_lat+0.25),
           grid_lon = case_when(substr(grid_lon,5,6)=="75" ~
                                  grid_lon+0.25,
                                substr(grid_lon,5,6)=="25" ~
                                  grid_lon-0.25),
           grid = paste0(grid_lat,grid_lon)) %>%
    group_by(grid, year) %>%
    summarise(across(where(is.numeric), mean))
  
}

# subset data by season
if(season == "spring"){dat <- dat %>% filter(month < 7)}
if(season == "fall"){dat <- dat %>% filter(month > 7)}
if(season == "both"){dat <- dat}

# make a data frame with all grid cells & years to ensure a row for each
dat_all_grid_years <- expand.grid(grid = unique(dat$grid), 
                                  year = train_years) %>%
  mutate(haulid = "temp",
         length = NA,
         spp = "scallop",
         number_at_length = NA,
         btemp = NA,
         lat = NA,
         lon = NA,
         grid_lat = NA,
         grid_lon = NA
  )

# add the "all grid-years" df to the end of dat
dat_augmented <- dat %>%
  bind_rows(dat_all_grid_years)

# choose which cells to include and deal with missing data

use_patches = dat_augmented %>%
  group_by(haulid, grid, year, btemp) %>% 
  summarise(dens = sum(number_at_length)) %>% # , na.rm = T not used on purpose to count NAs
  group_by(year, grid) %>% 
  summarise(mean_dens = mean(dens, na.rm = T),
            mean_sbt = mean(btemp, na.rm = T),
            .groups = "drop") %>%  # get mean density (all sizes) / haul for the patch*year combo 
  arrange(grid, year) %>%
  group_by(grid) %>%
  summarise(dens_NA = sum(is.na(mean_dens)),
            sbt_NA = sum(is.na(mean_sbt)),
            mean_dens = mean(mean_dens, na.rm = T)) %>%
  filter(dens_NA <= max_NA_dens,
         sbt_NA <= max_NA_sbt)

if(is.null(topn) == FALSE){
use_patches <- use_patches  %>%
    arrange(-mean_dens) %>%
    slice(1:topn)
}

if(use_custom_grid == 1){
use_patches <- use_patches %>%
  filter(grid %in% custom_grid_vector) # I've been making these in the scallop_custom_grids.R script
}

# dat filtered to have no more than threshold NAs
dat <- dat_augmented %>%
  filter(dat_augmented$grid %in% use_patches$grid) %>% 
  mutate(patch = as.integer(as.factor(grid)))

# prep dat length data
dat_lengths <- dat %>% 
  group_by(length, year, grid, patch) %>% 
  summarise(sum_num_at_length = sum(number_at_length, na.rm = T),
            .groups = "drop") %>%
  filter(is.na(length) == FALSE)

# prep dat length data excluding small scallops (< 80 mm)
dat_lengths80 <- dat %>% 
  filter(length >= 8) %>%
  group_by(length, year, grid, patch) %>% 
  summarise(sum_num_at_length = sum(number_at_length, na.rm = T),
            .groups = "drop") %>%
  filter(is.na(length) == FALSE)

# subset appropriate climate variable (e.g., SBT min, mean, or max)
if(btemp_meas == "mean"){
  clim_avg <- clim_avg %>%
    select(grid, year, mean_temp) %>%
    rename(climvar = mean_temp)
}

if(btemp_meas == "min"){
  clim_avg <- clim_avg %>%
    select(grid, year, min_temp) %>%
    rename(climvar = min_temp)
}

if(btemp_meas == "max"){
  clim_avg <- clim_avg %>%
    select(grid, year, max_temp) %>%
    rename(climvar = max_temp)
}

# prep dat dens data
dat_dens <- dat %>% 
  group_by(haulid, grid, patch, year, btemp) %>% 
  summarise(dens = sum(number_at_length, na.rm = F),
            .groups = "drop") %>% # get total no. scallops in each haul, of any size
  group_by(year, grid, patch) %>% 
  summarise(mean_dens = mean(dens, na.rm = T),
            btemp_morely = mean(btemp, na.rm = T),
            n_haul = sum(haulid!="temp"),
            .groups = "drop") %>% # get mean density (all sizes) / haul for the patch*year combo 
  left_join(clim_avg, by = c("grid", "year"))
  

# prep dat dens data, excluding smaller scallops (< 80 mm)
dat_dens80 <- dat %>% 
  filter(length >= 8) %>%
  group_by(haulid, grid, patch, year, btemp) %>% 
  summarise(dens = sum(number_at_length, na.rm = F),
            .groups = "drop") %>% # get total no. scallops in each haul, of any size
  group_by(year, grid, patch) %>% 
  summarise(mean_dens = mean(dens, na.rm = T),
            btemp_morely = mean(btemp, na.rm = T),
            n_haul = sum(haulid!="temp"),
            .groups = "drop") %>% # get mean density (all sizes) / haul for the patch*year combo 
  left_join(clim_avg, by = c("grid", "year"))

# create patch ID and count patches
  patches <- sort(unique(use_patches$grid))
  np = length(patches) 
  
  # add in lagged sea bottom temperature
  dat_dens <- dat_dens %>%
    arrange(grid, year) %>%
    mutate(lclimvar = lag(climvar)) %>%
    filter(year != (min(train_years)-1)) # this also removes all the incorrect values in year min(year)-1 that carried over from other grid cells in the lag
  
  # add in lagged sea bottom temperature (for version without small scallops < 80 mm)
  dat_dens80 <- dat_dens80 %>%
    arrange(grid, year) %>%
    mutate(lclimvar = lag(climvar)) %>%
    filter(year != (min(train_years)-1)) # this also removes all the incorrect values in year min(year)-1 that carried over from other grid cells in the lag
  
  # plot to see which patches were selected
dat_dens %>%
    group_by(grid, patch) %>%
    summarise(mean_dens = mean(mean_dens, na.rm = T),
              climvar = mean(climvar, na.rm = T)) %>%
    mutate(grid_lat = substr(grid, 1, 5),
           grid_lon = as.numeric(substr(grid, 6, 11))) %>%
  ggplot() +
    geom_tile(aes(x = -grid_lon, y = grid_lat, fill = mean_dens)) +
  geom_text(aes(x = -grid_lon, y = grid_lat, label = patch), color = "white") +
  viridis::scale_fill_viridis()
  # ggsave(here(plotsave, "AA_grid_map.png"), height = 6, width = 8, dpi = 400)

# plot to see which patches were selected (excluding small scallops < 80 mm)
dat_dens80 %>%
  group_by(grid, patch) %>%
  summarise(mean_dens = mean(mean_dens, na.rm = T),
            climvar = mean(climvar, na.rm = T)) %>%
  mutate(grid_lat = substr(grid, 1, 5),
         grid_lon = as.numeric(substr(grid, 6, 11))) %>%
  ggplot() +
  geom_tile(aes(x = -grid_lon, y = grid_lat, fill = mean_dens)) +
  geom_text(aes(x = -grid_lon, y = grid_lat, label = patch), color = "white") +
  viridis::scale_fill_viridis() +
  labs(title = "Excluding small scallops (< 80 mm)")
# ggsave(here(plotsave, "AA_grid_map_large_only.png"), height = 6, width = 8, dpi = 400)

# map where the NA values are
dat_dens %>% 
  left_join(select(use_patches, grid, dens_NA, sbt_NA), by = "grid") %>%
  group_by(grid, patch) %>%
  summarise(dens_NA = mean(dens_NA, na.rm = T),
            sbt_NA = mean(sbt_NA, na.rm = T)) %>%
  mutate(grid_lat = substr(grid, 1, 5),
         grid_lon = as.numeric(substr(grid, 6, 11))) %>%
  ggplot() +
  geom_tile(aes(x = -grid_lon, y = grid_lat, fill = dens_NA)) +
  geom_text(aes(x = -grid_lon, y = grid_lat, label = dens_NA), color = "white") +
  viridis::scale_fill_viridis() +
  theme_bw() +
  labs(title = season, fill = "NA dens\n(years)", x = "", y = "")
# ggsave(here(plotsave,
#             paste0("AAAA_densNA_years_grid_map_",
#               season, ".png")), height = 6, width = 8, dpi = 400)

# split into train and test data
dat_test_lengths <- dat_lengths %>% 
  filter(year %in% test_years)

dat_train_lengths <- dat_lengths %>% 
  filter(year %in% train_years)

dat_test_dens <- dat_dens %>% 
  filter(year %in% test_years)

dat_train_dens <- dat_dens %>% 
  filter(year %in% train_years)

# versions of dat_test/train_dens without small scallops for plotting
dat_test_dens80 <- dat_dens80 %>% 
  filter(year %in% test_years)

dat_train_dens80 <- dat_dens80 %>% 
  filter(year %in% train_years)

# get time dimension
years <- sort(unique(dat_train_lengths$year)) 
years_proj <- sort(unique(dat_test_lengths$year))
ny <- length(years)
ny_proj <- length(years_proj)

# get patch area 
patchdat <- dat %>% 
  select(grid, grid_lat, grid_lon) %>%
  distinct() %>%
  filter(is.na(grid_lat)==FALSE) %>%
  mutate(max_lon = grid_lon+0.5,
         min_lon = grid_lon-0.5) %>% # these are +0.25 when using 1/2 degree cells 
  rowwise() %>% 
  mutate(lon_dist = distGeo(p1=c(max_lon, grid_lat), p2=c(min_lon, grid_lat))/1000, # get distance between the furthest longitudes in km, at the midpoint of the lat band 
         patch_area_km2 = lon_dist * 111) %>%  # 1 degree latitude = 111 km 
  select(grid, patch_area_km2) %>% 
  filter(grid %in% patches) %>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(grid)))
meanpatcharea <- mean(patchdat$patch_area_km2)

# set fixed parameters from stock assessment
loo.gb <- mean(c(143.9, 140.0, 148.6, 121.1, 144.9, 152.5, 143.6, 143.3, 145.1))
loo.ma <- mean(c(133.3, 151.8))
loo = mean(c(loo.gb, loo.ma)) # mean from Hart & Chute: 142.6
k = 0.37 # means from Hart & Chute; GB: 0.33; MAB: 0.40
m = 0.35 # 2018 assessment: M = 0.35 for all years and ages.
# z = exp(-m-f) # THIS DISAPEARED FROM ORIGINAL CODE; IS IT NO LONGER NEEDED?
age_at_maturity = 4
t0=-.2 # NEED TO UPDATE FOR SCALLOP?
cv= 0.2 # guess
min_age = 1
max_age = 14 # based on age data in scallop dredge survey

if(time_varying_f==TRUE){
# now that we have the max_age, fill in f for years above 7 (since the f for age=7 is really for 7+)
# SCALLOP: don't need to do this as we input data for all ages
# older_ages <- expand_grid(age=seq(max(dat_f_age_prep$age)+1, max_age, 1), year= unique(dat_f_age_prep$year)) %>% 
#   left_join(dat_f_age_prep %>% filter(age==max(age)) %>% select(year, f))
# dat_f_age_prep %<>% bind_rows(older_ages)
}

# make length to age conversions
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

# prep f data
if(time_varying_f==TRUE){
  dat_f_age <- dat_f_age_prep %>% 
    filter(year %in% years) 
  
  dat_f_age_proj <- dat_f_age_prep %>% 
    filter(year %in% years_proj) %>%
    bind_rows(dat_f_age %>% filter(year==max(year))) # need final year of training data to initialize projection
} else {
  f_prep=0.31 # average of 2008-2019 from 2018 scallop stock Assess & 2020 update
}

lbins <- unique(length_at_age_key$length_bin)
# lbins <- sort(unique(dat_train_lengths$length))
n_lbins <- length(lbins) 

n_ages <- nrow(l_at_a_mat)

# selectivity_at_bin
selectivity_at_bin_df <- 
  read.csv(here("processed-data",
                "scallop_trawl_selectivity_at_length.csv")) %>%
  filter(bin %in% as.character(lbins+0.5)) %>%
  select(selectivity)

selectivity_at_bin <- selectivity_at_bin_df$selectivity
rm(selectivity_at_bin_df)

# now that years are defined above, convert them into indices in the datasets
# be sure all these dataframes have exactly the same year range! 
year_index <- data.frame(year = sort(unique(c(dat_train_dens$year, dat_test_dens$year)))) %>%
  mutate(index = as.integer(as.factor(year)))

dat_train_dens <- dat_train_dens %>%
  left_join(year_index) %>% mutate(year = index) %>% select(-index)

dat_train_dens80 <- dat_train_dens80 %>%
  left_join(year_index) %>% mutate(year = index) %>% select(-index)

dat_test_dens <- dat_test_dens %>%
  left_join(year_index) %>% mutate(year = index) %>% select(-index)

dat_test_dens80 <- dat_test_dens80 %>%
  left_join(year_index) %>% mutate(year = index) %>% select(-index)

dat_train_lengths <- dat_train_lengths %>%
  left_join(year_index) %>% mutate(year = index) %>% select(-index)

dat_test_lengths <- dat_test_lengths %>%
  left_join(year_index) %>% mutate(year = index) %>% select(-index)

if(time_varying_f==TRUE){
  dat_f_age$year = as.integer(as.factor(dat_f_age$year))
  dat_f_age_proj$year = as.integer(as.factor(dat_f_age_proj$year))
}

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

plot(len[20,,20])

# create mean density array for Stan
dens <- array(NA, dim=c(np, ny))
for(p in 1:np){
  print(p)
  for(y in 1:ny){
    tmp2 <- dat_train_dens %>% filter(patch==p, year==y) 
    dens[p,y] <- tmp2$mean_dens * meanpatcharea
  }
}

# create climate variable array for Stan (called sbt but could be o2 as well)
sbt <- array(NA, dim=c(np,ny))
for(p in 1:np){
  print(p)
  for(y in 1:ny){
    tmp3 <- dat_train_dens %>% filter(patch==p, year==y) 
    sbt[p,y] <- tmp3$climvar
  }
}

# or lagged sbt (1 yr)
if(lagged_sbt == 1){
sbt <- array(NA, dim=c(np,ny))
for(p in 1:np){
  print(p)
  for(y in 1:ny){
    tmp3 <- dat_train_dens %>% filter(patch==p, year==y) 
    sbt[p,y] <- tmp3$lclimvar
  }
}
} # close if lagged_sbt = 1

# make matrix identifying NA locations within dens or sbt
dens_notNA <- array(1, dim(dens))
dens_notNA[which(is.na(dens))] <- 0
dens_notNA[which(is.na(sbt))] <- 0

# assign a dummy numeric value to all NAs in dens and sbt
dens[is.na(dens)] <- 999999
sbt[is.na(sbt)] <- 999999

sbt_proj <- array(NA, dim=c(np,ny_proj))
for(p in 1:np){
  print(p)
  for(y in 1:ny_proj){
    tmp6 <- dat_test_dens %>% filter(patch==p, year==(y+ny)) 
    if(nrow(tmp6) != 0){
    sbt_proj[p,y] <- tmp6$climvar
    }
  }
}

sbt_proj[is.na(sbt_proj)] <- 999999

f <- array(NA, dim=c(n_ages,ny))
for(a in min_age:max_age){
  for(y in 1:ny){
    if(time_varying_f==TRUE){
      tmp4 <- dat_f_age %>% filter(age==a, year==y) 
      f[a,y] <- tmp4$f # Alexa had a+1 because matrix indexing starts at 1 not 0; but scallops has no age 0
    } else{
      f[a,y] <- f_prep
    }
  }
}

# f_proj - f values for projected years
f_proj <- array(NA, dim=c(n_ages,(ny_proj+1)))
for(a in min_age:max_age){
  for(y in 1:(ny_proj+1)){
    if(time_varying_f==TRUE){
      tmp5 <- dat_f_age_proj %>% filter(age==a, year==y)
      f_proj[a,y] <- tmp5$f # Alexa added 1 because ages started at 0 and matrix indexing starts at 1 not 0; but scallops have no age 0
    } else{
      f_proj[a,y] <-f_prep
    }
  }
}

######
# fit model
######

stan_data <- list(
  np=np,
  n_ages=n_ages,
  ny_train=ny,
  ny_proj=ny_proj,
  n_lbins=n_lbins,
  n_p_l_y = len,
  abund_p_y = dens,
  abund_p_y_notNA = dens_notNA,
  manual_selectivity = manual_selectivity,
  selectivity_at_bin = selectivity_at_bin,
  sbt = sbt,
  sbt_proj=sbt_proj,
  m=m,
  f=f,
  f_proj=f_proj,
  k=k,
  loo=loo,
  t0=t0,
  cv=cv,
  length_50_sel_guess=60, # based on trawl selectivity data from Dvora
  n_lbins = n_lbins, 
  bin_mids=lbins+0.5, # also not sure if this is the right way to calculate the midpoints
  sel_100 = 3, # based on age at length data 
  age_at_maturity = age_at_maturity,
  l_at_a_key = l_at_a_mat,
  wt_at_age = wt_at_age,
  do_dirichlet = do_dirichlet,
  eval_l_comps = eval_l_comps, # evaluate length composition data? 0=no, 1=yes
  T_dep_mortality = T_dep_mortality, 
  T_dep_recruitment = T_dep_recruitment, # think carefully before making more than one of the temperature dependencies true
  spawner_recruit_relationship = spawner_recruit_relationship, 
  run_forecast=run_forecast
)
saveRDS(stan_data, here("processed-data", "scallop_stan_data_20221011a.rds"))
# stan_data <- readRDS(here("processed-data", "scallop_stan_data_20221011a.rds"))

warmups <- 400
total_iterations <- 500
max_treedepth <-  10
n_chains <- 2
n_cores <- 2
n_thin <- 1
np <- stan_data$np
init_rec1 <- dat_train_dens %>% group_by(patch) %>% summarize(mean_dens = mean(mean_dens, na.rm = T))
init_rec <- log(1e-06 + init_rec1$mean_dens*100)

stan_model_fit <- stan(file = here::here("src","process_sdm_based_on_20221003.stan"), # check that it's the right model!
                       data = stan_data,
                       chains = n_chains,
                       warmup = warmups,
                       thin = n_thin,
                           init = list(list(log_mean_recruits = init_rec, #rep(log(1000), np),
                                            theta_d = 1,
                                           ssb0=1000000),
                                       list(log_mean_recruits = init_rec, #rep(log(1000), np),
                                            theta_d = 1,
                                            ssb0=1000000)),
                       iter = total_iterations,
                       cores = n_cores,
                       refresh = 100,
                       save_warmup = F,
                       save_dso = F,
                       control = list(max_treedepth = max_treedepth,
                                      adapt_delta = 0.85)
)

saveRDS(stan_model_fit, here("results","stan_model_fit_run20221011a.rds"))
stan_model_fit <- readRDS(here("results","stan_model_fit_run20221011a.rds"))

bayesplot::mcmc_pairs(stan_model_fit)

# library(shinystan)
launch_shinystan(stan_model_fit)

# examine rhats
summary(stan_model_fit)$summary$Rhat

post <- list(
T_adjust = rstan::extract(stan_model_fit, "T_adjust")$T_adjust,
T_adjust_proj = rstan::extract(stan_model_fit, "T_adjust")$T_adjust,
Topt = rstan::extract(stan_model_fit, "Topt")$Topt,
width = rstan::extract(stan_model_fit, "width")$width,
dens_p_y_hat80 = rstan::extract(stan_model_fit, "dens_p_y_hat80")$dens_p_y_hat80,
proj_dens_p_y_hat80 = rstan::extract(stan_model_fit, "proj_dens_p_y_hat80")$proj_dens_p_y_hat80,
dens_p_y_hat80_lambda = rstan::extract(stan_model_fit, "dens_p_y_hat80_lambda")$dens_p_y_hat80_lambda,
proj_dens_p_y_hat80_lambda = rstan::extract(stan_model_fit, "proj_dens_p_y_hat80_lambda")$proj_dens_p_y_hat80_lambda
)
saveRDS(post, "results/stan_model_posts_run20221006b.rds")

quantile(rstan::extract(stan_model_fit, "Topt")$Topt, c(0.025, 0.5, 0.975))
quantile(rstan::extract(stan_model_fit, "width")$width, c(0.025, 0.5, 0.975))
quantile(rstan::extract(stan_model_fit, "sigma_r")$sigma_r, c(0.025, 0.5, 0.975))
quantile(rstan::extract(stan_model_fit, "sigma_obs")$sigma_obs, c(0.025, 0.5, 0.975))
quantile(rstan::extract(stan_model_fit, "beta_obs")$beta_obs, c(0.025, 0.5, 0.975))
quantile(rstan::extract(stan_model_fit, "theta_d")$theta_d, c(0.025, 0.5, 0.975))
quantile(rstan::extract(stan_model_fit, "p_length_50_sel")$p_length_50_sel, c(0.025, 0.5, 0.975))

# examine T_adjust (modeled and projected)
test=apply(rstan::extract(stan_model_fit, "T_adjust")$T_adjust, c(2,3), median)
plot(test[1,], ylim = c(0,1), type = "l")
for(i in 1:np){
  points(test[i,], ylim = c(0,1), type = "l")
}
test=apply(rstan::extract(stan_model_fit, "T_adjust_proj")$T_adjust_proj, c(2,3), median)
plot(test[i,], ylim = c(0,1), type = "l")
for(i in 1:np){
  points(test[i,], ylim = c(0,1), type = "l")
}
rm(test)

# Plot Topt/width posterior predictive distribution
topt_plot = rstan::extract(stan_model_fit, c("Topt", "width"))
lowersbt <- median(topt_plot$Topt)-3*median(topt_plot$width)
uppersbt <- median(topt_plot$Topt)+3*median(topt_plot$width)

topt_plot2 <- data.frame(
  t_adjust = exp(-0.5 * ((seq(lowersbt,uppersbt, by = .5) - median(topt_plot$Topt))/median(topt_plot$width))^2),
  temp = seq(lowersbt,uppersbt, by = .5))

topt_plot3_CI = lapply(1:250, function(x){
  data.frame(temp = seq(lowersbt,uppersbt, by = .5), 
             topt = topt_plot$Topt[x], 
             width = topt_plot$width[x],
             iter = x,
             t_adjust = exp(-0.5 * ((seq(lowersbt,uppersbt, by = .5) - topt_plot$Topt[x])/topt_plot$width[x])^2))}) %>%
  do.call(rbind, .)

ggplot(topt_plot2) +
  geom_line(aes(x = temp, y = t_adjust, group = as.factor(iter)), 
            data = topt_plot3_CI,
            color = "firebrick",
            alpha = .1) +
  geom_line(aes(x = temp, y = t_adjust),
            size = 1) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Maximum montly sea bottom temperature (C)", y = "Relative adult survival") # "Recruitment suitability"
ggsave(here(plotsave, "AAA_plot_Topt&width_curve.png"))
rm(topt_plot, topt_plot2, topt_plot3_CI)

# a = rstan::extract(stan_model_fit, "theta_d")
# hist(a$sigma_obs)
# rstanarm::launch_shinystan(stan_model_fit)


# plot important parameters 
plot(stan_model_fit, pars=c('sigma_r','sigma_obs','d','width','Topt','alpha','beta_obs','theta_d'))
ggsave(here(plotsave, "A_plot_important_parameters1.png"))
plot(stan_model_fit, pars=c('sigma_r','sigma_obs','d','alpha','beta_obs','theta_d'))
ggsave(here(plotsave, "B_plot_important_parameters2.png"))

# assess abundance fits

hist(extract(stan_model_fit, "mean_recruits")$mean_recruits)
# save it to appropriate results folder: C_mean_recruits_hist2.png

abund_p_y <- dat_train_dens %>%
  mutate(abundance = mean_dens * meanpatcharea) %>%
  mutate(grid_lat = substr(grid, 1, 4),
         grid_lon = as.numeric(substr(grid, 5, 9))) # note: this is set up for 1x1 grid currently

abund_p_y80 <- dat_train_dens80 %>%
  mutate(abundance = mean_dens * meanpatcharea) %>%
  mutate(grid_lat = substr(grid, 1, 4),
         grid_lon = as.numeric(substr(grid, 5, 9))) # note: this is set up for 1x1 grid currently
  
# grab subset of cells from MA_9 group
# unique(abund_p_y$patch[abund_p_y$grid %in% MA_9$grid]) # c(25, 26, 30, 31, 32, 43, 44, 45)

abund_p_y_hat <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat[patch,year]) %>%
  left_join(distinct(select(abund_p_y, grid, patch)))

abund_p_y_hat80 <- tidybayes::spread_draws(stan_model_fit, dens_p_y_hat80[patch,year]) %>%
  left_join(distinct(select(abund_p_y, grid, patch)))

### Plot log mean abundance (all sizes) vs. time (observed and predicted)
i = 18
# look at raw density data for the patch if you want
# test <- abund_p_y %>% 
#   filter(grid %in% unique(abund_p_y_hat$grid)[i])
for( i in 1:np ) {
  print(i)
patch_filter <- unique(abund_p_y_hat$grid)[i]
abundance_v_time <- abund_p_y_hat %>% 
  filter(grid %in% patch_filter) %>%
  ggplot(aes(year, dens_p_y_hat+1)) + 
  stat_lineribbon() + 
  geom_point(data = filter(abund_p_y, grid %in% patch_filter), aes(year, abundance+1), color = "blue") +
  facet_wrap(~grid, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer() +
scale_y_log10()
abundance_v_time

ggsave(plot = abundance_v_time, filename=here(plotsave,
                                       paste0("F_log_density_v_time_no_length_comps_patch_",
                                       as.character(i), ".png")), width=5, height=5, dpi = 200)
}

### Plot log mean abundance (>= 80 mm) vs. time (observed and predicted)
i = 18
# look at raw density data for the patch if you want
# test <- abund_p_y %>% 
#   filter(grid %in% unique(abund_p_y_hat$grid)[i])
for( i in 1:np ) {
  print(i)
  patch_filter <- unique(abund_p_y_hat80$grid)[i]
  abundance_v_time80 <- abund_p_y_hat80 %>% 
    filter(grid %in% patch_filter) %>%
    ggplot(aes(year, dens_p_y_hat80+1)) + 
    stat_lineribbon() + 
    geom_point(data = filter(abund_p_y80, grid %in% patch_filter), aes(year, abundance+1), color = "blue") +
    facet_wrap(~grid, scales = "free_y") +
    labs(x="Year",y="Abundance") + 
    scale_fill_brewer() +
    scale_y_log10()
  abundance_v_time80
  
  ggsave(plot = abundance_v_time, filename=here(plotsave,
                                                paste0("F_log_density80_v_time_no_length_comps_patch_",
                                                       as.character(i), ".png")), width=5, height=5, dpi = 200)
}

# time series of predicted densities
abund_p_y_hat %>%
  group_by(year, grid, patch) %>%
  summarise(dens_p_y_hat = mean(dens_p_y_hat)) %>%
  ggplot() + 
  geom_line(aes(x = year, y = dens_p_y_hat, color = grid))

ggsave(filename=here(plotsave,
                     "G_dens_p_y_hat_timeseries.png"), 
       width=5, height=5, dpi = 200)

# time series of predicted densities (excluding < 80 mm scallops)
abund_p_y_hat80 %>%
  group_by(year, grid, patch) %>%
  summarise(dens_p_y_hat80 = mean(dens_p_y_hat80)) %>%
  ggplot() + 
  geom_line(aes(x = year, y = dens_p_y_hat80, color = grid))

ggsave(filename=here(plotsave,
                     "G_dens_p_y_hat80_timeseries.png"), 
       width=5, height=5, dpi = 200)

# time series of measured densities
abund_p_y %>%
  # filter(year < 35) %>%
  ggplot() + 
  geom_line(aes(x = year, y = abundance, color = grid))

ggsave(filename=here(plotsave,
                     "G2_dens_p_y_timeseries.png"), 
       width=5, height=5, dpi = 200)

# time series of measured densities excluding < 80 mm scallops
abund_p_y80 %>%
  # filter(year < 35) %>%
  ggplot() + 
  geom_line(aes(x = year, y = abundance, color = grid))

ggsave(filename=here(plotsave,
                     "G2_dens_p_y80_timeseries.png"), 
       width=5, height=5, dpi = 200)

# MAP predicted densities
i = 7
for( i in 1:ny ) {
  print(i)
  year_filter <- unique(abund_p_y_hat$year)[i]
abund_p_y_hat %>% 
    filter(year %in% year_filter) %>%
    group_by(grid) %>%
    summarize(dens_p_y_hat = mean(dens_p_y_hat),
              .groups = "drop") %>%
  mutate(grid_lat = substr(grid, 1, 4),
         grid_lon = as.numeric(substr(grid, 5, 9))) %>%
    ggplot() + 
    geom_tile(aes(x = grid_lon, y = grid_lat, fill = dens_p_y_hat)) + 
    geom_point(data = filter(abund_p_y, year %in% year_filter), 
               aes(x = grid_lon, y = grid_lat, color = abundance), 
               pch = 16) +
    labs(x="",y="", title = i) +
  viridis::scale_fill_viridis() +
  viridis::scale_color_viridis()

  ggsave(filename=here(plotsave,
                                                paste0("F2_dens_p_hat_map_year_",
                                                       as.character(i), ".png")), width=5, height=5, dpi = 200)
}

# MAP predicted densities of scallops at least 80 mm long
i = 7
for( i in 1:ny ) {
  print(i)
  year_filter <- unique(abund_p_y_hat80$year)[i]
  abund_p_y_hat80 %>% 
    filter(year %in% year_filter) %>%
    group_by(grid) %>%
    summarize(dens_p_y_hat80 = mean(dens_p_y_hat80),
              .groups = "drop") %>%
    mutate(grid_lat = substr(grid, 1, 4),
           grid_lon = as.numeric(substr(grid, 5, 9))) %>%
    ggplot() + 
    geom_tile(aes(x = grid_lon, y = grid_lat, fill = dens_p_y_hat80)) + 
    geom_point(data = filter(abund_p_y80, year %in% year_filter), 
               aes(x = grid_lon, y = grid_lat, color = abundance), 
               pch = 16) +
    labs(x="",y="", title = i) +
    viridis::scale_fill_viridis() +
    viridis::scale_color_viridis()
  
  ggsave(filename=here(plotsave,
                       paste0("F2_dens_p_hat80_map_year_",
                              as.character(i), ".png")), width=5, height=5, dpi = 200)
}

# observed vs. predicted mean abundance
abund_obs_pred <- abund_p_y_hat %>%
  group_by(grid, year) %>%
  summarise(dens_p_y_hat = mean(dens_p_y_hat),
            .groups = "drop") %>%
  left_join(abund_p_y, by = c("grid", "year"))
corr <- cor.test(x=abund_obs_pred$abundance, y=abund_obs_pred$dens_p_y_hat, method = 'spearman')

abund_obs_pred %>%
  ggplot() +
  geom_point(aes(x = abundance, y = dens_p_y_hat),
             color = "firebrick", size = 3,
             shape = 1, alpha = 0.75) +
  annotate("text", x = max(abund_obs_pred$abundance,na.rm = T)-100, 
           y = min(abund_obs_pred$dens_p_y_hat)+100, 
           label = paste0("SRC = ", round(corr$estimate,3)), hjust = 1, vjust = 0) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Observed abundance", y = "Predicted abundance")

ggsave(paste0(plotsave, "/H_observed_vs_predicted.png"))

# observed vs. predicted mean abundance excluding small (< 80 mm) scallops
abund_obs_pred80 <- abund_p_y_hat80 %>%
  group_by(grid, year) %>%
  summarise(dens_p_y_hat80 = mean(dens_p_y_hat80),
            .groups = "drop") %>%
  left_join(abund_p_y80, by = c("grid", "year"))
corr <- cor.test(x=abund_obs_pred80$abundance, y=abund_obs_pred80$dens_p_y_hat80, method = 'spearman')

abund_obs_pred80 %>%
  ggplot() +
  geom_point(aes(x = abundance, y = dens_p_y_hat80),
             color = "firebrick", size = 3,
             shape = 1, alpha = 0.75) +
  annotate("text", x = max(abund_obs_pred80$abundance,na.rm = T)-100, 
           y = min(abund_obs_pred80$dens_p_y_hat80)+100, 
           label = paste0("SRC = ", round(corr$estimate,3)), hjust = 1, vjust = 0) +
  # geom_abline(slope = 1, intercept = 0, color = "darkgray", lty = 2) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  labs(x = "Observed abundance", y = "Predicted abundance",
       title = "Large scallops only (80+ mm)")

ggsave(paste0(plotsave, "/H_observed80_vs_predicted80.png"))

# # assess length comp fits
# n_p_l_y_hat <- tidybayes::gather_draws(stan_model_fit, n_p_l_y_hat[year,patch,length])
# saveRDS(n_p_l_y_hat, "results/n_p_l_y_hat_run20220426a.rds")
# 
# n_p_l_y_hat <- readRDS("results/n_p_l_y_hat_run20220413a.rds")
# 
# # neff <- tidybayes::gather_draws(stan_model_fit, n_eff[patch,year], n = 500)
# 
# #neff <- tidybayes::gather_draws(stan_model_fit, n_eff[patch,year], n = 500)
# 
# dat_train_lengths <- dat_train_lengths %>% 
#   group_by(patch, year) %>% 
#   mutate(p_length = sum_num_at_length / sum(sum_num_at_length))
# 
# dat_train_lengths %>% 
#   group_by(year,patch) %>% 
#   summarise(n = sum(sum_num_at_length)) %>% 
#   # filter(year != 35) %>%
#   ggplot(aes(year, n, color =factor(patch))) + 
#   geom_point()
# ggsave(here(plotsave, "J_n_by_year_and_patch.png"))
# 
# p = 15
# 
# n_p_l_y_hat %>% 
#   ungroup() %>% 
#   filter(patch == p) %>% 
#   group_by(patch, year, .iteration) %>% 
#   mutate(pvalue = .value / sum(.value)) %>% 
#   ggplot(aes(length, pvalue)) + 
#   stat_lineribbon() + 
#   geom_point(data = dat_train_lengths %>% filter(patch == p), aes(length*10,p_length), color = "red", alpha = 0.2) +
#   facet_wrap(~year, scales = "free_y")
# 
# ggsave(here(plotsave, "K_p_vs_length_by_year_patch15_times10.png"), width = 10, height= 10, dpi =400)
# 
# # length frequency over time
# l_freq_time <- n_p_l_y_hat %>% 
#   ungroup() %>% 
#   ggplot(aes(x=length, y=..density.., weight=.value)) + 
#   geom_histogram(bins=50) +
#   facet_grid(patch~year)
# 
# ggsave(l_freq_time, filename=here(plotsave,"L_length_freq_time_scallop.png"), scale=1.5, width=15, height=10, dpi = 400)
# 
# # is there a temperature - recruitment relationship? 
# 
# ## first need to figure out which lengths correspond to age 0 
# length_at_age_key %>% 
#   ggplot(aes(x=length_bin, y=p_bin)) +
#   geom_line() + 
#   facet_wrap(~age)
# ggsave(here(plotsave, "M_length_at_age_curves.png"))
# 
# ## let's call recruits <5cm
# 
# # selectivity_at_bin_df <- gather_draws(stan_model_fit, selectivity_at_bin[n_lbins])
# 
# selectivity_at_bin_df <- data.frame(n_lbins = lbins,
#                                     .value = selectivity_at_bin)
# 
# mean_sel <- selectivity_at_bin_df %>% 
#   group_by(n_lbins) %>% 
#   summarise(mean_sel = mean(.value)) %>% 
#   rename(length = n_lbins)
# 
# # get proportion of recruits adjusted by selectivity
# recruits <- n_p_l_y_hat %>% 
#   left_join(mean_sel, by="length") %>% 
#   mutate(recruit = ifelse(length <=5, "yes", "no"),
#          val_post_sel = .value / mean_sel) %>% 
#   group_by(recruit, patch, year, .draw) %>% 
#   summarise(sumcount = sum(val_post_sel, na.rm = T)) %>% 
#   ungroup() %>% 
#   group_by(patch, year, .draw) %>% 
#   mutate(prop_recruit = sumcount / sum(sumcount, na.rm = T)) %>% 
#   filter(recruit == "yes")
# 
# dat_train_sbt <- dat_train_dens # SCALLOP: because I didn't define a stand-alone sbt data frame
# recruits %>% group_by(patch, year) %>% 
#   summarize(mean_prop_rec = mean(prop_recruit)) %>% 
#   left_join(dat_train_sbt, by=c('patch','year')) %>% 
#   ggplot(aes(x=mean_sbt, y=mean_prop_rec, color=year)) +
#   geom_point() +
#   viridis::scale_color_viridis()
# 
# ggsave(here(plotsave, "N_recruits.png"), height = 6, width = 6, dpi = 400)
# 
# # another way to look at recruits vs. sbt
# recruits %>% group_by(patch, year) %>% 
#   summarize(mean_prop_rec = mean(prop_recruit)) %>% 
#   left_join(dat_train_sbt, by=c('patch','year')) %>%
#   filter(is.na(mean_sbt)==F) %>%
#   mutate(sbt_cat = cut(mean_sbt, breaks = c(0,6,9,12,15,18,21,24,30))) %>%
#   ggplot() +
#   geom_boxplot(aes(x = sbt_cat, y = mean_prop_rec)) +
#   theme_bw() +
#   labs(x = "Sea bottom temperature (C)")
# 
# ggsave(here(plotsave, "N2_recruits_boxplot.png"), height = 6, width = 6, dpi = 400)
# 
# # plot recruitment deviates
# # note that the length of rec_dev is actually 34 not 35
# gather_draws(stan_model_fit, rec_dev[ny]) %>% 
#   ggplot(aes(x=ny, y=.value)) + 
#   stat_lineribbon() + 
#   scale_fill_brewer()
# 
# ggsave(here(plotsave, "O_rec_dev.png"), width = 6, height = 6, dpi = 400)
# 
# # plot raw
# raw_v_time <- gather_draws(stan_model_fit, raw[ny]) %>%
#   ggplot(aes(x=ny, y=.value)) + 
#   stat_lineribbon() + 
#   scale_fill_brewer()
# ggsave(raw_v_time, filename=here(plotsave,"P_raw_v_time.png"))
# 
# # plot actual recruitment
# n_p_a_y_hat <- gather_draws(stan_model_fit, n_p_a_y_hat[np, n_ages, ny]) %>%
#   left_join(distinct(dat_dens[,2:3]), by = c("np" = "patch"))
# 
# n_p_a_y_hat_sum <- n_p_a_y_hat %>%
#   group_by(grid, np, n_ages, ny) %>%
#   summarise(nhat = median(.value),
#             q2.5 = quantile(.value, 0.025),
#             q97.5 = quantile(.value, 0.975),
#             .groups = "drop")
# 
# i = 15
# for( i in unique(n_p_a_y_hat_sum$np)) {
#   print(i)
#   patch_filter <- unique(n_p_a_y_hat_sum$grid)[i]
# n_p_a_y_hat_sum %>%
#   filter(grid == patch_filter) %>%
#   ggplot(aes(x=ny, y=nhat)) +
#   geom_line(color = "firebrick") +
#   geom_ribbon(aes(x = ny, ymin = q2.5, ymax = q97.5),
#               fill = "firebrick", alpha = 0.25) +
#   facet_wrap(~n_ages) +
#   labs(y = "Individuals", title = patch_filter) + 
#   theme_bw()
# #scale_y_log10()
# ggsave(here(plotsave, 
#             paste0("Q_recruitment_n_age_by_patch_year",
#             "patch_", patch_filter, ".png")),
#        height = 10, width = 10, dpi = 400)
# }
# 
# # plot all ages over time
# n_p_a_y_hat %>% 
#   ggplot(aes(x=ny, y=.value)) +
#   stat_lineribbon(size = .5) +
#   scale_fill_brewer() +
#   facet_grid(np~n_ages)
# ggsave(here(plotsave, "R_n_p_a_y_hat_vs_age_year.png"),
#        width = 15, height = 15, dpi = 400)
# 
# # detection stats 
# detect <- spread_draws(stan_model_fit, theta[patch,year])
# 
# i = 15
# for( i in unique(n_p_a_y_hat_sum$np)) {
#   print(i)
# detect %>%
#    filter(patch == i) %>% 
#   ggplot(aes(x=year, y=theta)) +
#   stat_lineribbon() + 
#   facet_wrap(~patch) +
#   scale_fill_brewer()
# 
# ggsave(here(plotsave, 
#             paste0("S_detection_stats_patch_", 
#                    i, ".png")))
# }
# 
beta_obs <- extract(stan_model_fit, "beta_obs")$beta_obs
hist(beta_obs)

median(beta_obs)

theta_eq <- data.frame(dens = 1:45000,
                       theta = ((1/(1+exp(-median(beta_obs)*1:45000))) - 0.5)*2)
theta_eq %>%
  ggplot() +
  geom_line(aes(x = dens, y = theta, group = 1)) +
  theme_bw()
ggsave(here(plotsave,
            "S2_beta_obs_effect_graphed.png"))

# patch abundance observed vs. expected tile plots
est_patch_abund <- abund_p_y_hat %>% 
  group_by(year, patch) %>% 
  summarise(abundance = mean(dens_p_y_hat))

est_patch_abund80 <- abund_p_y_hat80 %>% 
  group_by(year, patch) %>% 
  summarise(abundance = mean(dens_p_y_hat80))

observed_abundance_tile <- abund_p_y %>% 
  # filter(abundance < 20000000) %>%
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 36, 4)) +
  scale_y_continuous(breaks=seq(0, 36, 4)) +
  labs(title="Observed", x="Year", y="Patch", fill="Abundance") +
  viridis::scale_fill_viridis()

observed_abundance_tile80 <- abund_p_y80 %>% 
  # filter(abundance < 20000000) %>%
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 36, 4)) +
  scale_y_continuous(breaks=seq(0, 36, 4)) +
  labs(title="Observed", x="Year", y="Patch", fill="Abundance") +
  viridis::scale_fill_viridis()

estimated_abundance_tile <- est_patch_abund %>% 
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 36, 4)) +
  scale_y_continuous(breaks=seq(0, 36, 4)) +
  labs(title="Estimated", x="Year", y="Patch", fill="Abundance") +
  viridis::scale_fill_viridis()

estimated_abundance_tile80 <- est_patch_abund80 %>% 
  ggplot(aes(x=year, y=patch, fill=abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(0, 36, 4)) +
  scale_y_continuous(breaks=seq(0, 36, 4)) +
  labs(title="Estimated", x="Year", y="Patch", fill="Abundance") +
  viridis::scale_fill_viridis()

ggsave(observed_abundance_tile, filename=here(plotsave,"U_abundance_v_time_observed_tileplot_no_length_comps2.png"))
ggsave(estimated_abundance_tile, filename=here(plotsave,"V_abundance_v_time_estimated_tileplot_no_length_comps.png"))
ggsave(observed_abundance_tile80, filename=here(plotsave,"U_abundance80_v_time_observed_tileplot_no_length_comps2.png"))
ggsave(estimated_abundance_tile80, filename=here(plotsave,"V_abundance80_v_time_estimated_tileplot_no_length_comps.png"))

# plot temperature difference from optimum over space and time
Topt <- extract(stan_model_fit, "Topt")$Topt
hist(Topt) # note, not normally distributed

dat_train_sbt <- dat_train_dens # SCALLOP: because I didn't define a stand-alone sbt data frame
dat_train_sbt %>%
  mutate(Tdiff = mean_sbt - median(Topt)) %>%
  ggplot(aes(x=year, y=patch, fill=Tdiff)) +
  geom_tile() + 
  scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=0)

ggsave(here(plotsave, "W_Tdiff_by_patch_yearB.png"))





##################################
## ... not using the code below yet



########
# evaluate forecast
########
proj_dens_p_y_hat <- gather_draws(stan_model_fit, proj_dens_p_y_hat[np, ny_proj])

proj_abund_p_y <- dat_test_dens %>%
  mutate(abundance = mean_dens * meanpatcharea)
# left_join(patchdat, by = c("lat_floor", "patch")) %>% 
# group_by(patch, year) %>% 
# summarise(abundance = sum(mean_dens *patch_area_km2)) %>% 
#  ungroup()

proj_observed_abundance_tile <- proj_abund_p_y %>% 
  mutate(Year = (year + min(years_proj) - 1), Latitude = (patch + min(patches) - 1), Abundance=abundance) %>% 
  ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Observed")


proj_est_patch_abund <- proj_dens_p_y_hat %>% 
  group_by(ny_proj, np) %>% 
  summarise(abundance = mean(.value))

proj_estimated_abundance_tile <- proj_est_patch_abund %>% 
  mutate(Year = (ny_proj + min(years_proj) - 1), Latitude = (np + min(patches) - 1), Abundance=abundance) %>% 
  ggplot(aes(x=Year, y=Latitude, fill=Abundance)) +
  geom_tile() +
  theme_bw() +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1)) +
  scale_fill_continuous(labels = scales::comma) + # fine to comment this out if you don't have the package installed, it just makes the legend pretty
  labs(title="Estimated")
ggsave(proj_estimated_abundance_tile, filename=here("results","proj_estimated_abundance_v_time_tileplot.png"), scale=0.9)
ggsave(proj_observed_abundance_tile, filename=here("results","proj_observed_abundance_v_time_tileplot.png"), scale=0.9)

proj_abundance_v_time <- proj_dens_p_y_hat %>% 
  rename(patch=np) %>% 
  ggplot(aes(ny_proj, .value)) + 
  stat_lineribbon() + 
  geom_point(data = proj_abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_x_continuous(breaks=seq(0, 10, 2), limits=c(0, 10)) +
  theme(legend.position="none") +
  scale_fill_brewer()
ggsave(proj_abundance_v_time, filename=here("results","proj_density_v_time.png"), width=7, height=4)

# centroid position by year 
dat_centroid_proj <- proj_abund_p_y %>% 
  group_by(year) %>% 
  summarise(centroid_lat = weighted.mean(x=patch, w=abundance)) %>%
  mutate(Year = (year + min(years_proj) - 1), Latitude = (centroid_lat + min(patches) - 1))

# model fit centroid -- should eventually estimate in model for proper SE -- just exploring here
est_centroid_proj <- proj_dens_p_y_hat %>% 
  group_by(ny_proj, .iteration) %>%  # IS THIS SUPPOSED TO BE .ITERATION? CHECK WHEN MODEL IS RUN FOR LONGER 
  summarise(centroid_lat = weighted.mean(x=np, w=.value)) %>% 
  ungroup()

gg_centroid_proj_prep <- dat_centroid_proj %>% 
  ggplot(aes(Year, Latitude)) + 
  geom_point(color = "red") +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(37.6, 39.2, 0.2), limits=c(37.6, 39.2)) +
  theme(legend.position = "none") +
  labs(title="Centroid Position") 

gg_centroid_proj <- est_centroid_proj %>% 
  mutate(Year = (ny_proj + min(years_proj) - 1), Latitude = (centroid_lat + min(patches) - 1)) %>% 
  ggplot(aes(Year, Latitude)) + 
  stat_lineribbon() + 
  scale_fill_brewer() +
  geom_point(data = dat_centroid_proj, aes(Year, Latitude), color = "red") +
  scale_x_continuous(breaks=seq(min(years_proj), max(years_proj), 1)) +
  scale_y_continuous(breaks=seq(37.6, 39.2, 0.2)) +
  theme(legend.position = "none") +
  labs(title="Centroid Position") 
ggsave(gg_centroid_proj_prep, filename=here("results","proj_centroid_v_time_prep.png"), scale=0.9)
ggsave(gg_centroid_proj, filename=here("results","proj_centroid_v_time.png"), scale=0.9)
