# check DRM data has all haulids that are in OceanAdapt data
OA_fall <- readRDS(here("data_noaa/trawl/",
                          "dat_explodedNortheast US Fall.rds")) %>%
  filter(substr(haulid, 1, 4) %in% 1963:2015)

OA_spring <- readRDS(here("data_noaa/trawl/",
                        "dat_explodedNortheast US Spring.rds")) %>%
  filter(substr(haulid, 1, 4) %in% 1963:2015)

length(unique(OA_fall$haulid)) # 17923
length(unique(OA_spring$haulid)) # 16329
# sum = 34252 or 477 more than the DRM data. Need to check.
length(unique(hauls_drm$haulid)) # 33775


show increasing proportion of hauls w/ scallops 1963-2015
test = hauls_drm %>% group_by(haulid, year) %>% summarize(n = sum(number_at_length), .groups = "drop")
test2 = table(test$year, test$n>0) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rename(year = 1, pres = 2, Freq = 3) %>%
  pivot_wider(id_cols = c(year), names_from = pres,
              values_from = Freq) %>%
  rename(year = 1, not_pres = 2, pres = 3) %>%
  mutate(p = pres/(pres + not_pres))