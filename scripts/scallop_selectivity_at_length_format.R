library(here)
library(dplyr)

# Scallops
# read in trawl selectivity data in 5 cm length bins from Dvora
sel <- 
  read.csv(here("data_noaa/trawl/", 
      "scallop_trawl_selectivity_at_length_original_bins.csv"))

# make a 1 cm bin midpoint vector of all possible lengths
mids <- seq(1.5, 249.5, by = 1)

# make a column of 5 cm bin mids to match up with the 1 cm bins
bin_lookup <- c(rep(sel$bin[1], 4),
  do.call(c, lapply(sel$bin[2:28], function(x)rep(x, 5))),
  rep(sel$bin[28], 110)
)

# make a dataframe combining both
sel_length <- data.frame(
  mids, bin_lookup) %>%
  left_join(sel, by = c("bin_lookup" = "bin")) %>%
  select(bin = mids, selectivity = mean_sel)

# write to csv file
write.csv(sel_length, here("data_noaa/trawl/", 
                "scallop_trawl_selectivity_at_length.csv"),
          row.names = F)
