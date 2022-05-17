library(dplyr)

sf <- read.csv("processed-data/scallop_digitized_average_F.csv")
n_years <- nrow(sf)
n_ages <- 14

year_list <- list()
age_list <- list()
f_list <- list()
for(i in 1:n_years){

year_list[[i]] <- rep(i,n_ages)
age_list[[i]] <- 1:n_ages  
f_list[[i]] <- rep(sf[i,2], n_ages)
}

sf_final <- data.frame(
  Year = do.call(c, year_list),
  Age = do.call(c, age_list),
  F = do.call(c, f_list)
)

# attempt at correcting for fishing selectivity based on Hart&Shank
sf_final_cor <- sf_final %>%
  mutate(F = case_when(Age == 1 ~ F*0.1,
                       Age == 2 ~ F*0.5,
                       TRUE ~ F))

write.csv(sf_final_cor, "processed-data/scallop_F_by_age.csv", 
          row.names = F)  
