# latitude vs. temp.
plot(hd_df_step1$lat, hd_df_step1$mean_bottemp, col = "darkgray",
     xlab = "lattitude", ylab = "mean bottom temp")
points(hd_df_step1[hd_df_step1$mean_tot_num>5,]$lat, 
       hd_df_step1[hd_df_step1$mean_tot_num>5,]$mean_bottemp, 
       col = "orange", pch = 16)
points(hd_df_step1[hd_df_step1$mean_tot_num>500,]$lat, 
       hd_df_step1[hd_df_step1$mean_tot_num>500,]$mean_bottemp, 
       col = "firebrick", pch = 16)
