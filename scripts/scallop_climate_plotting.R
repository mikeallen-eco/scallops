# latitude vs. temp.
plot(hd_df$lat, hd_df$mean_bottemp, col = "darkgray",
     xlab = "lattitude", ylab = "mean bottom temp")
points(hd_df[hd_df$mean_tot_num>5,]$lat, 
       hd_df[hd_df$mean_tot_num>5,]$mean_bottemp, 
       col = "orange", pch = 16)
points(hd_df[hd_df$mean_tot_num>500,]$lat, 
       hd_df[hd_df$mean_tot_num>500,]$mean_bottemp, 
       col = "firebrick", pch = 16)

hd_df %>%
  ggplot(aes(x = lat, y = mean_bottemp)) +
  geom_point(color = "darkgray", pch = 1) +
  geom_point(data = filter(hd_df, mean_tot_num>5), 
             color = "orange") +
  geom_point(data = filter(hd_df, mean_tot_num>500), 
             color = "firebrick") +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  labs(x = "Latitude", y = "Mean bottom temp.")

ggsave("figures/scallop_bottom_temps.png", dpi = 400, height = 4, width = 6)

library(gganimate)
library(gifski)
test <- hd_df %>%
    ggplot(aes(x = lat, y = mean_bottemp)) +
    geom_point(color = "darkgray", pch = 1) +
    geom_point(data = filter(hd_df, mean_tot_num>5), 
               color = "orange") +
    geom_point(data = filter(hd_df, mean_tot_num>500), 
               color = "firebrick") +
    theme_bw() +
    theme(text = element_text(size = 15)) +
    labs(x = "Latitude", y = "Mean bottom temp.", 
         title = 'Year: {frame_time}') +
  transition_time(year)# +
  # enter_fade() + 
  # exit_shrink() +
  # ease_aes('sine-in-out') 

  animate(test, gifski_renderer(file = "test.gif"), nframes = 47, start_pause = 1)
