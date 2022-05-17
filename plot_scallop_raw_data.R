rbind(dat_train_dens40_plus, dat_test_dens40_plus) %>%
  mutate(type = c(rep("train", nrow(dat_train_dens)), 
                  rep("test", nrow(dat_test_dens)))) %>%
  ggplot() +
  geom_line(aes(x = year, y = mean_dens, group = grid, color = type)) +
  ggtitle("40mm and above")
ggsave("C:\\Users\\Mike\\mike_files\\Research\\scallops\\figures\\drm_explore\\9patch_40plus.png",
       height = 6, width = 6, dpi = 400)

rbind(dat_train_dens_below40, dat_test_dens_below40) %>%
  mutate(type = c(rep("train", nrow(dat_train_dens)), 
                  rep("test", nrow(dat_test_dens)))) %>%
  ggplot() +
  geom_line(aes(x = year, y = mean_dens, group = grid, color = type)) +
  ggtitle("below 40mm")
ggsave("C:\\Users\\Mike\\mike_files\\Research\\scallops\\figures\\drm_explore\\9patch_below40.png",
       height = 6, width = 6, dpi = 400)