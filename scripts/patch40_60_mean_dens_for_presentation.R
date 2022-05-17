patch_filter <- unique(abund_p_y_hat$grid)[40:60]
abund_p_y %>% 
  filter(grid %in% patch_filter) %>%
  group_by(year) %>%
  summarise(mean_dens = mean(mean_dens, na.rm = T), 
            .groups = "drop") %>%
  ggplot() + 
  geom_line(aes(year, mean_dens)) +
  geom_point(aes(year, mean_dens), size = 3) +
  geom_smooth(aes(year, mean_dens), color = "red", size = 2) +
  # facet_wrap(~grid, scales = "free_y") +
  labs(x="Year",y="Abundance") + 
  scale_fill_brewer() +
  scale_y_log10() +
  theme_bw() +
  theme(text = element_text(size = 14))

ggsave(here(plotsave, 
            paste0("FORPRES_patches40-60_mean_density_v_time.png")), 
       width=5, height=5, dpi = 400)
