# drake strip comparison

source("R/drake_data_prep.R")

long_d <- xdata %>%
  pivot_longer(cols = dplyr::select_if(.,is.numeric) %>% names(),
               names_to = "variable",
               values_to = "value")

p<-long_d %>%
  ggplot(aes(x = strip_type, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable, scales="free")
  

ggsave(p, filename = "figs/big_facet.png",
       width=20,height=20, bg="white")

all_cover <- plant_cover %>%
  group_by(strip_type, plot, species_code) %>%
  summarise(cover_pct = sum(cover_pct, na.rm=T)/n_subplots) %>%
  ungroup()

occurrences <- plant_cover %>%
  group_by(plot, species_code) %>%
  summarise(prevalence = n()/n_subplots) %>%
  ungroup() %>%
  unique() %>%
  pivot_wider(id_cols = "plot",
              names_from = "species_code",
              values_from = "prevalence",
              values_fill = 0) %>%
  pivot_longer(cols = names(.)[2:ncol(.)],
               names_to = "species_code",
               values_to = "prevalence") %>%
  mutate(strip_type = str_sub(plot, 7,10)) %>%
  group_by(strip_type, species_code) %>%
  summarise(prevalence = sum(prevalence)) %>%
  ungroup()
  

p1 <- occurrences %>%
  ggplot(aes(x = strip_type, y = prevalence, fill=strip_type)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~species_code)

ggsave(p1, filename = "figs/plants_facet.png",
       width=10,height=10, bg="white")
