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
  ungroup() %>%
  unique()

occurrences <- plant_cover %>%
  group_by(plot, species_code) %>%
  summarise(prevalence = n()) %>%
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
  ungroup() %>%
  group_by(species_code) %>%
  mutate(totalprev = sum(prevalence)) %>%
  ungroup() %>%
  filter(totalprev>0)
  
sp_list <- read_csv("data/drake_veg_data_2022 - species_list.csv") %>%
  dplyr::select(species_code = code,group, introduced, perennial, woody, graminoid,
                rhizomatous, pp = photosynthetic_pathway) %>%
  mutate(fg = paste0(ifelse(introduced=="yes", "I", "N"),
                     ifelse(perennial == "yes", "P", "A"),
                     ifelse(graminoid == "yes", "G", "D"))) %>%
  dplyr::select(species_code, fg)

p1 <- occurrences %>%
  left_join(sp_list) %>%
  filter(str_sub(fg, 1,1) == "N") %>%
  ggplot(aes(x = strip_type, y = prevalence, fill=strip_type)) +
  geom_bar(stat = 'identity') +
  facet_wrap(fg~species_code, ncol=6) +
  theme_classic();p1

p2 <- occurrences %>%
  left_join(sp_list) %>%
  filter(str_sub(fg, 1,1) == "I") %>%
  ggplot(aes(x = strip_type, y = prevalence, fill=strip_type)) +
  geom_bar(stat = 'identity') +
  facet_wrap(fg~species_code, ncol=5) +
  theme_classic();p2

ggsave(p1, filename = "figs/plants_facet.png",
       width=10,height=10, bg="white")
